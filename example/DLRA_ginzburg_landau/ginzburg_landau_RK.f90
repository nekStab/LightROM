module Ginzburg_Landau_RK
   ! Standard Library.
   use stdlib_optval, only : optval
   use stdlib_linalg, only : diag, svd, svdvals
   use stdlib_io_npy, only : save_npy, load_npy
   ! LightKrylov for linear algebra.
   use LightKrylov
   use LightKrylov, only : wp => dp
   use LightKrylov_AbstractVectors ! linear_combination
   use LightKrylov_Utils ! svd, sqrtm
   ! LightROM
   use LightROM_Utils ! Balancing_Transformation
   ! Lyapunov Solver
   use LightROM_LyapunovSolvers
   use LightROM_LyapunovUtils
   ! Riccati Solver
   use LightROM_RiccatiSolvers
   use LightROM_RiccatiUtils
   ! Ginzburg Landau
   use Ginzburg_Landau_Base
   use Ginzburg_Landau_Operators
   use Ginzburg_Landau_RK_Lyapunov
   use Ginzburg_Landau_Utils
   use Ginzburg_Landau_Tests
   !use fortime
   implicit none

   private :: this_module
   public  :: run_lyap_reference_RK

   character*128, parameter :: this_module = 'Ginzburg_Landau_RK'

contains

   subroutine run_lyap_reference_RK(LTI, Xref_BS, Xref_RK, U0, S0, Tend, nsteps, iref, ifload, ifverb)
      ! LTI system
      type(lti_system),              intent(inout) :: LTI
      ! Reference solution (BS)
      real(wp),                      intent(in)    :: Xref_BS(N,N)
      ! Reference solution (RK)
      real(wp),                      intent(out)   :: Xref_RK(N,N)
      ! Initial condition
      type(state_vector),            intent(inout) :: U0(:)
      real(wp),                      intent(inout) :: S0(:,:)
      real(wp),                      intent(in)    :: Tend
      integer,                       intent(in)    :: nsteps
      integer,                       intent(in)    :: iref
      ! Optional
      logical, optional,             intent(in)    :: ifload
      logical                                      :: if_load_npy
      logical, optional,             intent(in)    :: ifverb
      logical                                      :: verb

      ! Internals
      type(rk_lyapunov),             allocatable   :: RK_propagator
      type(state_matrix)                           :: X_mat(2)
      real(wp),                      allocatable   :: X_RKlib(:,:,:)
      real(wp)                                     :: Bmat(N,2)
      real(wp)                                     :: U0_mat(N,N)
      integer                                      :: irep, printstep
      real(wp)                                     :: etime, Tstep
      ! OUTPUT
      real(wp)                                     :: U_out(N,N)
      real(wp)                                     :: X_out(N,N)
      real(wp),           allocatable              :: U_load(:,:)
      real(wp)                                     :: X_mat_flat(N**2)
      real(wp)                                     :: res_flat(N**2)
      character*128      :: oname, note
      ! timer
      integer            :: clock_rate, clock_start, clock_stop
      logical            :: if_save_npy

      if_load_npy  = optval(ifload, .false.)
      verb         = optval(ifverb, .false.)

      call system_clock(count_rate=clock_rate)

      printstep = 10
      Tstep = Tend/nsteps
      allocate(X_RKlib(N, N, nsteps))
      if (if_load_npy) then
         print *, ''
         print *, 'Load the solution to the differential Lyapunov equation computed with RKlib:'
         print *, ''
      else
         if_save_npy = .true.
         print *, ''
         print *, 'Compute the solution to the differential Lyapunov equation with RKlib:'
         print *, ''
         ! initialize exponential propagator
         RK_propagator = RK_lyapunov(Tstep)
         ! Set random initial condition
         call get_state(U0_mat(:,:rk_X0), U0)
         if (if_save_npy) then
            X_out = matmul( matmul( U0_mat(:,:rk_X0), matmul( S0, transpose(U0_mat(:,:rk_X0)) ) ), weight_mat)
            ! Save forcing RK
            write(oname,'("local/RK/data_GL_lyapconv_BBTW_RK_n",I4.4,".npy")') nx
            call save_data(oname, reshape(BBTW_flat, (/ N,N /)))
            ! Save initial condition
            write(oname,'("local/RK/data_GL_lyapconv_X0_RK_n",I4.4,".npy")') nx
            call save_data(oname, X_out)
            ! Save forcing DLRA
            Bmat = 0.0_wp
            call get_state(Bmat, LTI%B(:rk_b))
            X_out = matmul(Bmat, transpose(Bmat))
            write(oname,'("local/RK/data_GL_lyapconv_BBTW_DLRA_n",I4.4,".npy")') nx
            call save_data(oname, X_out)
         end if
         X_out = matmul(matmul( U0_mat(:,:rk_X0), matmul( S0, transpose(U0_mat(:,:rk_X0)) ) ), weight_mat)
         ! Set initial condition for RK
         call set_state(X_mat(1:1), X_out)
      end if
      write(*,'(A7,A10,A19,A19,A19,A12)') ' RKlib:','Tend','|| X_R ||_2/N', '|| X_R-X_B ||_2/N', '|| res_R ||_2/N','etime'
      write(*,*) '-------------------------------------------------------------------------------------'
      call get_state(U0_mat(:,:rk_X0), U0)
      X_out = matmul(matmul( U0_mat(:,:rk_X0), matmul( S0, transpose(U0_mat(:,:rk_X0)) ) ), weight_mat)
      call CALE(res_flat, reshape(X_out, shape(res_flat)), BBTW_flat, .false.)
      write(*,'(I7,F10.2,3(1X,F18.12),F10.2," s",A)') 0, 0.0, norm2(X_out)/N, &
                        & norm2(X_out - Xref_BS)/N, norm2(res_flat)/N, 0.0, ''
      do irep = 1, nsteps
         write(oname,'("local/RK/data_GL_lyapconv_X_RK_n",I4.4,"_r",I3.3,".npy")') nx, irep
         if (if_load_npy) then
            if (mod(irep,printstep) == 0 .or. irep==iref) then
               call load_data(oname, U_load)
               X_RKlib(:,:,irep) = U_load
            end if
            etime = 0.0_wp
         else
            call system_clock(count=clock_start)     ! Start Timer
            ! integrate
            call RK_propagator%matvec(X_mat(1), X_mat(2))
            call system_clock(count=clock_stop)      ! Stop Timer
            etime = real(clock_stop-clock_start)/real(clock_rate)
            ! recover output
            call get_state(X_RKlib(:,:,irep), X_mat(2:2))
            ! save output
            if (if_save_npy) call save_data(oname, X_RKlib(:,:,irep))
            ! replace input
            call set_state(X_mat(1:1), X_RKlib(:,:,irep))
         end if
         ! compute residual
         if (.not.if_load_npy .or. (mod(irep,printstep) == 0  .or. irep==iref)) then
            if (irep == iref) then
               write(note,*) '   < reference'
            else
               write(note,*) ''
            end if
            call CALE(res_flat, reshape(X_RKlib(:,:,irep), shape(res_flat)), BBTW_flat, .false.)
            write(*,'(I7,F10.2,3(1X,F18.12),F10.2," s",A)') irep, irep*Tstep, norm2(X_RKlib(:,:,irep))/N, &
                        & norm2(X_RKlib(:,:,irep)-Xref_BS)/N, norm2(res_flat)/N, etime, trim(note) 
         end if
      enddo
      Xref_RK(:,:) = X_RKlib(:,:,iref)
      return
   end subroutine run_lyap_reference_RK

end module Ginzburg_Landau_RK