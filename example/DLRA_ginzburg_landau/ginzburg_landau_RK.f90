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
   public  :: run_lyap_RK

   character*128, parameter :: this_module = 'Ginzburg_Landau_RK'

contains

   subroutine run_lyap_RK(LTI, U0, S0, Tend, ifsave, ifverb)
      ! LTI system
      type(lti_system),              intent(inout) :: LTI
      ! Initial condition
      type(state_vector),            intent(inout) :: U0(:)
      real(wp),                      intent(inout) :: S0(:,:)
      real(wp),                      intent(in)    :: Tend
      ! Optional
      logical, optional,             intent(in)    :: ifsave
      logical                                      :: if_save_npy
      logical, optional,             intent(in)    :: ifverb
      logical                                      :: verb

      ! Internals
      type(LR_state),                allocatable   :: X_state
      type(rk_lyapunov),             allocatable   :: RK_propagator
      type(state_matrix)                           :: X_mat(2)
      real(wp),                      allocatable   :: X_RKlib(:,:,:)
      real(wp)                                     :: X_mat_ref(N,N)
      real(wp)                                     :: U0_mat(N,N)
      integer                                      :: info
      integer                                      :: i, j, ito, rk, nsteps, torder
      integer                                      :: irep, nrep
      real(wp)                                     :: etime, tau
      ! OUTPUT
      real(wp)                                     :: U_out(N,N)
      real(wp)                                     :: X_out(N,N)
      real(wp)                                     :: Bmat(N,2)
      real(wp)                                     :: tmp(N,2)
      real(wp),           allocatable              :: U_load(:,:)
      real(wp)                                     :: X_mat_flat(N**2)
      real(wp)                                     :: res_flat(N**2)
      character*128      :: oname
      character*128      :: onameU
      character*128      :: onameS
      integer            :: iostatus
      ! timer
      integer            :: clock_rate, clock_start, clock_stop
      logical            :: existfile, if_load_data
      ! DLRA opts
      type(dlra_opts)                              :: opts

      if_save_npy  = optval(ifsave, .false.)
      verb         = optval(ifverb, .false.)

      if_load_data = .true.

      call system_clock(count_rate=clock_rate)

      nrep = 80
      allocate(X_RKlib(N, N, nrep))
      if (if_load_data) then
         print *, ''
         print *, 'Load the solution to the differential Lyapunov equation computed with RKlib:'
         print *, ''
      else
         print *, ''
         print *, 'Compute the solution to the differential Lyapunov equation with RKlib:'
         print *, ''
         ! initialize exponential propagator
         RK_propagator = RK_lyapunov(Tend)
         ! Set random initial condition
         call get_state(U0_mat(:,1:rk_X0), U0)
         if (if_save_npy) then
            X_out = matmul( U0_mat(:,1:rk_X0), matmul( S0, transpose(U0_mat(:,1:rk_X0)) ) )
            ! Save forcing RK
            write(oname,'("local/data_GL_lyapconv_BBTW_RK_n",I4.4,".npy")') nx
            call save_data(oname, reshape(BBTW_flat, (/ N,N /)))
            ! Save initial condition
            write(oname,'("local/data_GL_lyapconv_X0_RK_n",I4.4,".npy")') nx
            call save_data(oname, X_out)
            ! Save forcing DLRA
            Bmat = 0.0_wp
            call get_state(Bmat, LTI%B(1:rk_b))
            X_out = matmul(Bmat, transpose(Bmat))
            write(oname,'("local/data_GL_lyapconv_BBTW_DLRA_n",I4.4,".npy")') nx
            call save_data(oname, X_out)
         end if
         X_out = matmul( U0_mat(:,1:rk_X0), matmul( S0, transpose(U0_mat(:,1:rk_X0)) ) )
         ! Set initial condition for RK
         call set_state(X_mat(1:1), X_out)
      end if 
      write(*,'(A8,A18,A18,A18,A18)') ' RKlib:','Tend','|| X_RK ||_2/N', '|| res ||_2/N','Elapsed time'
      write(*,*) ' ------------------------------------------------------------------------------'
      do irep = 1, nrep
         write(oname,'("local/data_GL_lyapconv_X_RK_n",I4.4,"_r",I3.3,".npy")') nx, irep
         if (if_load_data) then
            call load_data(oname, U_load)
            X_RKlib(:,:,irep) = U_load
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
         call CALE(res_flat, reshape(X_RKlib(:,:,irep), shape(res_flat)), BBTW_flat, .false.)
         write(*,'(I8,F18.4,E18.8,E18.8,F16.4," s")') irep, irep*Tend, norm2(X_RKlib(:,:,irep))/N, &
                        & norm2(res_flat)/N, etime 
      enddo
      return
   end subroutine run_lyap_RK

end module Ginzburg_Landau_RK