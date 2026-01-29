module Ginzburg_Landau_Tests_Riccati
   ! Standard Library.
   use stdlib_optval, only : optval
   use stdlib_linalg, only : diag, svd, svdvals
   use stdlib_io_npy, only : save_npy, load_npy
   ! LightKrylov for linear algebra.
   use LightKrylov
   use LightKrylov_AbstractVectors ! linear_combination
   use LightKrylov_Utils ! svd, sqrtm
   ! LightROM
   use LightROM_Utils ! Balancing_Transformation
   ! Riccati Solver
   use LightROM_RiccatiSolvers
   use LightROM_RiccatiUtils
   ! Ginzburg Landau
   use Ginzburg_Landau_Base
   use Ginzburg_Landau_Operators
   use Ginzburg_Landau_RKlib
   use Ginzburg_Landau_Utils
   implicit none

   integer,       parameter :: rkmax = 80
   integer,       parameter :: rk_X0_riccati = 10

   private :: this_module
   public  :: rkmax, rk_X0_riccati

   character(len=*), parameter :: this_module = 'Ginzburg_Landau_Tests_Riccati'

contains

   !-------------------------------------------------------------------------------------------
   !
   ! RICCATI EQUATION
   !
   !-------------------------------------------------------------------------------------------

   subroutine run_riccati_reference_RK(LTI, Xref, Xref_RK, U0, S0, Tend, nrep, iref)
      ! LTI system
      type(lti_system),              intent(inout) :: LTI
      ! Reference solution (SD)
      real(dp),                      intent(in)    :: Xref(N,N)
      ! Reference solution (RK)
      real(dp),                      intent(out)   :: Xref_RK(N,N)
      ! Initial condition
      type(state_vector),            intent(inout) :: U0(:)
      real(dp),                      intent(inout) :: S0(:,:)
      real(dp),                      intent(in)    :: Tend
      integer,                       intent(in)    :: nrep
      integer,                       intent(in)    :: iref

      ! Internals
      type(rk_riccati) ,             allocatable   :: RK_propagator
      type(state_matrix)                           :: X_mat(2)
      real(dp),                      allocatable   :: X_RK(:,:,:)
      real(dp)                                     :: U0_mat(N,N)
      integer                                      :: irep
      real(dp)                                     :: etime, Tstep
      ! OUTPUT
      real(dp)                                     :: X_out(N,N)
      character(len=128)                           :: note
      ! timer
      integer :: clock_rate, clock_start, clock_stop

      call system_clock(count_rate=clock_rate)

      Tstep = Tend/nrep
      allocate(X_RK(N, N, nrep))
      ! initialize exponential propagator
      RK_propagator = RK_riccati(Tstep)
      ! Set initial condition for RK
      call reconstruct_solution(X_out, U0, S0)
      call set_state(X_mat(1:1), X_out, 'Set initial condition')
      write(*,'(A7,A10,A19,A19,A19,A12)') ' RKlib:','Tend','| X_RK |/N', '| X_RK - X_SD |/N', '| res_RK |/N','etime'
      write(*,*) '-------------------------------------------------------------------------------------'
      write(*,'(I7,F10.4,3(1X,E18.6),F10.4," s",A)') 0, 0.0, norm2(X_out)/N, norm2(X_out - Xref)/N, &
                                                            & norm2(CARE(X_out, CTQcCW, BRinvBTW))/N, 0.0, ''
      do irep = 1, nrep
         call system_clock(count=clock_start)     ! Start Timer
         ! integrate
         call RK_propagator%matvec(X_mat(1), X_mat(2))
         call system_clock(count=clock_stop)      ! Stop Timer
         etime = real(clock_stop-clock_start)/real(clock_rate)
         ! recover output
         call get_state(X_RK(:,:,irep), X_mat(2:2), 'Extract RK solution')
         ! replace input
         call set_state(X_mat(1:1), X_RK(:,:,irep), 'Reset initial condition')
         ! compute residual
         if (irep == iref) then
            write(note,*) '   < reference'
         else
            write(note,*) ''
         end if
         write(*,'(I7,F10.4,3(1X,E18.6),F10.4," s",A)') irep, irep*Tstep, norm2(X_RK(:,:,irep))/N, &
                     & norm2(X_RK(:,:,irep)-Xref)/N, norm2(CARE(X_RK(:,:,irep), CTQcCW, BRinvBTW))/N, etime, trim(note) 
      enddo
      Xref_RK(:,:) = X_RK(:,:,iref)
      print *, ''
      print '(A,F16.12)', '  |  X_RK  |/N = ', norm2(Xref_RK)/N
      print '(A,F16.12)', '  | res_RK |/N = ', norm2(CARE(Xref_RK, CTQcCW, BRinvBTW))/N
   end subroutine run_riccati_reference_RK

   subroutine run_riccati_DLRA_test(LTI, Xref, Xref_RK, U0, S0, Tend, dtv, rkv, TOv, nprint, home, if_save_output)
      ! LTI system
      type(lti_system),              intent(inout) :: LTI
      ! Reference solution (BS)
      real(dp),                      intent(in)    :: Xref(N,N)
      ! Reference solution (RK)
      real(dp),                      intent(in)    :: Xref_RK(N,N)
      ! Initial condition
      type(state_vector),            intent(inout) :: U0(:)
      real(dp),                      intent(inout) :: S0(:,:)
      real(dp),                      intent(in)    :: Tend
      ! vector of dt values
      real(dp),                      intent(in)    :: dtv(:)
      ! vector of rank values
      integer,                       intent(in)    :: rkv(:)
      ! vector of torders
      integer,                       intent(in)    :: TOv(:)
      ! number of singular values to print
      integer,                       intent(in)    :: nprint
      character(len=128),            intent(in)    :: home
      logical,                       intent(in)    :: if_save_output

      ! Internals
      type(LR_state),                allocatable   :: X
      type(state_matrix)                           :: X_mat(2)
      integer                                      :: info, i, j, k, rk, torder, irep, nrep, nsteps
      real(dp)                                     :: etime, tau, Tstep, Ttot
      ! OUTPUT
      real(dp)                                     :: X_out(N,N)

      ! SVD                         
      integer                                      :: is, ie
      integer,                         parameter   :: irow = 8
      real(dp),                        allocatable :: svals(:)
      character(len=128)                           :: note
      ! DLRA options
      type(dlra_opts)                              :: opts
      ! timer
      integer                                      :: clock_rate, clock_start, clock_stop
      
      Tstep = one_rdp

      note = 'P'

      ! basic opts
      opts = dlra_opts(chktime=one_rdp, inc_tol=atol_dp, if_rank_adaptive=.false.)

      call system_clock(count_rate=clock_rate)

      print '(A16,A8,A4,A10,A8,A10,4(A19),A12)', 'DLRA:','  rk',' TO','dt','steps','Tend', &
                     & '| X_D |/N', '| X_D - X_RK |/N', '| X_D - X_SD |/N', '| res_D |/N', 'etime'
      X = LR_state()
      do i = 1, size(rkv)
         rk = rkv(i)
         do j = 1, size(TOv)
            torder = TOv(j)
            do k = 1, size(dtv)
               tau = dtv(k)
               nsteps = nint(Tend/tau)
               
               ! set solver options
               opts%mode = torder

               ! Initialize low-rank representation with rank rk
               call X%initialize_LR_state(U0, S0, rk, rkmax, .false.)

               ! run integrator
               call system_clock(count=clock_start)     ! Start Timer
               call projector_splitting_DLRA_riccati_integrator(X, LTI%prop, LTI%B, LTI%CT, Qc, Rinv, Tend, tau, info, &
                                                            & exptA=exptA, iftrans=.true., options=opts)
               call system_clock(count=clock_stop)      ! Stop Timer
               etime = real(clock_stop-clock_start)/real(clock_rate)
               ! Reconstruct solution
               call reconstruct_solution(X_out, X)
               write(*,'(I4," ",A4,1X,A6,I8," TO",I1,F10.6,I8,F10.4,4(E19.8),F10.2," s")') 1, note, 'OUTPUT', &
               & rk, torder, tau, nsteps, Tend, &
               & norm2(X_out)/N, norm2(X_out - Xref_RK)/N, norm2(X_out - Xref)/N, &
               & norm2(CARE(X_out, CTQcCW, BRinvBTW))/N, etime
               deallocate(X%U); deallocate(X%S)
            end do
            print *, ''
         end do
      end do
      if (nprint > 0) then
         svals = svdvals(Xref)
         do i = 1, ceiling(nprint*one_rdp/irow)
            is = (i-1)*irow+1; ie = i*irow
            print '(1X,A,I2,A,I2,*(1X,F16.12))', 'SVD(X_SD) ', is, '-', ie, ( svals(j), j = is, ie )
         end do
         print *, ''
         svals = svdvals(Xref_RK)
         do i = 1, ceiling(nprint*one_rdp/irow)
            is = (i-1)*irow+1; ie = i*irow
            print '(1X,A,I2,A,I2,*(1X,F16.12))', 'SVD(X_RK) ', is, '-', ie, ( svals(j), j = is, ie )
         end do
         print *, ''
         svals = svdvals(X_out)
         do i = 1, ceiling(nprint*one_rdp/irow)
            is = (i-1)*irow+1; ie = i*irow
            print '(1X,A,I2,A,I2,*(1X,F16.12))', 'SVD(X_D ) ', is, '-', ie, ( svals(j), j = is, ie )
         end do
      end if
   end subroutine run_riccati_DLRA_test

   subroutine run_riccati_DLRArk_test(LTI, Xref, Xref_RK, U0, S0, Tend, dtv, TOv, tolv, nprint, home, if_save_output)
      ! LTI system
      type(lti_system),              intent(inout) :: LTI
      ! Reference solution (BS)
      real(dp),                      intent(in)    :: Xref(N,N)
      ! Reference solution (RK)
      real(dp),                      intent(in)    :: Xref_RK(N,N)
      ! Initial condition
      type(state_vector),            intent(inout) :: U0(:)
      real(dp),                      intent(inout) :: S0(:,:)
      real(dp),                      intent(in)    :: Tend
      ! vector of dt values
      real(dp),                      intent(in)    :: dtv(:)
      ! vector of torders
      integer,                       intent(in)    :: TOv(:)
      ! vector of adaptation tolerance
      real(dp),                      intent(in)    :: tolv(:)
      ! number of singular values to print
      integer,                       intent(in)    :: nprint
      character(len=128),            intent(in)    :: home
      logical,                       intent(in)    :: if_save_output

      ! Internals
      type(LR_state),                allocatable   :: X
      type(state_matrix)                           :: X_mat(2)
      integer                                      :: info, i, j, k, rk, torder, irep, nsteps
      real(dp)                                     :: etime, tau
      ! OUTPUT
      real(dp)                                     :: X_out(N,N)
      ! timer
      integer                                      :: is, ie
      integer,                         parameter   :: irow = 8
      real(dp),                        allocatable :: svals(:)
      character(len=128)                           :: note
      ! DLRA options
      type(dlra_opts)                              :: opts
      ! timer
      integer                                      :: clock_rate, clock_start, clock_stop

      ! basic opts
      opts = dlra_opts(chktime=one_rdp, inc_tol=atol_dp, if_rank_adaptive=.true.)

      call system_clock(count_rate=clock_rate)

      note = 'P'
    
      print '(A16,A8,A4,A10,A8,A10,4(A19),A12)', 'DLRA:','rk_end',' TO','dt','steps','Tend', &
               & '| X_D |/N', '| X_D - X_RK |/N', '| X_D - X_SD |/N', '| res_D |/N', 'etime'
      rk = rk_X0_riccati
      X = LR_state()
      do i = 1, size(tolv)
         opts%tol = tolv(i)
         print '(A,E9.2)', ' SVD tol = ', opts%tol
         print *, ''
         do j = 1, size(TOv)
            torder = TOv(j)
            do k = 1, size(dtv)
               tau = dtv(k)
               nsteps = nint(Tend/tau)
               
               ! set solver options
               opts%mode = torder

               ! Initialize low-rank representation with rank rk
               call X%initialize_LR_state(U0, S0, rk, rkmax, opts%if_rank_adaptive)
               X%tot_time = 0.0_dp
               X%time     = 0.0_dp
               X%step     = 0

               ! run integrator
               call system_clock(count=clock_start)     ! Start Timer
               call projector_splitting_DLRA_riccati_integrator(X, LTI%prop, LTI%B, LTI%CT, Qc, Rinv, Tend, tau, info, &
                                                               & exptA=exptA, iftrans=.true., options=opts)
               call system_clock(count=clock_stop)      ! Stop Timer
               etime = real(clock_stop-clock_start)/real(clock_rate)
               ! Reconstruct solution
               call reconstruct_solution(X_out, X)
               write(*,'(I4," ",A4,1X,A6,I8," TO",I1,F10.6,I8,F10.4,4(E19.8),F10.2," s")') 1, note, 'OUTPUT', &
                                 & X%rk, torder, tau, nsteps, Tend, &
                                 & norm2(X_out)/N, norm2(X_out - Xref_RK)/N, norm2(X_out - Xref)/N, &
                                 & norm2(CARE(X_out, CTQcCW, BRinvBTW))/N, etime
               deallocate(X%U); deallocate(X%S)
            end do
            print *, ''
         end do
         svals = svdvals(Xref)
         do k = 1, ceiling(nprint*one_rdp/irow)
            is = (k-1)*irow+1; ie = k*irow
            print '(1X,A,I2,A,I2,*(1X,F16.12))', 'SVD(X_BS) ', is, '-', ie, ( svals(j), j = is, ie )
         end do
         print *, ''
         svals = svdvals(Xref_RK)
         do k = 1, ceiling(nprint*one_rdp/irow)
            is = (k-1)*irow+1; ie = k*irow
            print '(1X,A,I2,A,I2,*(1X,F16.12))', 'SVD(X_RK) ', is, '-', ie, ( svals(j), j = is, ie )
         end do
         print *, ''
         svals = svdvals(X_out)
         do k = 1, ceiling(nprint*one_rdp/irow)
            is = (k-1)*irow+1; ie = k*irow
            print '(1X,A,I2,A,I2,*(1X,F16.12))', 'SVD(X_D ) ', is, '-', ie, ( svals(j), j = is, ie )
         end do
      end do

   end subroutine run_riccati_DLRArk_test

end module Ginzburg_Landau_Tests_Riccati