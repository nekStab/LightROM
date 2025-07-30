module Ginzburg_Landau_Tests
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
   use Ginzburg_Landau_RKlib
   use Ginzburg_Landau_Utils
   implicit none

   integer,       parameter :: rkmax = 80
   integer,       parameter :: rk_X0 = 10

   private :: this_module
   public  :: rkmax, rk_X0

   character(len=*), parameter :: this_module = 'Ginzburg_Landau_Tests'

contains

   !
   ! LYAPUNOV EQUATION
   !

   subroutine run_lyap_reference_RK(LTI, Xref, Xref_RK, U0, S0, Tend, nrep, iref, adjoint)
      ! LTI system
      type(lti_system),              intent(inout) :: LTI
      ! Reference solution (BS)
      real(wp),                      intent(in)    :: Xref(N,N)
      ! Reference solution (RK)
      real(wp),                      intent(out)   :: Xref_RK(N,N)
      ! Initial condition
      type(state_vector),            intent(inout) :: U0(:)
      real(wp),                      intent(inout) :: S0(:,:)
      real(wp),                      intent(in)    :: Tend
      integer,                       intent(in)    :: nrep
      integer,                       intent(in)    :: iref
      logical,                       intent(in)    :: adjoint

      ! Internals
      type(rk_lyapunov),             allocatable   :: RK_propagator
      type(state_matrix)                           :: X_mat(2)
      real(wp),                      allocatable   :: X_RK(:,:,:)
      real(wp)                                     :: U0_mat(N,N)
      integer                                      :: irep
      real(wp)                                     :: etime, Tstep
      ! OUTPUT
      real(wp)                                     :: X_out(N,N)
      character(len=128)                           :: note
      ! timer
      integer :: clock_rate, clock_start, clock_stop

      call system_clock(count_rate=clock_rate)

      Tstep = Tend/nrep
      allocate(X_RK(N, N, nrep))
      ! initialize exponential propagator
      RK_propagator = RK_lyapunov(Tstep)
      ! Set initial condition for RK
      call reconstruct_solution(X_out, U0, S0)
      call set_state(X_mat(1:1), X_out, 'Set initial condition')
      write(*,'(A7,A10,A19,A19,A19,A12)') ' RKlib:','Tend','| X_RK |/N', '| X_RK - X_BS |/N', '| res_RK |/N','etime'
      write(*,*) '-------------------------------------------------------------------------------------'
      write(*,'(I7,F10.4,3(1X,E18.6),F10.4," s",A)') 0, 0.0, norm2(X_out)/N, norm2(X_out - Xref)/N, &
                                                            & norm2(CALE(X_out, adjoint))/N, 0.0, ''
      do irep = 1, nrep
         call system_clock(count=clock_start)     ! Start Timer
         ! integrate
         if (adjoint) then
            call RK_propagator%rmatvec(X_mat(1), X_mat(2))
         else
            call RK_propagator%matvec(X_mat(1), X_mat(2))
         end if
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
                     & norm2(X_RK(:,:,irep)-Xref)/N, norm2(CALE(X_RK(:,:,irep), adjoint))/N, etime, trim(note) 
      enddo
      Xref_RK(:,:) = X_RK(:,:,iref)
      print *, ''
      print '(A,F16.12)', '  |  X_RK  |/N = ', norm2(Xref_RK)/N
      print '(A,F16.12)', '  | res_RK |/N = ', norm2(CALE(Xref_RK, adjoint))/N
      return
   end subroutine run_lyap_reference_RK

   subroutine run_lyap_DLRA_test(LTI, Xref, Xref_RK, U0, S0, Tend, dtv, rkv, TOv, nprint, adjoint, home, if_save_output)
      ! LTI system
      type(lti_system),              intent(inout) :: LTI
      ! Reference solution (BS)
      real(wp),                      intent(in)    :: Xref(N,N)
      ! Reference solution (RK)
      real(wp),                      intent(in)    :: Xref_RK(N,N)
      ! Initial condition
      type(state_vector),            intent(inout) :: U0(:)
      real(wp),                      intent(inout) :: S0(:,:)
      real(wp),                      intent(in)    :: Tend
      ! vector of dt values
      real(wp),                      intent(in)    :: dtv(:)
      ! vector of rank values
      integer,                       intent(in)    :: rkv(:)
      ! vector of torders
      integer,                       intent(in)    :: TOv(:)
      ! number of singular values to print
      integer,                       intent(in)    :: nprint
      logical,                       intent(in)    :: adjoint
      character(len=128),            intent(in)    :: home
      logical,                       intent(in)    :: if_save_output

      ! Internals
      type(LR_state),                allocatable   :: X
      type(state_matrix)                           :: X_mat(2)
      integer                                      :: info, i, j, k, rk, torder, irep, nrep, nsteps
      real(wp)                                     :: etime, tau, Tstep, Ttot
      ! OUTPUT
      real(wp)                                     :: X_out(N,N)

      ! SVD                         
      integer                                      :: is, ie
      integer,                         parameter   :: irow = 8
      real(wp),                        allocatable :: svals(:)
      character(len=128)                           :: note
      ! DLRA options
      type(dlra_opts)                              :: opts
      ! timer
      integer                                      :: clock_rate, clock_start, clock_stop
      
      Tstep = 1.0_wp

      if (adjoint) then
         note = 'Yobs'
      else
         note = 'Xctl'
      end if

      ! basic opts
      opts = dlra_opts(chktime=1.0_wp, inc_tol=atol_dp, if_rank_adaptive=.false.)

      call system_clock(count_rate=clock_rate)

      print '(A16,A8,A4,A10,A8,A10,4(A19),A12)', 'DLRA:','  rk',' TO','dt','steps','Tend', &
                     & '| X_D |/N', '| X_D - X_RK |/N', '| X_D - X_BS |/N', '| res_D |/N', 'etime'
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
               if (adjoint) then
                  call projector_splitting_DLRA_lyapunov_integrator(X, LTI%prop, LTI%CT, Tend, tau, info, &
                                                            & exptA=exptA, iftrans=.true., options=opts)
               else
                  call projector_splitting_DLRA_lyapunov_integrator(X, LTI%prop, LTI%B, Tend, tau, info, &
                                                            & exptA=exptA, iftrans=.false., options=opts)
               end if
               call system_clock(count=clock_stop)      ! Stop Timer
               etime = real(clock_stop-clock_start)/real(clock_rate)
               ! Reconstruct solution
               call reconstruct_solution(X_out, X)
               write(*,'(I4," ",A4,1X,A6,I8," TO",I1,F10.6,I8,F10.4,4(E19.8),F10.2," s")') 1, note, 'OUTPUT', &
                                 & rk, torder, tau, nsteps, Tend, &
                                 & norm2(X_out)/N, norm2(X_out - Xref_RK)/N, norm2(X_out - Xref)/N, &
                                 & norm2(CALE(X_out, adjoint))/N, etime
               deallocate(X%U); deallocate(X%S)
            end do
            print *, ''
         end do
      end do
      if (nprint > 0) then
         svals = svdvals(Xref)
         do i = 1, ceiling(nprint*1.0_wp/irow)
            is = (i-1)*irow+1; ie = i*irow
            print '(1X,A,I2,A,I2,*(1X,F16.12))', 'SVD(X_BS) ', is, '-', ie, ( svals(j), j = is, ie )
         end do
         print *, ''
         svals = svdvals(Xref_RK)
         do i = 1, ceiling(nprint*1.0_wp/irow)
            is = (i-1)*irow+1; ie = i*irow
            print '(1X,A,I2,A,I2,*(1X,F16.12))', 'SVD(X_RK) ', is, '-', ie, ( svals(j), j = is, ie )
         end do
         print *, ''
         svals = svdvals(X_out)
         do i = 1, ceiling(nprint*1.0_wp/irow)
            is = (i-1)*irow+1; ie = i*irow
            print '(1X,A,I2,A,I2,*(1X,F16.12))', 'SVD(X_D ) ', is, '-', ie, ( svals(j), j = is, ie )
         end do
      end if

   end subroutine run_lyap_DLRA_test

   subroutine run_lyap_DLRArk_test(LTI, Xref, Xref_RK, U0, S0, Tend, dtv, TOv, tolv, nprint, adjoint, home, if_save_output)
      ! LTI system
      type(lti_system),              intent(inout) :: LTI
      ! Reference solution (BS)
      real(wp),                      intent(in)    :: Xref(N,N)
      ! Reference solution (RK)
      real(wp),                      intent(in)    :: Xref_RK(N,N)
      ! Initial condition
      type(state_vector),            intent(inout) :: U0(:)
      real(wp),                      intent(inout) :: S0(:,:)
      real(wp),                      intent(in)    :: Tend
      ! vector of dt values
      real(wp),                      intent(in)    :: dtv(:)
      ! vector of torders
      integer,                       intent(in)    :: TOv(:)
      ! vector of adaptation tolerance
      real(wp),                      intent(in)    :: tolv(:)
      ! number of singular values to print
      integer,                       intent(in)    :: nprint
      logical,                       intent(in)    :: adjoint
      character(len=128),            intent(in)    :: home
      logical,                       intent(in)    :: if_save_output

      ! Internals
      type(LR_state),                allocatable   :: X
      type(state_matrix)                           :: X_mat(2)
      integer                                      :: info, i, j, k, rk, torder, irep, nsteps
      real(wp)                                     :: etime, tau
      ! OUTPUT
      real(wp)                                     :: X_out(N,N)
      ! timer
      integer                                      :: is, ie
      integer,                         parameter   :: irow = 8
      real(wp),                        allocatable :: svals(:)
      character(len=128)                           :: note
      ! DLRA options
      type(dlra_opts)                              :: opts
      ! timer
      integer                                      :: clock_rate, clock_start, clock_stop

      ! basic opts
      opts = dlra_opts(chktime=1.0_wp, inc_tol=atol_dp, if_rank_adaptive=.true.)

      call system_clock(count_rate=clock_rate)

      if (adjoint) then
         note = 'Yobs'
      else
         note = 'Xctl'
      end if
 
      print '(A16,A8,A4,A10,A8,A10,4(A19),A12)', 'DLRA:','rk_end',' TO','dt','steps','Tend', &
               & '| X_D |/N', '| X_D - X_RK |/N', '| X_D - X_BS |/N', '| res_D |/N', 'etime'
      rk = rk_X0
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
               if (adjoint) then
                  call projector_splitting_DLRA_lyapunov_integrator(X, LTI%prop, LTI%CT, Tend, tau, info, &
                                                               & exptA=exptA, iftrans=.true., options=opts)
               else
                  call projector_splitting_DLRA_lyapunov_integrator(X, LTI%prop, LTI%B, Tend, tau, info, &
                                                               & exptA=exptA, iftrans=.false., options=opts)
               end if
               call system_clock(count=clock_stop)      ! Stop Timer
               etime = real(clock_stop-clock_start)/real(clock_rate)
               ! Reconstruct solution
               call reconstruct_solution(X_out, X)
               write(*,'(I4," ",A4,1X,A6,I8," TO",I1,F10.6,I8,F10.4,4(E19.8),F10.2," s")') 1, note, 'OUTPUT', &
                                 & X%rk, torder, tau, nsteps, Tend, &
                                 & norm2(X_out)/N, norm2(X_out - Xref_RK)/N, norm2(X_out - Xref)/N, &
                                 & norm2(CALE(X_out, adjoint))/N, etime
               deallocate(X%U); deallocate(X%S)
            end do
            print *, ''
         end do
         svals = svdvals(Xref)
         do k = 1, ceiling(nprint*1.0_wp/irow)
            is = (k-1)*irow+1; ie = k*irow
            print '(1X,A,I2,A,I2,*(1X,F16.12))', 'SVD(X_BS) ', is, '-', ie, ( svals(j), j = is, ie )
         end do
         print *, ''
         svals = svdvals(Xref_RK)
         do k = 1, ceiling(nprint*1.0_wp/irow)
            is = (k-1)*irow+1; ie = k*irow
            print '(1X,A,I2,A,I2,*(1X,F16.12))', 'SVD(X_RK) ', is, '-', ie, ( svals(j), j = is, ie )
         end do
         print *, ''
         svals = svdvals(X_out)
         do k = 1, ceiling(nprint*1.0_wp/irow)
            is = (k-1)*irow+1; ie = k*irow
            print '(1X,A,I2,A,I2,*(1X,F16.12))', 'SVD(X_D ) ', is, '-', ie, ( svals(j), j = is, ie )
         end do
      end do

   end subroutine run_lyap_DLRArk_test

   !
   ! RICCATI EQUATION
   !

   subroutine run_ricc_reference_RK(LTI, Xref, Xref_RK, U0, S0, Tend, nrep, iref, adjoint)
      ! LTI system
      type(lti_system),              intent(inout) :: LTI
      ! Reference solution (SD)
      real(wp),                      intent(in)    :: Xref(N,N)
      ! Reference solution (RK)
      real(wp),                      intent(out)   :: Xref_RK(N,N)
      ! Initial condition
      type(state_vector),            intent(inout) :: U0(:)
      real(wp),                      intent(inout) :: S0(:,:)
      real(wp),                      intent(in)    :: Tend
      integer,                       intent(in)    :: nrep
      integer,                       intent(in)    :: iref
      logical,                       intent(in)    :: adjoint

      ! Internals
      type(rk_riccati),              allocatable   :: RK_propagator
      type(state_matrix)                           :: X_mat(2)
      real(wp),                      allocatable   :: X_RK(:,:,:)
      real(wp)                                     :: U0_mat(N,N)
      integer                                      :: irep
      real(wp)                                     :: etime, Tstep
      ! OUTPUT
      real(wp)                                     :: X_out(N,N)
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
      write(*,'(A7,A10,A19,A19,A19,A12)') ' RKlib:','Tend','| X_RK |/N', '| X_RK - X_BS |/N', '| res_RK |/N','etime'
      write(*,*) '-------------------------------------------------------------------------------------'
      write(*,'(I7,F10.4,3(1X,E18.6),F10.4," s",A)') 0, 0.0, norm2(X_out)/N, norm2(X_out - Xref)/N, &
                                                            & norm2(CARE(X_out, CTQcCW, BRinvBTW))/N, 0.0, ''
      do irep = 1, nrep
         call system_clock(count=clock_start)     ! Start Timer
         ! integrate
         if (adjoint) then
            call RK_propagator%rmatvec(X_mat(1), X_mat(2))
         else
            call RK_propagator%matvec(X_mat(1), X_mat(2))
         end if
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
   end subroutine run_ricc_reference_RK

   subroutine run_ricc_DLRA_test(LTI, Xref, Xref_RK, U0, S0, Tend, dtv, rkv, TOv, nprint, adjoint, home, if_save_output)
      ! LTI system
      type(lti_system),              intent(inout) :: LTI
      ! Reference solution (BS)
      real(wp),                      intent(in)    :: Xref(N,N)
      ! Reference solution (RK)
      real(wp),                      intent(in)    :: Xref_RK(N,N)
      ! Initial condition
      type(state_vector),            intent(inout) :: U0(:)
      real(wp),                      intent(inout) :: S0(:,:)
      real(wp),                      intent(in)    :: Tend
      ! vector of dt values
      real(wp),                      intent(in)    :: dtv(:)
      ! vector of rank values
      integer,                       intent(in)    :: rkv(:)
      ! vector of torders
      integer,                       intent(in)    :: TOv(:)
      ! number of singular values to print
      integer,                       intent(in)    :: nprint
      logical,                       intent(in)    :: adjoint
      character(len=128),            intent(in)    :: home
      logical,                       intent(in)    :: if_save_output

      ! Internals
      type(LR_state),                allocatable   :: X
      type(state_matrix)                           :: X_mat(2)
      integer                                      :: info, i, j, k, rk, torder, irep, nrep, nsteps
      real(wp)                                     :: etime, tau, Tstep, Ttot
      ! OUTPUT
      real(wp)                                     :: X_out(N,N)

      ! SVD                         
      integer                                      :: is, ie
      integer,                         parameter   :: irow = 8
      real(wp),                        allocatable :: svals(:)
      character(len=128)                           :: note
      ! DLRA options
      type(dlra_opts)                              :: opts
      ! timer
      integer                                      :: clock_rate, clock_start, clock_stop
      
      Tstep = 1.0_wp

      note = 'P'

      ! basic opts
      opts = dlra_opts(chktime=1.0_wp, inc_tol=atol_dp, if_rank_adaptive=.false.)

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
                                                            & exptA=exptA, iftrans=adjoint, options=opts)
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
         do i = 1, ceiling(nprint*1.0_wp/irow)
            is = (i-1)*irow+1; ie = i*irow
            print '(1X,A,I2,A,I2,*(1X,F16.12))', 'SVD(X_BS) ', is, '-', ie, ( svals(j), j = is, ie )
         end do
         print *, ''
         svals = svdvals(Xref_RK)
         do i = 1, ceiling(nprint*1.0_wp/irow)
            is = (i-1)*irow+1; ie = i*irow
            print '(1X,A,I2,A,I2,*(1X,F16.12))', 'SVD(X_RK) ', is, '-', ie, ( svals(j), j = is, ie )
         end do
         print *, ''
         svals = svdvals(X_out)
         do i = 1, ceiling(nprint*1.0_wp/irow)
            is = (i-1)*irow+1; ie = i*irow
            print '(1X,A,I2,A,I2,*(1X,F16.12))', 'SVD(X_D ) ', is, '-', ie, ( svals(j), j = is, ie )
         end do
      end if
   end subroutine run_ricc_DLRA_test

end module Ginzburg_Landau_Tests