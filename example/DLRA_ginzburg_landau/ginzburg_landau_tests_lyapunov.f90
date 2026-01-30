module Ginzburg_Landau_Tests_Lyapunov
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
   ! Lyapunov Solver
   use LightROM_LyapunovSolvers
   use LightROM_LyapunovUtils
   ! Ginzburg Landau
   use Ginzburg_Landau_Base
   use Ginzburg_Landau_Operators
   use Ginzburg_Landau_RKlib
   use Ginzburg_Landau_Utils
   implicit none

   integer,       parameter :: rkmax = 80
   integer,       parameter :: rk_X0_lyapunov = 10

   private :: this_module
   public  :: rkmax, rk_X0_lyapunov

   character(len=*), parameter :: this_module = 'Ginzburg_Landau_Tests_Lyapunov'

contains

   !-------------------------------------------------------------------------------------------
   !
   ! LYAPUNOV EQUATION
   !
   !-------------------------------------------------------------------------------------------

   subroutine run_lyapunov_reference_RK(LTI, Xref, Xref_RK, U0, S0, Tend, nrep, iref, adjoint)
      ! LTI system
      type(lti_system),              intent(inout) :: LTI
      ! Reference solution (BS)
      real(dp),                      intent(in)    :: Xref(N,N)
      ! Reference solution (RK)
      real(dp),                      intent(out)   :: Xref_RK(N,N)
      ! Initial condition
      type(state_vector),            intent(inout) :: U0(:)
      real(dp),                      intent(inout) :: S0(:,:)
      real(dp),                      intent(in)    :: Tend
      integer,                       intent(in)    :: nrep
      integer,                       intent(in)    :: iref
      logical,                       intent(in)    :: adjoint

      ! Internals
      type(rk_lyapunov),             allocatable   :: RK_propagator
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

      character(len=*), parameter :: case = 'RKLIB'
      character(len=*), parameter :: eq   = 'lyap'

      call system_clock(count_rate=clock_rate)

      Tstep = Tend/nrep
      allocate(X_RK(N, N, nrep))
      ! initialize exponential propagator
      RK_propagator = RK_lyapunov(Tstep)
      ! Set initial condition for RK
      call reconstruct_solution(X_out, U0, S0)
      call set_state(X_mat(1:1), X_out, 'Set initial condition')
      call print_header(case, eq)
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
         ! print information
         write(note,*) merge('   < reference', '              ', irep == iref)
         call print_rklib_output(eq, irep, Tstep, X_RK, Xref, etime, note, adjoint)
      enddo
      Xref_RK(:,:) = X_RK(:,:,iref)
      print *, ''
      print '(A,F16.12)', '  |  X_RK  |/N = ', norm2(Xref_RK)/N
      print '(A,F16.12)', '  | res_RK |/N = ', norm2(CALE(Xref_RK, adjoint))/N
      return
   end subroutine run_lyapunov_reference_RK

   subroutine run_lyapunov_DLRA_test(LTI, Xref, Xref_RK, U0, S0, Tend, dtv, rkv, TOv, nprint, adjoint, home, if_save_output)
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
      logical,                       intent(in)    :: adjoint
      character(len=128),            intent(in)    :: home
      logical,                       intent(in)    :: if_save_output
      ! IO
      real(wp),                      allocatable   :: U_load(:,:)
      logical :: exist_file
      
      ! Internals
      character(len=256)                           :: fname
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

      character(len=*), parameter :: case = 'DLRA_FIXED'
      character(len=*), parameter :: eq   = 'lyap'
      
      Tstep = one_rdp

      write(note,*) merge('Yobs', 'Xctl', adjoint)

      ! basic opts
      opts = dlra_opts(chktime=one_rdp, inc_tol=atol_dp, if_rank_adaptive=.false.)

      call system_clock(count_rate=clock_rate)

      call print_header(case, eq)
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

               fbase = make_filename(home, case, eq, note, rk, torder, tau, Tend)
               fname = fbase(:index(fbase,'.npy')-1)//'_U.npy'
               inquire(file=trim(fname), exist=exist_file)
               if (exist_file) then
                  ! U
                  call load_npy(trim(fname), U_load)
                  X%rk = size(U_load, 2)
                  allocate(X%U(X%rk), source=U0(1))
                  call set_state(X%U, U_load, 'load_from_file')
                  ! S
                  fname = fbase(:index(fbase,'.npy')-1)//'_S.npy'
                  call load_npy(trim(fname), U_load)
                  allocate(X%S(X%rk,X%rk), source=U_load)
               else
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
                  ! save output
                  if (if_save_output) then
                     fname = make_filename(home, case, eq, note, rk, torder, tau, Tend)
                     call save_LR_state_npy(fname, X, weight_mat)
                  end if
               end if
               ! Reconstruct solution
               call reconstruct_solution(X_out, X)
               ! print information
               call print_dlra_output(eq, rk, torder, tau, nsteps, Tend, X_out, Xref_RK, Xref, etime, note, adjoint)
               deallocate(X%U); deallocate(X%S)
            end do
            print *, ''
         end do
      end do
      if (nprint > 0) then
         call print_svdvals(Xref,     'X_BS', nprint, irow)
         call print_svdvals(Xref_RK,  'X_RK', nprint, irow)
         call print_svdvals(X_out,    'X_D ', nprint, irow)
      end if

   end subroutine run_lyapunov_DLRA_test

   subroutine run_lyapunov_DLRArk_test(LTI, Xref, Xref_RK, U0, S0, Tend, dtv, TOv, tolv, nprint, adjoint, home, if_save_output)
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
      logical,                       intent(in)    :: adjoint
      character(len=128),            intent(in)    :: home
      logical,                       intent(in)    :: if_save_output
      ! IO
      real(wp),                      allocatable   :: U_load(:,:)
      logical :: exist_file

      ! Internals
      character(len=256)                           :: fname
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

      character(len=*), parameter :: case = 'DLRA_ADAPT'
      character(len=*), parameter :: eq   = 'lyap'

      ! basic opts
      opts = dlra_opts(chktime=one_rdp, inc_tol=atol_dp, if_rank_adaptive=.true.)

      call system_clock(count_rate=clock_rate)

      write(note,*) merge('Yobs', 'Xctl', adjoint)
 
      call print_header(case, eq)
      rk = rk_X0_lyapunov
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

               fbase = make_filename(home, case, eq, note, rk, torder, tau, Tend)
               fname = fbase(:index(fbase,'.npy')-1)//'_U.npy'
               inquire(file=trim(fname), exist=exist_file)
               if (exist_file) then
                  ! U
                  call load_npy(trim(fname), U_load)
                  X%rk = size(U_load, 2)
                  allocate(X%U(X%rk), source=U0(1))
                  call set_state(X%U, U_load, 'load_from_file')
                  ! S
                  fname = fbase(:index(fbase,'.npy')-1)//'_S.npy'
                  call load_npy(trim(fname), U_load)
                  allocate(X%S(X%rk,X%rk), source=U_load)
               else
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
                  ! save output
                  if (if_save_output) then
                     fname = make_filename(home, case, eq, note, rk, torder, tau, Tend, opts%tol)
                     call save_LR_state_npy(fname, X, weight_mat)
                  end if
               end if
               ! Reconstruct solution
               call reconstruct_solution(X_out, X)
               ! print information
               call print_dlra_output(eq, X%rk, torder, tau, nsteps, Tend, X_out, Xref_RK, Xref, etime, note, adjoint)
               deallocate(X%U); deallocate(X%S)
            end do
            print *, ''
         end do
         call print_svdvals(Xref,     'X_BS', nprint, irow)
         call print_svdvals(Xref_RK,  'X_RK', nprint, irow)
         call print_svdvals(X_out,    'X_D ', nprint, irow)
      end do

   end subroutine run_lyapunov_DLRArk_test

end module Ginzburg_Landau_Tests_Lyapunov