program demo
   ! Standard Library
   use stdlib_optval, only : optval 
   use stdlib_linalg, only : eye, svdvals, eye
   use stdlib_math, only : logspace
   use stdlib_logger, only : information_level, warning_level, debug_level, error_level, none_level
   ! LightKrylov for linear algebra
   use LightKrylov
   use LightKrylov_Logger
   use LightKrylov_ExpmLib
   use LightKrylov_Utils
   ! LightROM
   use LightROM_AbstractLTIsystems
   use LightROM_Utils
   use LightROM_Timing
   use LightROM_DLRAIntegrators
   ! Laplacian
   use Laplacian2D_LTI_Lyapunov_Base
   use Laplacian2D_LTI_Lyapunov_Operators
   use Laplacian2D_LTI_Lyapunov_RKlib
   use Laplacian2D_LTI_Lyapunov_Utils
   implicit none

   character(len=*), parameter :: this_module = 'Laplacian2D_LTI_Lyapunov_Main'
   character(len=*), parameter :: home = 'example/DLRA_laplacian2D_lti_lyapunov/local/'

   !----------------------------------------------------------
   !-----     LYAPUNOV EQUATION FOR LAPLACE OPERATOR     -----
   !----------------------------------------------------------

   ! DLRA
   integer, parameter :: rkmax = 11
   integer, parameter :: rk_X0 = 2
   ! rk_B is set in laplacian2D.f90
   type(dlra_opts) :: opts

   integer  :: rk, torder
   real(dp) :: dt, tol, Tend, Tstep
   ! vector of dt values
   real(dp), allocatable :: dtv(:)
   ! vector of rank values
   integer,  allocatable :: rkv(:)
   ! vector of temporal order
   integer, allocatable  :: TOv(:)
   ! vector of tolerances
   real(dp), allocatable :: tolv(:)

   ! Exponential propagator (RKlib).
   type(rklib_lyapunov_mat), allocatable :: RK_propagator

   ! LTI system
   type(lti_system)                      :: LTI

   ! Laplacian
   type(laplace_operator),   allocatable :: A

   ! LR representation
   type(LR_state)                  :: X
   type(state_vector), allocatable :: U(:)
   real(dp),           allocatable :: S(:,:)
   
   ! STATE MATRIX (RKlib)
   type(state_matrix)              :: X_mat_RKlib(2)
   real(dp),           allocatable :: X_RKlib(:,:,:)
   real(dp)                        :: X_RKlib_ref(N,N)

   ! Initial condition
   type(state_vector)              :: U0(rkmax)
   real(dp)                        :: S0(rkmax,rkmax)
   ! Matrix
   real(dp)                        :: X0(N,N)

   ! OUTPUT
   real(dp)                        :: U_out(N,rkmax)
   real(dp)                        :: X_out(N,N)

   integer                         :: info, i, j, k, irep, nrep

   ! PROBLEM DEFINITION
   real(dp)  :: Adata(N,N)
   real(dp)  :: Bdata(N,N)
   real(dp)  :: BBTdata(N,N)
   
   ! LAPACK SOLUTION
   real(dp)  :: Xref(N,N)
   real(dp)  :: svals(N)

   ! Misc
   integer   :: clock_rate, clock_start, clock_stop
   integer, parameter :: irow = 8 ! how many numbers to print per row

   !--------------------------------
   ! Define which examples to run:
   !
   logical, parameter :: short_test = .true.
   !
   ! Skip the computations with small dt/small tolerance to speed up test
   !
   logical, parameter :: run_fixed_rank_short_integration_time_test   = .true.
   !
   ! Integrate the same initial condition for a short time with Runge-Kutta and DLRA.
   !
   ! The solution will be far from steady state (the residual will be large) for both methods.
   ! This test shows the convergence of the method as a function of the step size, the rank
   ! and the temporal order of DLRA.
   !
   logical, parameter :: run_fixed_rank_long_integration_time_test    = .true.
   !
   ! Integrate the same initial condition to steady state with Runge-Kutta and DLRA.
   !
   ! As the steady state is approached, the error/residual for Runge-Kutta goes to zero.
   ! Similarly, the test shows the effect of step size, rank and temporal order on the solution
   ! using DLRA
   !
   logical, parameter :: run_rank_adaptive_long_integration_time_test = .true.
   !
   ! Integrate the same initial condition to steady state with Runge-Kutta and DLRA using an 
   ! adaptive rank.
   !
   ! The DLRA algorthm automatically determines the rank necessary to integrate the equations
   ! such that the error on the singular values does not exceed a chosen tolerance. This rank
   ! depends on the tolerance but also the chosen time-step.
   !
   logical, parameter :: run_steady_state_convergence_test = .true.
   !
   ! Integrate the same initial condition to steady state with DLRA using the optimal rank.
   !
   ! Convergence is determined using the change in the resolved singular values, in particular
   ! the Frobenius norm of the low-rank solution. As expected, the minimum error scales with 
   ! the timestep and the temporal order of the scheme
   !
   !--------------------------------

   call logger_setup(logfile=trim(home)//'lightkrylov.log', log_level=error_level, log_stdout=.false., log_timestamp=.true.)

   ! Initialize timers for LightKrylov and LightROM
   call initialize_timers()
   call global_lightROM_timer%add_timer('DLRA Laplacian Lyapunov example', start=.true.)
   call global_lightROM_timer%add_timer('Direct solution (LAPACK)', start=.true.)
   ! Enumerate timers to check proper initialization
   call enumerate_timers()

   call system_clock(count_rate=clock_rate)

   print *, '#########################################################################'
   print *, '#                                                                       #'
   print *, '#               DYNAMIC LOW-RANK APPROXIMATION  -  DLRA                 #'
   print *, '#                                                                       #'
   print *, '#########################################################################'
   print *, ''
   print *, '          LYAPUNOV EQUATION FOR THE 2D LAPLACE OPERATOR:'
   print *, ''
   print *, '                   Algebraic Lyapunov equation:'
   print *, '                     0 = A @ X + X @ A.T + B @ B.T'
   print *, ''               
   print *, '                 Differential Lyapunov equation:'
   print *, '                   \dot{X} = A @ X + X @ A.T + B @ B.T'
   print *, ''
   print '(13X,A14,I4,"x",I4)', 'Problem size: ', N, N
   print *, ''
   print *, '            Initial condition: rank(X0) =', rk_X0
   print *, '            Inhomogeneity:     rank(B)  =', rk_B
   print *, ''
   print *, '#########################################################################'
   print *, ''

   ! Define RHS B
   do i = 1, rk_b
      call B(i)%rand(ifnorm = .false.)   
   end do
   call get_state(Bdata(:,:rk_b), B, 'Get B')
   BBTdata = matmul(Bdata(:,:rk_b), transpose(Bdata(:,:rk_b)))
   BBT(:N**2) = reshape(BBTdata, shape(BBT))

   ! Define LTI system
   LTI = lti_system()
   A = laplace_operator()
   call LTI%initialize_lti_system(A, B, B)
   call zero_basis(LTI%CT)

   ! Define initial condition of the form X0 + U0 @ S0 @ U0.T SPD
   print '(4X,A)', 'Define initial condition.'
   print *, ''
   call generate_random_initial_condition(U0(:rk_X0), S0(:rk_X0,:rk_X0), rk_X0)
   call get_state(U_out(:,:rk_X0), U0(:rk_X0), 'Get initial condition')
   
   ! Compute the full initial condition X0 = U_in @ S0 @ U_in.T
   X0 = matmul( U_out(:,:rk_X0), matmul(S0(:rk_X0,:rk_X0), transpose(U_out(:,:rk_X0))))
   
   print *, 'SVD X0'
   svals = svdvals(X0)
   print '(1X,A16,2X,*(F15.12,1X))', 'SVD(X0)[1-8]:', svals(:irow)

   !------------------
   ! COMPUTE EXACT SOLUTION OF THE LYAPUNOV EQUATION WITH LAPACK
   !------------------
   print *, ''
   print *, '#########################################################################'
   print *, '#                                                                       #'
   print *, '#    I.   Exact solution of the algebraic Lyapunov equation (LAPACK)    #'
   print *, '#                                                                       #'
   print *, '#########################################################################'
   print *, ''
   call system_clock(count=clock_start)     ! Start Timer
   ! Explicit 2D laplacian
   call build_operator(Adata)
   ! Solve Lyapunov equation
   call solve_lyapunov(Xref, Adata, BBTdata)
   call system_clock(count=clock_stop)      ! Stop Timer

   print '(A40,F10.4," s")', '--> X_ref.    Elapsed time:', real(clock_stop-clock_start)/real(clock_rate)
   print *, ''

   ! Explicit 2D laplacian
   call build_operator(Adata)
   ! sanity check
   X0 = CALE(Xref, Adata, BBTdata)
   print '(4X,A,E15.7)', 'Direct problem: | res(X_ref) |/N = ', norm2(X0)/N
   print *, ''
   ! compute svd
   svals = svdvals(Xref)
   print '(1X,A16,2X,*(F15.12,1X))', 'SVD(X_ref)[1-8]:', svals(:irow)

   call global_lightROM_timer%stop('Direct solution (LAPACK)')
   call global_lightROM_timer%add_timer('Short time: Runge-Kutta', start=.true.)

   !------------------
   ! COMPUTE SOLUTION WITH RK FOR DIFFERENT INTEGRATION TIMES AND COMPARE TO STUART-BARTELS
   !------------------
   print *, ''
   print *, '#########################################################################'
   print *, '#                                                                       #'
   print *, '#    IIa.  Solution using Runge-Kutta over a short time horizon         #'
   print *, '#                                                                       #'
   print *, '#########################################################################'
   print *, ''
   ! initialize exponential propagator
   Tend = 0.001_dp
   RK_propagator = rklib_lyapunov_mat(Tend)

   allocate(X_RKlib(N, N, 1))
   call get_state(U_out(:,:rk_X0), U0(:rk_X0), 'Get initial condition')
   X0 = matmul( U_out(:,:rk_X0), matmul(S0(:rk_X0,:rk_X0), transpose(U_out(:,:rk_X0))))
   call set_state(X_mat_RKlib(1:1), X0, 'Set RK X0')
   print '(A10,A26,A26,A20)', 'RKlib:','Tend','| X_RK - X_ref |/N', 'Elapsed time'
   print *, '         ------------------------------------------------------------------------'
   call system_clock(count=clock_start)     ! Start Timer
   ! integrate
   call RK_propagator%matvec(X_mat_RKlib(1), X_mat_RKlib(2))
   ! recover output
   call get_state(X_RKlib(:,:,1), X_mat_RKlib(2:2), 'Get RK solution')
   ! replace input
   call set_state(X_mat_RKlib(1:1), X_RKlib(:,:,1), 'Reset RK X0')
   call system_clock(count=clock_stop)      ! Stop Timer
   print '(I10,F26.4,E26.8,F18.4," s")', 1, Tend, &
                  & norm2(X_RKlib(:,:,1) - Xref)/N, &
                  & real(clock_stop-clock_start)/real(clock_rate)

   ! Choose relevant reference case from RKlib
   X_RKlib_ref = X_RKlib(:,:,1)

   deallocate(X_RKlib)

   print *, ''
   svals = svdvals(X_RKlib_ref)
   print '(1X,A16,2X,*(F15.12,1X))', 'SVD(X_RK )[1-8]:', svals(:irow)

   call global_lightROM_timer%stop('Short time: Runge-Kutta')
   call global_lightROM_timer%add_timer('Short time: DLRA', start=.true.)

   !------------------
   ! COMPUTE DLRA FOR SHORTEST INTEGRATION TIMES FOR DIFFERENT DT AND COMPARE WITH RK SOLUTION
   !------------------
   print *, ''
   print *, '#########################################################################'
   print *, '#                                                                       #'
   print *, '#    IIb.  Solution using fixed-rank DLRA                               #'
   print *, '#                                                                       #'
   print *, '#########################################################################'
   print *, ''

   if (run_fixed_rank_short_integration_time_test) then
      call logger%log_message('#', module=this_module)
      call logger%log_message('# Fixed-rank short time-horizon integration', module=this_module)
      call logger%log_message('#', module=this_module)
      print '(A10,A8,A4,A10,A8,3(A20),A20)', 'DLRA:','rk',' TO','dt','Tend','| X_LR - X_RK |/N', &
         & '| X_LR - X_ref |/N','| res_LR |/N', 'Elapsed time'
      write(*,'(A)', ADVANCE='NO') '         ------------------------------------------------'
      print '(A)', '--------------------------------------------------------------'
      
      ! Choose input ranks and integration steps
      rkv = [ 2, 3, 4 ]
      if (short_test) then
         dtv = logspace(-4.0_dp, -3.0_dp, 2, 10)
      else
         dtv = logspace(-6.0_dp, -3.0_dp, 4, 10)
      end if
      dtv = dtv(size(dtv):1:-1)
      TOv  = [ 1, 2 ]

      irep = 0
      X = LR_state()
      do i = 1, size(rkv)
         rk = rkv(i)
         do j = 1, size(Tov)
            torder = TOv(j)
            do k = 1, size(dtv)
               irep = irep + 1
               dt = dtv(k)

               ! Reset input
               call X%init(U0, S0, rk)

               ! run step
               opts = dlra_opts(mode=torder, if_rank_adaptive=.false.)
               call system_clock(count=clock_start)     ! Start Timer
               call Lyapunov_integrator(X, LTI%A, LTI%B, &
                                       & Tend, dt, info, &
                                       & exptA=exptA, options=opts)
               call system_clock(count=clock_stop)      ! Stop Timer

               ! Reconstruct solution
               call get_state(U_out(:,:rk), X%U, 'Reconstruct solution')
               X_out = matmul(U_out(:,:rk), matmul(X%S, transpose(U_out(:,:rk))))
               X0 = CALE(X_out, Adata, BBTdata)
               print '(A10,I8," TO",I1,F10.6,F8.4,3(E20.8),F18.4," s")', &
                                 & 'OUTPUT', rk, torder, dt, Tend, &
                                 & norm2(X_RKlib_ref - X_out)/N, norm2(X_out - Xref)/N, &
                                 & norm2(X0)/N, real(clock_stop-clock_start)/real(clock_rate)
               deallocate(X%U, X%S)
            end do
            print *, ''
         end do
         print *, ''
      end do
      nrep = irep

      svals = svdvals(X_RKlib_ref)
      print '(1X,A16,2X,*(F15.12,1X))', 'SVD(X_RK)[1-8]:', svals(:irow)
      svals = svdvals(X_out)
      print '(1X,A16,2X,*(F15.12,1X))', 'SVD(X_LR)[1-8]:', svals(:irow)

      call reset_solver(LyapunovSolver_counter, 'Lyapunov', 'Lyapunov_')
   else
      print *, 'Skip.'
   end if

   ! Reset timers
   call global_lightROM_timer%stop('Short time: DLRA')
   call A%reset_timer()
   call RK_propagator%reset_timer()
   call reset_timers()
   call global_lightROM_timer%add_timer('Steady-State: Runge-Kutta', start=.true.)

   !------------------
   ! COMPUTE SOLUTION WITH RK FOR DIFFERENT INTEGRATION TIMES AND COMPARE TO STUART-BARTELS
   !------------------
   print *, ''
   print *, '#########################################################################'
   print *, '#                                                                       #'
   print *, '#    IIIa.  Solution using Runge-Kutta to steady state                  #'
   print *, '#                                                                       #'
   print *, '#########################################################################'
   print *, ''
   ! initialize exponential propagator
   nrep  = 10
   Tstep = 0.1_dp
   RK_propagator = rklib_lyapunov_mat(Tstep)
   Tend = nrep*Tstep

   allocate(X_RKlib(N, N, nrep))
   call get_state(U_out(:,:rk_X0), U0(:rk_X0), 'Get initial conditions')
   X0 = matmul( U_out(:,:rk_X0), matmul(S0(:rk_X0,:rk_X0), transpose(U_out(:,:rk_X0))))
   call set_state(X_mat_RKlib(1:1), X0, 'Set RK X0')
   print '(A10,A26,A26,A20)', 'RKlib:','Tend','| X_RK - X_ref |/N', 'Elapsed time'
   print *, '         ------------------------------------------------------------------------'
   do irep = 1, nrep
      call system_clock(count=clock_start)     ! Start Timer
      ! integrate
      call RK_propagator%matvec(X_mat_RKlib(1), X_mat_RKlib(2))
      ! recover output
      call get_state(X_RKlib(:,:,irep), X_mat_RKlib(2:2), 'Get RK solution')
      ! replace input
      call set_state(X_mat_RKlib(1:1), X_RKlib(:,:,irep), 'Reset RK X0')
      call system_clock(count=clock_stop)      ! Stop Timer
      print '(I10,F26.4,E26.8,F18.4," s")', irep, irep*Tstep, &
                     & norm2(X_RKlib(:,:,irep) - Xref)/N, &
                     & real(clock_stop-clock_start)/real(clock_rate)
   end do

   ! Choose relevant reference case from RKlib
   X_RKlib_ref = X_RKlib(:,:,nrep)

   deallocate(X_RKlib)

   print *, ''
   svals = svdvals(Xref)
   print '(1X,A16,2X,*(F15.12,1X))', 'SVD(X_ref)[1-8]:', svals(:irow)
   svals = svdvals(X_RKlib_ref)
   print '(1X,A16,2X,*(F15.12,1X))', 'SVD(X_RK )[1-8]:', svals(:irow)

   call global_lightROM_timer%stop('Steady-State: Runge-Kutta')
   call global_lightROM_timer%add_timer('Steady-State: DLRA', start=.true.)

   !------------------
   ! COMPUTE DLRA FOR LONGER INTEGRATION TIMES FOR DIFFERENT DT AND COMPARE WITH RK SOLUTION
   !------------------
   print *, ''
   print *, '#########################################################################'
   print *, '#                                                                       #'
   print *, '#    IIIb.  Solution using fixed-rank DLRA                              #'
   print *, '#                                                                       #'
   print *, '#########################################################################'
   print *, ''
   if (run_fixed_rank_long_integration_time_test) then
      call logger%log_message('#', module=this_module)
      call logger%log_message('# Fixed-rank long time-horizon integration', module=this_module)
      call logger%log_message('#', module=this_module)
      print '(A10,A8,A4,A10,A8,3(A20),A20)', 'DLRA:','rk',' TO','dt','Tend','| X_LR - X_RK |/N', &
         & '| X_LR - X_ref |/N','| res_LR |/N', 'Elapsed time'
      write(*,'(A)', ADVANCE='NO') '         ------------------------------------------------'
      print '(A)', '--------------------------------------------------------------'
      
      ! Choose input ranks and integration steps
      rkv = [ 4,  8 ]
      if (short_test) then
         dtv = logspace(-2.0_dp, -1.0_dp, 2, 10)
      else
         dtv = logspace(-4.0_dp, -1.0_dp, 4, 10)
      end if
      dtv = dtv(size(dtv):1:-1)
      TOv  = [ 1, 2 ]

      irep = 0
      X = LR_state()
      do i = 1, size(rkv)
         rk = rkv(i)
         do j = 1, size(Tov)
            torder = TOv(j)
            do k = 1, size(dtv)
               irep = irep + 1
               dt = dtv(k)

               ! Reset input
               call X%init(U0, S0, rk)

               ! run step
               opts = dlra_opts(mode=torder, if_rank_adaptive=.false.)
               call system_clock(count=clock_start)     ! Start Timer
               call Lyapunov_integrator(X, LTI%A, LTI%B, &
                                       & Tend, dt, info, &
                                       & exptA=exptA, options=opts)
               call system_clock(count=clock_stop)      ! Stop Timer

               ! Reconstruct solution
               call get_state(U_out(:,:rk), X%U, 'Reconstruct solution')
               X_out = matmul(U_out(:,:rk), matmul(X%S, transpose(U_out(:,:rk))))
               X0 = CALE(X_out, Adata, BBTdata)
               print '(A10,I8," TO",I1,F10.6,F8.4,3(E20.8),F18.4," s")', &
                                 & 'OUTPUT', rk, torder, dt, Tend, &
                                 & norm2(X_RKlib_ref - X_out)/N, norm2(X_out - Xref)/N, &
                                 & norm2(X0)/N, real(clock_stop-clock_start)/real(clock_rate)
               deallocate(X%U, X%S)
            end do
            print *, ''
         end do
         print *, ''
      end do

      svals = svdvals(Xref)
      print '(1X,A16,2X,*(F15.12,1X))', 'SVD(X_ref)[1-8]:', svals(:irow)
      svals = svdvals(X_RKlib_ref)
      print '(1X,A16,2X,*(F15.12,1X))', 'SVD(X_RK )[1-8]:', svals(:irow)
      svals = svdvals(X_out)
      print '(1X,A16,2X,*(F15.12,1X))', 'SVD(X_LR )[1-8]:', svals(:irow)
      print *, ''
      call reset_solver(LyapunovSolver_counter, 'Lyapunov', 'Lyapunov_')
   else
      print *, 'Skip.'
   end if

   call global_lightROM_timer%stop('Steady-State: DLRA')
   call global_lightROM_timer%add_timer('Steady-State: rank-adaptive DLRA', start=.true.)

   print *, '#########################################################################'
   print *, '#                                                                       #'
   print *, '#    IIIc.  Solution using rank-adaptive DLRA                           #'
   print *, '#                                                                       #'
   print *, '#########################################################################'
   print *, ''
   if (run_rank_adaptive_long_integration_time_test) then
      call logger%log_message('#', module=this_module)
      call logger%log_message('# Rank-adaptive long time-horizon integration', module=this_module)
      call logger%log_message('#', module=this_module)
      ! Choose input ranks and integration step
      rk = 8 ! This is for initialisation, but the algorithm will choose the appropriate rank automatically
      if (short_test) then
         dtv = logspace(-2.0_dp, -1.0_dp, 2, 10)
         tolv = logspace(-8.0_dp, -4.0_dp, 2, 10)
      else
         dtv = logspace(-4.0_dp, -1.0_dp, 4, 10)
         tolv = logspace(-12.0_dp, -4.0_dp, 3, 10)
      end if
      dtv = dtv(size(dtv):1:-1)
      tolv = tolv(size(tolv):1:-1)
      TOv  = [ 1, 2 ]

      irep = 0
      X = LR_state()
      do i = 1, size(tolv)
         tol = tolv(i)
         print '(A,E9.2)', ' SVD tol = ', tol
         print *, ''
         print '(A10,A8,A4,A10,A8,3(A20),A20)', 'DLRA:','rk_end',' TO','dt','Tend','| X_LR - X_RK |/N', &
            & '| X_LR - X_ref |/N','| res_LR |/N', 'Elapsed time'
         write(*,'(A)', ADVANCE='NO') '         ------------------------------------------------'
         print '(A)', '--------------------------------------------------------------'
         do j = 1, size(Tov)
            torder = TOv(j)
            do k = 1, size(dtv)
               irep = irep + 1
               dt = dtv(k)

               ! Reset input
               call X%init(U0, S0, rk, rkmax, if_rank_adaptive=.true.)

               ! run step
               opts = dlra_opts(mode=torder, if_rank_adaptive=.true., tol=tol)
               call system_clock(count=clock_start)     ! Start Timer
               call Lyapunov_integrator(X, LTI%A, LTI%B, &
                                       & Tend, dt, info, &
                                       & exptA=exptA, options=opts)
               call system_clock(count=clock_stop)      ! Stop Timer
               rk = X%rk

               ! Reconstruct solution
               call get_state(U_out(:,:rk), X%U(:rk), 'Reconstruct solution')
               X_out = matmul(U_out(:,:rk), matmul(X%S(:rk,:rk), transpose(U_out(:,:rk))))
               X0 = CALE(X_out, Adata, BBTdata)
               print '(A10,I8," TO",I1,F10.6,F8.4,3(E20.8),F18.4," s")', &
                                    & 'OUTPUT', X%rk, torder, dt, Tend, &
                                    & norm2(X_RKlib_ref - X_out)/N, norm2(X_out - Xref)/N, &
                                    & norm2(X0)/N, real(clock_stop-clock_start)/real(clock_rate)
               deallocate(X%U, X%S)
            end do
            print *, ''
         end do
         svals = svdvals(Xref)
         print '(1X,A16,2X,*(F15.12,1X))', 'SVD(X_ref)[1-8]:', svals(:irow)
         svals = svdvals(X_RKlib_ref)
         print '(1X,A16,2X,*(F15.12,1X))', 'SVD(X_RK )[1-8]:', svals(:irow)
         svals = svdvals(X_out)
         print '(1X,A16,2X,*(F15.12,1X))', 'SVD(X_LR )[1-8]:', svals(:irow)
         print *, ''
         print *, '#########################################################################'
         print *, ''
      end do
      call reset_solver(LyapunovSolver_counter, 'Lyapunov', 'Lyapunov_')
   else
      print *, 'Skip.'
   end if

   ! Reset timers
   call global_lightROM_timer%stop('Steady-State: rank-adaptive DLRA')
   call A%reset_timer()
   call RK_propagator%reset_timer()
   call reset_timers()
   call global_lightROM_timer%add_timer('Steady-State: opt. rank DLRA', start=.true.)

   print *, '#########################################################################'
   print *, '#                                                                       #'
   print *, '#    IVa.  DLRA with optimal fixed rank to steady state                 #'
   print *, '#                                                                       #'
   print *, '#########################################################################'
   print *, ''
   if (run_steady_state_convergence_test) then
      call logger%log_message('#', module=this_module)
      call logger%log_message('# Optimum rank long time-horizon integration', module=this_module)
      call logger%log_message('#', module=this_module)
      ! Choose input ranks and integration steps
      rk = 7
      dt = 1e-4_dp
      !torder = 2
      if (short_test) then
         tolv = logspace(-6.0_dp, -4.0_dp, 2, 10)
         TOv  = [ 1 ]
      else
         tolv = logspace(-8.0_dp, -4.0_dp, 3, 10)
         TOv  = [ 1, 2 ]
      end if
      tolv = tolv(size(tolv):1:-1)
      Tend = 1.0_dp

      do i = 1, size(tolv)
         tol = tolv(i)
         print '(A,E9.2)', ' Increment tol = ', tol
         print '(A10,A8,A4,A10,A8,3(A20),A20)', 'DLRA:','rk',' TO','dt','Tend','| X_LR - X_RK |/N', &
                        & '| X_LR - X_ref |/N','| res_LR |/N', 'Elapsed time'
         write(*,'(A)', ADVANCE='NO') '         ------------------------------------------------'
         print '(A)', '--------------------------------------------------------------'
         do j = 1, size(TOv)
            torder = TOv(j)
            ! set options
            opts = dlra_opts( mode             = torder, &
                              ! convergence check
                           &  chkctrl_time     = .true., &
                           &  chktime          = 0.01_dp, &
                           &  inc_tol          = tol, &
                              ! rank-adaptive settings     
                           &  if_rank_adaptive = .false., &
                           &  tol              = 1e-6_dp)

            ! Reset input
            call X%init(U0, S0, rk)

            ! run step
            call system_clock(count=clock_start)     ! Start Timer
            call Lyapunov_integrator(X, LTI%A, LTI%B, &
                                    & Tend, dt, info, &
                                    & exptA=exptA, options=opts)
            call system_clock(count=clock_stop)      ! Stop Timer

            ! Reconstruct solution
            call get_state(U_out(:,:rk), X%U, 'Reconstruct solution')
            X_out = matmul(U_out(:,:rk), matmul(X%S, transpose(U_out(:,:rk))))
            X0 = CALE(X_out, Adata, BBTdata)
            print '(A10,I8," TO",I1,F10.6,F8.4,3(E20.8),F18.4," s")', &
                              & 'OUTPUT', rk, torder, dt, Tend, &
                              & norm2(X_RKlib_ref - X_out)/N, norm2(X_out - Xref)/N, &
                              & norm2(X0)/N, real(clock_stop-clock_start)/real(clock_rate)
            deallocate(X%U, X%S)
         end do
         print *, ''
         if (X%is_converged) then
            print '(A,F6.2,A,I0,A)', 'Solution converged at t= ', X%time, ' after ', X%step, ' steps.'
         else
            print '(A,F6.2,A,I0,A)', 'Solution not converged at t= ', X%time, 'after ', X%step, ' steps.'
         end if
         print *, ''
         svals = svdvals(Xref)
         print '(1X,A16,2X,*(F15.12,1X))', 'SVD(X_ref)[1-8]:', svals(:irow)
         svals = svdvals(X_RKlib_ref)
         print '(1X,A16,2X,*(F15.12,1X))', 'SVD(X_RK )[1-8]:', svals(:irow)
         svals = svdvals(X_out)
         print '(1X,A16,2X,*(F15.12,1X))', 'SVD(X_LR )[1-8]:', svals(:irow)
         print *, ''
      end do
      call reset_solver(LyapunovSolver_counter, 'Lyapunov', 'Lyapunov_')
   else
      print *, 'Skip.'
      print *, ''
   end if

   ! Compute and print timer summary
   call finalize_timers()

end program demo
