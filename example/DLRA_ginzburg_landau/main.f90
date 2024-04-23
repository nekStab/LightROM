program demo
   use LightKrylov
   use LightKrylov_expmlib
   use LightKrylov_utils

   use LightROM_AbstractLTIsystems
   use LightROM_utils

   use LightROM_LyapunovSolvers
   use LightROM_LyapunovUtils

   use Ginzburg_Landau_Base
   use Ginzburg_Landau_RKlib

   use stdlib_optval, only : optval 
   use stdlib_linalg, only : eye
   use stdlib_math, only : all_close, logspace
   use stdlib_io_npy, only : save_npy
   implicit none

   ! DLRA
   integer, parameter :: rkmax = 14
   integer, parameter :: rk_X0 = 10
   logical, parameter :: verb  = .true.
   logical, parameter :: save  = .false.
   character*128      :: oname
   ! rk_B & rk_C are set in ginzburg_landau_base.f90

   integer  :: nrk, ntau, rk,  torder
   real(kind=wp) :: tau, Tend
   ! vector of dt values
   real(kind=wp), allocatable :: dtv(:)
   ! vector of rank values
   integer, allocatable :: rkv(:)

   ! Exponential propagator (RKlib).
   type(exponential_prop), allocatable       :: A

   ! LTI system
   type(lti_system)                          :: LTI

   ! LR representation
   type(LR_state)                            :: X
   type(state_vector), allocatable           :: U(:)
   real(kind=wp) , allocatable               :: S(:,:)
   
   ! Initial condition
   real(kind=wp)                             :: U0(2*nx, rkmax)
   real(kind=wp)                             :: S0(rkmax,rkmax)
   !real(kind=wp)                            :: X0(nx,nx)
!   
   !! OUTPUT
   !real(kind=wp)                        :: U_out(N,rkmax)
   !real(kind=wp)                        :: X_out(N,N)
!
   !!> Information flag.
   integer                                   :: info
   integer                                   :: i, j, k, irep, nrep
!
   !! PROBLEM DEFINITION
   !real(kind=wp)  :: Adata(N,N)
   !real(kind=wp)  :: Bdata(N,rkmax)
   !real(kind=wp)  :: BBTdata(N,N)
   !
   !! LAPACK SOLUTION
   !real(kind=wp)  :: Xref(N,N)
   !! DSYTD2
   !real(kind=wp)  :: Dm(N), work(N), wr(N), wi(N)
   !real(kind=wp)  :: E(N-1), tw(N-1)
   !real(kind=wp)  :: T(N,N), Q(N,N), Z(N,N), Vdata(N,N), Wdata(N,N), Ydata(N,N)
   !real(kind=wp)  :: scale
   !integer   :: isgn
   !! SVD
   real(kind=wp)  :: U_svd(2*nx,2*nx)
   real(kind=wp)  :: S_svd(rkmax)
   real(kind=wp)  :: V_svd(rkmax,rkmax)
!
   ! timer
   integer   :: clock_rate, clock_start, clock_stop

   !----------------------------------
   !-----     INITIALIZATION     -----
   !----------------------------------

   ! Initialize mesh and system parameters B, CT
   if (verb) write(*,*) 'Initialize parameters'
   call initialize_parameters()

   ! Initialize propagator
   if (verb) write(*,*) 'Initialize exponential propagator'
   A = exponential_prop(1.0_wp)

   ! Initialize LTI system
   if (verb) write(*,*) 'Initialize LTI system (A, B, CT, _)'
   LTI = lti_system()
   call LTI%initialize_lti_system(A, B, CT)

   ! Define initial condition
   if (verb) write(*,*) 'Define initial condition'
   allocate(U(1:rk_X0))
   call init_rand(U, .false.)
   call get_state(U0(:,1:rk_X0), U)
   ! Compute SVD to get low-rank representation
   call svd(U0(:,1:rk_X0), U_svd(:,1:2*nx), S_svd(1:rk_X0), V_svd(1:rk_X0,1:rk_X0))
   S0 = 0.0_wp
   do i = 1,rk_X0
      S0(i,i) = S_svd(i)
   end do
   call set_state(U, U_svd(:,1:rk_X0))

   ! Initialize low-rank representation with rank rk
   if (verb) write(*,*) 'Initialize LR state'
   rk = 6
   X = LR_state()
   call X%initialize_LR_state(U, S0, rk)

   tau  = 0.05_wp
   Tend = 1.0_wp

   ! run integrator
   if (verb) write(*,*) 'Integrate ... t0 =', 0.0_wp, ', dt = ', tau
   call system_clock(count=clock_start)     ! Start Timer
   call numerical_low_rank_splitting_lyapunov_integrator(X, LTI%A, LTI%B, Tend, tau, 1, info, &
                                                       & exptA=exptA, iftrans=.false., ifverb=verb)
   call system_clock(count=clock_stop)      ! Stop Timer
   if (verb) write(*,*) 'done. Tend =', Tend

   return
end program demo