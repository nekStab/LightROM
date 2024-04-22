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
   logical, parameter :: verb  = .false.
   logical, parameter :: save  = .false.
   character*128      :: oname
   ! rk_B & rk_C are set in ginzburg_landau_base.f90

   integer  :: nrk, ndt, rk,  torder
   real(wp) :: dt, Tend
   ! vector of dt values
   real(wp), allocatable :: dtv(:)
   ! vector of rank values
   integer,  allocatable :: rkv(:)

   ! Exponential propagator (RKlib).
   type(exponential_prop), allocatable :: prop

   ! LTI system
   type(lti_system)                    :: LTI
   integer                             :: p

   ! LR representation
   type(LR_state)                      :: X
   type(state_vector), allocatable     :: U(:)
   real(wp) , allocatable              :: S(:,:)
   
   ! Initial condition
   real(wp)                            :: U0(2*nx, rkmax)
   real(wp)                            :: S0(rkmax,rkmax)
   !real(wp)                            :: X0(nx,nx)
!   
   !! OUTPUT
   !real(wp)                        :: U_out(N,rkmax)
   !real(wp)                        :: X_out(N,N)
!
   !!> Information flag.
   integer                             :: info
   integer                             :: i, j, k, irep, nrep
!
   !! PROBLEM DEFINITION
   !real(wp)  :: Adata(N,N)
   !real(wp)  :: Bdata(N,rkmax)
   !real(wp)  :: BBTdata(N,N)
   !
   !! LAPACK SOLUTION
   !real(wp)  :: Xref(N,N)
   !! DSYTD2
   !real(wp)  :: Dm(N), work(N), wr(N), wi(N)
   !real(wp)  :: E(N-1), tw(N-1)
   !real(wp)  :: T(N,N), Q(N,N), Z(N,N), Vdata(N,N), Wdata(N,N), Ydata(N,N)
   !real(wp)  :: scale
   !integer   :: isgn
   !! SVD
   real(wp)  :: U_svd(2*nx,2*nx)
   real(wp)  :: S_svd(rkmax)
   real(wp)  :: V_svd(rkmax,rkmax)
!
   ! timer
   integer   :: clock_rate, clock_start, clock_stop

   !----------------------------------
   !-----     INITIALIZATION     -----
   !----------------------------------

   ! Initialize parameters (B, CT)
   call initialize_parameters()

   ! Define initial condition
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
   rk = 6
   X = LR_state()
   call X%initialize_LR_state(U,S0,rk)

   !LTI = lti_system(prop)

   ! Initialize propagator.
   !A = exponential_prop(tau)



end program demo