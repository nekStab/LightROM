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
   use stdlib_linalg, only : eye, diag
   use stdlib_math, only : all_close, logspace
   use stdlib_io_npy, only : save_npy
   implicit none

   ! DLRA
   integer, parameter :: rkmax = 14
   integer, parameter :: rk_X0 = 14
   logical, parameter :: verb  = .true.
   logical, parameter :: save  = .true.
   character*128      :: oname
   ! rk_B & rk_C are set in ginzburg_landau_base.f90

   integer  :: nrk, ntau, rk,  torder
   real(kind=wp) :: tau, Tend, Ttot
   ! vector of dt values
   real(kind=wp), allocatable :: tauv(:)
   ! vector of rank values
   integer, allocatable :: rkv(:)

   ! Exponential propagator (RKlib).
   type(exponential_prop), allocatable       :: A

   ! LTI system
   type(lti_system)                          :: LTI

   ! LR representation
   type(LR_state)                            :: X
   
   ! Initial condition
   real(kind=wp)                             :: U0_in(2*nx, rkmax)
   real(kind=wp)                             :: S0(rkmax,rkmax)
   type(state_vector), allocatable           :: U0(:)
   
   ! OUTPUT
   real(kind=wp)                             :: U_out(2*nx,rkmax)
   real(kind=wp)                             :: X_out(2*nx,2*nx)

   ! Information flag.
   integer                                   :: info
   integer                                   :: i, j, k, irep, nrep, nsteps
!
   ! SVD
   real(kind=wp)  :: U_svd(2*nx,2*nx)
   real(kind=wp)  :: S_svd(rkmax)
   real(kind=wp)  :: V_svd(rkmax,rkmax)
!
   ! timer
   integer   :: clock_rate, clock_start, clock_stop
   real(kind=wp) :: out(1,1)

   call system_clock(count_rate=clock_rate)

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
   allocate(U0(1:rk_X0))
   call init_rand(U0, .false.)
   call get_state(U0_in(:,1:rk_X0), U0)
   ! Compute SVD to get low-rank representation
   call svd(U0_in(:,1:rk_X0), U_svd(:,1:2*nx), S_svd(1:rk_X0), V_svd(1:rk_X0,1:rk_X0))
   S0 = 0.0_wp
   do i = 1,rk_X0
      S0(i,i) = S_svd(i)
   end do
   call set_state(U0, U_svd(:,1:rk_X0))

   tau  = 0.1_wp
   Tend = 1.0_wp
   nrep = 60

   nrk  = 4; allocate(rkv(1:nrk));   rkv  = (/ 2, 6, 10, 14 /)
   ntau = 2; allocate(tauv(1:ntau)); tauv = (/ 0.1, 0.01/)

   write(*,*) '---------------------'
   write(*,*) '   CONTROLABILITY'
   write(*,*) '---------------------'
   do torder = 1, 2
      do i = 1, nrk
         rk = rkv(i)
         do j = 1, ntau
            tau = tauv(j)
            ! Initialize low-rank representation with rank rk
            if (verb) write(*,*) 'Initialize LR state, rk =', rk
            X = LR_state()
            call X%initialize_LR_state(U0, S0, rk)
            ! Reset time
            Ttot = 0.0_wp
            write(*,'(A16,A4,A4,A10,A6,A8,A16,A20)') 'DLRA:','  rk',' TO','dt','steps','Tend','|| X_DLRA ||_2', 'Elapsed time'
            nsteps = nint(Tend/tau)
            do irep = 1, nrep
               ! run integrator
               call system_clock(count=clock_start)     ! Start Timer
               call numerical_low_rank_splitting_lyapunov_integrator(X, LTI%A, LTI%B, Tend, tau, torder, info, &
                                                                  & exptA=exptA, iftrans=.false., ifverb=.false.)
               call system_clock(count=clock_stop)      ! Stop Timer

               ! Reconstruct solution
               call get_state(U_out(:,1:rk), X%U)
               X_out = matmul(U_out(:,1:rk), matmul(X%S, transpose(U_out(:,1:rk))))

               Ttot = Ttot + Tend
               write(*,'(I4," ",A11,I4," TO",I1,F10.6,I6,F8.4,E16.8,F18.4," s")') irep, 'Xctl OUTPUT', &
                                 & rk, torder, tau, nsteps, Ttot, &
                                 & norm2(X_out), real(clock_stop-clock_start)/real(clock_rate)
            end do
            if (rk == 14 .and. j == 3) then
               if (save) then
                  write(oname,'("example/DLRA_ginzburg_landau/GL_Xdata_TO",I1,"_rk",I2.2,"_t",I1,".npy")') torder, rk, j
                  call save_npy(oname, X_out)
               end if
            end if
            deallocate(X%U)
            deallocate(X%S)
         end do
      end do
   end do

   write(*,*) '--------------------'
   write(*,*) '   OBSERVABILITY'
   write(*,*) '--------------------'
   do torder = 1, 2
      do i = 1, nrk
         rk = rkv(i)
         do j = 1, ntau
            tau = tauv(j)
            ! Initialize low-rank representation with rank rk
            if (verb) write(*,*) 'Initialize LR state, rk =', rk
            X = LR_state()
            call X%initialize_LR_state(U0, S0, rk)
            ! Reset time
            Ttot = 0.0_wp
            write(*,'(A16,A4,A4,A10,A6,A8,A16,A20)') 'DLRA:','  rk',' TO','dt','steps','Tend','|| X_DLRA ||_2', 'Elapsed time'
            nsteps = nint(Tend/tau)
            do irep = 1, nrep
               ! run integrator
               call system_clock(count=clock_start)     ! Start Timer
               call numerical_low_rank_splitting_lyapunov_integrator(X, LTI%A, LTI%CT, Tend, tau, torder, info, &
                                                                  & exptA=exptA, iftrans=.true., ifverb=.false.)
               call system_clock(count=clock_stop)      ! Stop Timer

               ! Reconstruct solution
               call get_state(U_out(:,1:rk), X%U)
               X_out = matmul(U_out(:,1:rk), matmul(X%S, transpose(U_out(:,1:rk))))

               Ttot = Ttot + Tend
               write(*,'(I4," ",A11,I4," TO",I1,F10.6,I6,F8.4,E16.8,F18.4," s")') irep, 'Yobs OUTPUT', &
                                 & rk, torder, tau, nsteps, Ttot, &
                                 & norm2(X_out), real(clock_stop-clock_start)/real(clock_rate)
            end do
            if (rk == 14 .and. j == 3) then
               if (save) then
                  write(oname,'("example/DLRA_ginzburg_landau/GL_Ydata_TO",I1,"_rk",I2.2,"_t",I1,".npy")') torder, rk, j
                  call save_npy(oname, X_out)
               end if
            end if
            deallocate(X%U)
            deallocate(X%S)
         end do
      end do
   end do

   return
end program demo