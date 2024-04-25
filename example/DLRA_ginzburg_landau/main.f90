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
   use stdlib_io_npy, only : save_npy, load_npy

   implicit none

   ! DLRA
   integer, parameter :: rkmax = 14
   integer, parameter :: rk_X0 = 14
   logical, parameter :: verb  = .true.
   logical, parameter :: run_DLRA_test = .false.
   logical, parameter :: run_BT_test = .true.
   logical, parameter :: save  = .true.
   logical, parameter :: load  = .true.
   character*128      :: oname
   character*128      :: onameU
   character*128      :: onameS
   character*128, parameter :: basepath = 'local/'
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
   type(LR_state)                            :: X     ! Controllability
   type(LR_state)                            :: Y     ! Observability

   ! BT
   type(state_vector), allocatable           :: T(:)
   type(state_vector), allocatable           :: Tinv(:)
   real(kind=wp),      allocatable           :: S(:)
   real(kind=wp),      allocatable           :: U_load(:,:)
   real(kind=wp),      allocatable           :: S_load(:,:)

   ! ROM
   real(kind=wp),      allocatable           :: Ahat(:,:)
   real(kind=wp),      allocatable           :: Bhat(:,:)
   real(kind=wp),      allocatable           :: Chat(:,:)
   real(kind=wp),      allocatable           :: D(:,:)
   
   ! Initial condition
   real(kind=wp)                             :: U0_in(2*nx, rkmax)
   real(kind=wp)                             :: S0(rkmax,rkmax)
   type(state_vector), allocatable           :: U0(:)
   type(state_vector), allocatable           :: Utmp(:)
   
   ! OUTPUT
   real(kind=wp)                             :: U_out(2*nx,rkmax)
   real(kind=wp)                             :: X_out(2*nx,2*nx)

   ! Information flag.
   integer                                   :: info
   integer                                   :: i, j, k, irep, nrep, nsteps
   logical                                   :: existU, existS
!
   ! SVD
   real(kind=wp)  :: U_svd(2*nx,2*nx)
   real(kind=wp)  :: S_svd(rkmax)
   real(kind=wp)  :: V_svd(rkmax,rkmax)
!
   ! timer
   integer   :: clock_rate, clock_start, clock_stop
   real(kind=wp) :: Tmax

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

   if (run_DLRA_test) then
      Tend = 1.0_wp
      nrep = 50

      nrk  = 4; allocate(rkv(1:nrk));   rkv  = (/ 2, 6, 10, 14 /)
      ntau = 2; allocate(tauv(1:ntau)); tauv = (/ 0.1, 0.01 /)

      write(*,*) '----------------------'
      write(*,*) '   CONTROLLABILITY'
      write(*,*) '----------------------'
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
               if (save) then
                  write(oname,'("GL_Xdata_TO",I1,"_rk",I2.2,"_t",I1,".npy")') torder, rk, j
                  call save_npy(trim(basepath)//trim(oname), X_out)
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
               Y = LR_state()
               call Y%initialize_LR_state(U0, S0, rk)
               ! Reset time
               Ttot = 0.0_wp
               write(*,'(A16,A4,A4,A10,A6,A8,A16,A20)') 'DLRA:','  rk',' TO','dt','steps','Tend','|| X_DLRA ||_2', 'Elapsed time'
               nsteps = nint(Tend/tau)
               do irep = 1, nrep
                  ! run integrator
                  call system_clock(count=clock_start)     ! Start Timer
                  call numerical_low_rank_splitting_lyapunov_integrator(Y, LTI%A, LTI%CT, Tend, tau, torder, info, &
                                                                     & exptA=exptA, iftrans=.true., ifverb=.false.)
                  call system_clock(count=clock_stop)      ! Stop Timer

                  ! Reconstruct solution
                  call get_state(U_out(:,1:rk), Y%U)
                  X_out = matmul(U_out(:,1:rk), matmul(Y%S, transpose(U_out(:,1:rk))))

                  Ttot = Ttot + Tend
                  write(*,'(I4," ",A11,I4," TO",I1,F10.6,I6,F8.4,E16.8,F18.4," s")') irep, 'Yobs OUTPUT', &
                                    & rk, torder, tau, nsteps, Ttot, &
                                    & norm2(X_out), real(clock_stop-clock_start)/real(clock_rate)
               end do
               if (save) then
                  write(oname,'("GL_Ydata_TO",I1,"_rk",I2.2,"_t",I1,".npy")') torder, rk, j
                  call save_npy(trim(basepath)//trim(oname), X_out)
               end if
               deallocate(Y%U)
               deallocate(Y%S)
            end do
         end do
      end do
   end if

   if (run_BT_test) then
      ! Set parameters
      rk     = 12
      tau    = 0.1_wp
      torder = 2
      ! integration time
      Tmax   = 50.0_wp
      nrep   = 5
      Tend   = Tmax/nrep
      nsteps = nint(Tend/tau)

      write(*,*)
      write(*,*) '----------------------'
      write(*,*) '   CONTROLLABILITY'
      write(*,*) '----------------------'
      write(*,*)

      onameU = "GL_Xctl_U.npy"
      onameS = "GL_Xctl_S.npy"

      X = LR_state()
      if (load) then
         inquire(file=trim(onameU), exist=existU)
         inquire(file=trim(onameS), exist=existS)
         if (existU .and. existS) then
            write(*,*) 'Load data from file:'
            write(*,*) '    ', trim(basepath)//trim(onameU)
            write(*,*) '    ', trim(basepath)//trim(onameS)
            call load_npy(trim(basepath)//trim(onameU), U_load)
            call load_npy(trim(basepath)//trim(onameU), S_load)
         else
            write(*,*) 'Files to load X not found.'
            STOP 1
         end if
         if (.not.allocated(Utmp)) allocate(Utmp(1:rk), source=U0(1))
         call set_state(Utmp, U_load(:,1:rk))
         call X%initialize_LR_state(Utmp, S_load, rk)
      else
         ! Initialize low-rank representation with rank rk
         if (verb) write(*,*) 'Initialize LR state, rk =', rk
         call X%initialize_LR_state(U0, S0, rk)
         ! Reset time
         Ttot = 0.0_wp
         write(*,'(A16,A4,A4,A10,A6,A8,A16,A20)') 'DLRA:','  rk',' TO','dt','steps','Tend','|| X_DLRA ||_2', 'Elapsed time'
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
         if (save) then
            call save_npy(onameU, U_out(:,1:rk))
            call save_npy(onameS, X%S(1:rk,1:rk))
         end if
      end if

      write(*,*)
      write(*,*) '--------------------'
      write(*,*) '   OBSERVABILITY'
      write(*,*) '--------------------'
      write(*,*)

      onameU = "GL_Yobs_U.npy"
      onameS = "GL_Yobs_S.npy"

      Y = LR_state()
      if (load) then
         inquire(file=trim(onameU), exist=existU)
         inquire(file=trim(onameS), exist=existS)
         if (existU .and. existS) then
            write(*,*) 'Load data from file:'
            write(*,*) '    ', trim(basepath)//trim(onameU)
            write(*,*) '    ', trim(basepath)//trim(onameS)
            call load_npy(trim(basepath)//trim(onameU), U_load)
            call load_npy(trim(basepath)//trim(onameU), S_load)
         else
            write(*,*) 'Files to load Y not found.'
            STOP 1
         end if
         if (.not.allocated(Utmp)) allocate(Utmp(1:rk), source=U0(1))
         call set_state(Utmp, U_load(:,1:rk))
         call Y%initialize_LR_state(Utmp, S_load, rk)
      else
         ! Initialize low-rank representation with rank rk
         if (verb) write(*,*) 'Initialize LR state, rk =', rk
         call Y%initialize_LR_state(U0, S0, rk)
         ! Reset time
         Ttot = 0.0_wp
         write(*,'(A16,A4,A4,A10,A6,A8,A16,A20)') 'DLRA:','  rk',' TO','dt','steps','Tend','|| X_DLRA ||_2', 'Elapsed time'
         do irep = 1, nrep
            ! run integrator
            call system_clock(count=clock_start)     ! Start Timer
            call numerical_low_rank_splitting_lyapunov_integrator(Y, LTI%A, LTI%CT, Tend, tau, torder, info, &
                                                               & exptA=exptA, iftrans=.true., ifverb=.false.)
            call system_clock(count=clock_stop)      ! Stop Timer

            ! Reconstruct solution
            call get_state(U_out(:,1:rk), Y%U)
            X_out = matmul(U_out(:,1:rk), matmul(Y%S, transpose(U_out(:,1:rk))))

            Ttot = Ttot + Tend
            write(*,'(I4," ",A11,I4," TO",I1,F10.6,I6,F8.4,E16.8,F18.4," s")') irep, 'Yobs OUTPUT', &
                              & rk, torder, tau, nsteps, Ttot, &
                              & norm2(X_out), real(clock_stop-clock_start)/real(clock_rate)
         end do
         if (save) then
            call save_npy(trim(basepath)//trim(onameU), U_out(:,1:rk))
            call save_npy(trim(basepath)//trim(onameS), Y%S(1:rk,1:rk))
         end if
      end if 
   
      write(*,*)
      write(*,*) '------------------------------'
      write(*,*) '   BALANCING TRANSFORMATION'
      write(*,*) '------------------------------'
      write(*,*)

      allocate(T(1:rk), source=U0(1))
      allocate(Tinv(1:rk), source=U0(1))
      allocate(S(1:rk))
      call Balanced_Transformation(T, S, Tinv, X, Y)
      call ROM_Petrov_Galerkin_Projection(Ahat, Bhat, Chat, D, LTI, T, Tinv)
      if (save) then
         write(onameU,'("example/DLRA_ginzburg_landau/GL_Ahat.npy")')
         call save_npy(onameU, Ahat)
         write(onameU,'("example/DLRA_ginzburg_landau/GL_Bhat.npy")')
         call save_npy(onameU, Bhat)
         write(onameU,'("example/DLRA_ginzburg_landau/GL_Chat.npy")')
         call save_npy(onameU, Chat)
      end if
   end if

   return
end program demo