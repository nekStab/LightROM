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
   logical, parameter :: run_DLRA_test = .true.
   logical, parameter :: run_BT_test = .false.
   logical, parameter :: save  = .false.
   logical, parameter :: load  = .false.
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
   type(GL_operator),      allocatable       :: A
   type(exponential_prop), allocatable       :: prop

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
   real(kind=wp),      allocatable           :: vecs(:,:)
   complex(kind=wp),   allocatable           :: vals(:)

   ! ROM
   real(kind=wp),      allocatable           :: Swrk(:,:)
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

   ! Counters
   integer                                   :: i, j, k, irep, nrep, nsteps

   ! IO
   integer                                   :: iostatus
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

   ! Initialize mesh and system parameters A, B, CT
   if (verb) write(*,*) 'Initialize parameters'
   call initialize_parameters()

   ! Initialize propagator
   if (verb) write(*,*) 'Initialize exponential propagator'
   prop = exponential_prop(1.0_wp)

   ! Initialize LTI system
   A = GL_operator()
   if (verb) write(*,*) 'Initialize LTI system (A, prop, B, CT, _)'
   LTI = lti_system()
   call LTI%initialize_lti_system(A, prop, B, CT)

   ! Define initial condition
   if (verb) write(*,*) 'Define initial condition'
   allocate(U0(1:rk_X0))
   call init_rand(U0, .false.)
   call get_state(U0_in(:,1:rk_X0), U0)
   ! Compute SVD to get low-rank representation
   call svd(U0_in(:,1:rk_X0), U_svd(:,1:2*nx), S_svd(1:rk_X0), V_svd(1:rk_X0,1:rk_X0))
   S0 = 0.0_wp
   S0(1:rk_X0,1:rk_X0) = diag(S_svd(1:rk_X0))
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
                  call numerical_low_rank_splitting_lyapunov_integrator(X, LTI%prop, LTI%B, Tend, tau, torder, info, &
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
                  call save_npy(oname, X_out)
                  if (iostatus /= 0) then; write(*,*) "Error loading file", trim(oname); STOP 2; end if
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
                  call numerical_low_rank_splitting_lyapunov_integrator(Y, LTI%prop, LTI%CT, Tend, tau, torder, info, &
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
                  call save_npy(oname, X_out)
                  if (iostatus /= 0) then; write(*,*) "Error loading file", trim(oname); STOP 2; end if
               end if
               deallocate(Y%U)
               deallocate(Y%S)
            end do
         end do
      end do
   end if

   if (run_BT_test) then
      ! Set parameters
      rk     = 14
      tau    = 0.1_wp
      torder = 2
      ! integration time
      Tmax   = 50.0_wp
      nrep   = 10
      Tend   = Tmax/nrep
      nsteps = nint(Tend/tau)
      allocate(vals(1:rk))
      allocate(vecs(1:rk,1:rk))

      write(*,*) ''
      write(*,*) '----------------------'
      write(*,*) '   CONTROLLABILITY'
      write(*,*) '----------------------'
      write(*,*) ''

      onameU = trim(basepath)//"GL_Xctl_U.npy"
      onameS = trim(basepath)//"GL_Xctl_S.npy"

      X = LR_state()
      if (load) then
         write(*,*) 'Load data from file:'
         write(*,*) '    ', trim(onameU)
         write(*,*) '    ', trim(onameS)
         inquire(file=onameU, exist=existU)
         inquire(file=onameS, exist=existS)
         if (existU .and. existS) then
            call load_npy(onameU, U_load, iostatus)
            if (iostatus /= 0) then; write(*,*) "Error loading file", trim(onameU); STOP 2; end if
            call load_npy(onameS, S_load, iostatus)
            if (iostatus /= 0) then; write(*,*) "Error loading file", trim(onameS); STOP 2; end if
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
            call numerical_low_rank_splitting_lyapunov_integrator(X, LTI%prop, LTI%B, Tend, tau, torder, info, &
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
            write(*,*) 'Save data to file:'
            write(*,*) '    ', trim(onameU)
            write(*,*) '    ', trim(onameS)
            call save_npy(onameU, U_out(:,1:rk), iostatus)
            if (iostatus /= 0) then; write(*,*) "Error saving file", trim(onameU); STOP 2; end if
            call save_npy(onameS, X%S(1:rk,1:rk), iostatus)
            if (iostatus /= 0) then; write(*,*) "Error saving file", trim(onameS); STOP 2; end if
         end if
      end if

      write(*,*) ''
      write(*,*) '--------------------'
      write(*,*) '   OBSERVABILITY'
      write(*,*) '--------------------'
      write(*,*) ''

      onameU = trim(basepath)//"GL_Yobs_U.npy"
      onameS = trim(basepath)//"GL_Yobs_S.npy"

      Y = LR_state()
      if (load) then
         write(*,*) 'Load data from file:'
         write(*,*) '    ', trim(onameU)
         write(*,*) '    ', trim(onameS)
         inquire(file=onameU, exist=existU)
         inquire(file=onameS, exist=existS)
         if (existU .and. existS) then
            call load_npy(onameU, U_load, iostatus)
            if (iostatus /= 0) then; write(*,*) "Error loading file", trim(onameU); STOP 2; end if
            call load_npy(onameS, S_load, iostatus)
            if (iostatus /= 0) then; write(*,*) "Error loading file", trim(onameS); STOP 2; end if
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
            call numerical_low_rank_splitting_lyapunov_integrator(Y, LTI%prop, LTI%CT, Tend, tau, torder, info, &
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
            write(*,*) 'Save data to file:'
            write(*,*) '    ', trim(onameU)
            write(*,*) '    ', trim(onameS)
            call save_npy(onameU, U_out(:,1:rk), iostatus)
            if (iostatus /= 0) then; write(*,*) "Error saving file", trim(onameU); STOP 2; end if
            call save_npy(onameS, Y%S(1:rk,1:rk), iostatus)
            if (iostatus /= 0) then; write(*,*) "Error saving file", trim(onameS); STOP 2; end if
         end if
      end if
         
      write(*,*) ''
      write(*,*) '------------------------------'
      write(*,*) '   BALANCING TRANSFORMATION'
      write(*,*) '------------------------------'
      write(*,*) ''

      allocate(Swrk(1:rk,1:rk))
      if (.not.allocated(Utmp)) allocate(Utmp(1:rk), source=U0(1))

      ! compute sqrt of coefficient matrix X%S and right-multiply it to X%U
      Swrk = 0.0_wp
      call sqrtm(Swrk(1:rk,1:rk), X%S)
      call mat_mult(Utmp, X%U, Swrk(1:rk,1:rk))
      call get_state(U0_in(:,1:rk), Utmp)
      ! compute SVD of updated X%U
      call svd(U0_in(:,1:rk), U_svd(:,1:2*nx), S_svd(1:rk), V_svd(1:rk,1:rk))
      call set_state(X%U, matmul(U_svd(:,1:rk), diag(S_svd(1:rk))))

      ! compute sqrt of coefficient matrix Y%S and right-multiply it to Y%U
      Swrk = 0.0_wp
      call sqrtm(Swrk(1:rk,1:rk), Y%S)
      call mat_mult(Utmp, Y%U, Swrk(1:rk,1:rk))
      call get_state(U0_in(:,1:rk), Utmp)
      ! compute SVD of updated Y%U
      call svd(U0_in(:,1:rk), U_svd(:,1:2*nx), S_svd(1:rk), V_svd(1:rk,1:rk)) 
      call set_state(Y%U, matmul(U_svd(:,1:rk), diag(S_svd(1:rk))))   

      ! compute balancing transformation based on SVD of Gramians
      allocate(T(1:rk), source=U0(1)); allocate(Tinv(1:rk), source=U0(1)); allocate(S(1:rk))
      call Balancing_Transformation(T, S, Tinv, X%U, Y%U)
      
      call ROM_Petrov_Galerkin_Projection(Ahat, Bhat, Chat, D, LTI, T, Tinv)

      if (save) then
         write(*,*) 'Save data to file:'
         onameU = trim(basepath)//"GL_Ahat.npy"
         write(*,*) '    ', trim(onameU)
         call save_npy(onameU, Ahat)
         if (iostatus /= 0) then; write(*,*) "Error saving file", trim(onameU); STOP 2; end if
         onameU = trim(basepath)//"GL_Bhat.npy"
         write(*,*) '    ', trim(onameU)
         call save_npy(onameU, Bhat)
         if (iostatus /= 0) then; write(*,*) "Error saving file", trim(onameU); STOP 2; end if
         onameU = trim(basepath)//"GL_Chat.npy"
         write(*,*) '    ', trim(onameU)
         call save_npy(onameU, Chat)
         if (iostatus /= 0) then; write(*,*) "Error saving file", trim(onameU); STOP 2; end if
      end if
      write(*,*) ''
   end if

   return
end program demo