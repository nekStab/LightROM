module Ginzburg_Landau_Utils
   ! LightKrylov for linear algebra.
   use LightKrylov
   use LightKrylov_utils
   ! LightROM
   use LightROM_AbstractLTIsystems
   use LightROM_utils
   ! Ginzburg Landau
   use Ginzburg_Landau_Base
   use Ginzburg_Landau_Operators

   use LightROM_LyapunovSolvers
   use LightROM_LyapunovUtils

   ! Standard Library.
   use stdlib_optval, only : optval
   use stdlib_linalg, only : eye, diag
   use stdlib_io_npy, only : save_npy, load_npy
   implicit none

   ! IO
   integer,       parameter :: iunit1 = 1
   integer,       parameter :: iunit2 = 2
   character*128, parameter :: basepath = 'local/'
   integer,       parameter :: rkmax = 14
   integer,       parameter :: rk_X0 = 14

   private
   public :: iunit1, iunit2, basepath, rkmax, rk_X0
   public :: run_DLRA_test, run_BT_test
   public :: stamp_logfile_header

contains

   subroutine run_DLRA_test(LTI, U0, S0, rkv, tauv, Tend, nrep, ifsave, ifverb, iflogs)
      ! LTI system
      type(lti_system),              intent(inout) :: LTI
      ! Initial condition
      type(state_vector),            intent(in) :: U0(:)
      real(kind=wp),                 intent(in) :: S0(:,:)
      ! vector of dt values
      real(kind=wp),                 intent(in) :: tauv(:)
      ! vector of rank values
      integer,                       intent(in) :: rkv(:)
      real(kind=wp),                 intent(in) :: Tend
      integer,                       intent(in) :: nrep
      ! Optional
      logical, optional,             intent(in) :: ifsave
      logical                                   :: if_save_npy
      logical, optional,             intent(in) :: ifverb
      logical                                   :: verb
      logical, optional,             intent(in) :: iflogs
      logical                                   :: if_save_logs
      
      ! Internal variables
      type(LR_state),     allocatable           :: X     ! Controllability
      type(LR_state),     allocatable           :: Y     ! Observability
      real(kind=wp)                             :: U_out(2*nx,rkmax)
      real(kind=wp)                             :: X_out(2*nx,2*nx)
      real(kind=wp),      allocatable           :: vecs(:,:)
      real(kind=wp),      allocatable           :: vals(:)
      real(kind=wp)                             :: sfro
      real(kind=wp)                             :: tau, Ttot, etime
      integer                                   :: i, j, k, rk, irep, nsteps
      integer                                   :: info, torder, iostatus
      real(kind=wp)                             :: lagsvd(rkmax)
      character*128      :: oname
      character*128      :: onameU
      character*128      :: onameS
      integer            :: clock_rate, clock_start, clock_stop

      if_save_npy  = optval(ifsave, .false.)
      verb         = optval(ifverb, .false.)
      if_save_logs = optval(iflogs, .false.)

      call system_clock(count_rate=clock_rate)

      write(*,*) '----------------------'
      write(*,*) '   CONTROLLABILITY'
      write(*,*) '----------------------'
      do torder = 1, 2
         do i = 1, size(rkv)
            rk = rkv(i)
            if (allocated(vals)) deallocate(vals)
            if (allocated(vecs)) deallocate(vecs)
            allocate(vals(1:rk))
            allocate(vecs(1:rk,1:rk))
            do j = 1, size(tauv)
               tau = tauv(j)
               ! Initialize low-rank representation with rank rk
               if (verb) write(*,*) 'Initialize LR state, rk =', rk
               X = LR_state()
               call X%initialize_LR_state(U0, S0, rk)
               ! Reset time
               Ttot = 0.0_wp
               lagsvd = 0.0_wp
               if (verb) write(*,*) 'Run DRLA'
               if (if_save_logs) then
                  write(oname,'("output_GL_X_norm__n",I4.4,"_TO",I1,"_rk",I2.2,"_t",E8.2,".txt")') nx, torder, rk, tau
                  open(unit=iunit1, file=trim(basepath)//oname)
                  call stamp_logfile_header(iunit1, 'Controllability Gramian', rk, tau, Tend, torder)
                  write(iunit1,'(A16,A4,A10,A16,A20)') 'DLRA:','  rk',' Tend','|| X_DLRA ||_2', 'Elapsed time'
                  write(oname,'("output_GL_X_sigma_n",I4.4,"_TO",I1,"_rk",I2.2,"_t",E8.2,".txt")') nx, torder, rk, tau
                  open(unit=iunit2, file=trim(basepath)//oname)
                  call stamp_logfile_header(iunit2, 'Controllability Gramian', rk, tau, Tend, torder)
                  write(iunit2,*) 'DLRA: T    sigma_i    d(sigma-i)/sigma-1    d(sigma_i)/sigma_i    ||Sum(sigma_i)||_2'
               end if
               write(*,'(A16,A4,A4,A10,A6,A8,A16,A20)') 'DLRA:','  rk',' TO','dt','steps','Tend','|| X_DLRA ||_2', 'Elapsed time'
               nsteps = nint(Tend/tau)
               do irep = 1, nrep
                  ! run integrator
                  etime = 0.0_wp
                  call system_clock(count=clock_start)     ! Start Timer
                  call numerical_low_rank_splitting_lyapunov_integrator(X, LTI%prop, LTI%B, Tend, tau, torder, info, &
                                                                     & exptA=exptA, iftrans=.false., ifverb=.false.)
                  call system_clock(count=clock_stop)      ! Stop Timer
                  etime = etime + real(clock_stop-clock_start)/real(clock_rate)
                  ! Compute LR basis spectrum
                  call dsval(X%S, vals)
                  if (if_save_logs) then
                     write(iunit2,'("sigma ",F8.4)',ADVANCE='NO') Ttot
                     do k = 1, rk
                        write(iunit2,'(E14.6)', ADVANCE='NO') vals(k)
                     end do
                     write (iunit2,'(A)', ADVANCE='NO') ' | '
                     do k = 1, rk
                        write(iunit2,'(E14.6)', ADVANCE='NO') abs(vals(k) - lagsvd(k))/lagsvd(1)
                     end do
                     write (iunit2,'(A)', ADVANCE='NO') ' | '
                     do k = 1, rk
                        write(iunit2,'(E14.6)', ADVANCE='NO') abs(vals(k) - lagsvd(k))/lagsvd(k)
                     end do
                     write (iunit2,'(A)', ADVANCE='NO') ' | '
                     lagsvd(1:rk) = lagsvd(1:rk) - vals
                     sfro = 0.0_wp
                     do k = 1, rk
                        sfro = sfro + lagsvd(k)**2
                     end do
                     sfro = sqrt(sfro)
                     write(iunit2,'(E14.6)'), sfro           
                  end if
                  lagsvd(1:rk) = vals

                  ! Reconstruct solution
                  call get_state(U_out(:,1:rk), X%U)
                  X_out = matmul(U_out(:,1:rk), matmul(X%S, transpose(U_out(:,1:rk))))

                  Ttot = Ttot + Tend
                  write(*,'(I4," ",A11,I4," TO",I1,F10.6,I6,F8.4,E16.8,F18.4," s")') irep, 'Xctl OUTPUT', &
                                    & rk, torder, tau, nsteps, Ttot, norm2(X_out), etime
                  if (if_save_logs) then
                     write(iunit1,'(I4," ",A11,I6,F8.4,E16.8,F18.4," s")') irep, 'Xctl OUTPUT', nsteps, Ttot, norm2(X_out), etime
                  end if
               end do
               if (if_save_logs) then
                  close(iunit1)
                  close(iunit2)
               end if
               if (if_save_npy) then
                  write(oname,'("GL_Xdata_TO",I1,"_rk",I2.2,"_t",I1,".npy")') torder, rk, j
                  call save_npy(trim(basepath)//oname, X_out, iostatus)
                  if (iostatus /= 0) then; write(*,*) "Error saving file ", trim(oname); STOP 2; end if
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
         do i = 1, size(rkv)
            rk = rkv(i)
            do j = 1, size(tauv)
               tau = tauv(j)
               ! Initialize low-rank representation with rank rk
               if (verb) write(*,*) 'Initialize LR state, rk =', rk
               Y = LR_state()
               call Y%initialize_LR_state(U0, S0, rk)
               ! Reset time
               Ttot = 0.0_wp
               if (verb) write(*,*) 'Run DRLA'
               if (if_save_logs) then
                  write(oname,'("output_GL_Y_norm__n",I4.4,"_TO",I1,"_rk",I2.2,"_t",E8.2,".txt")') nx, torder, rk, tau
                  open(unit=iunit1, file=trim(basepath)//oname)
                  call stamp_logfile_header(iunit1, 'Observability Gramian', rk, tau, Tend, torder)
                  write(iunit1,'(A16,A4,A10,A16,A20)') 'DLRA:','  rk',' Tend','|| X_DLRA ||_2', 'Elapsed time'
                  write(oname,'("output_GL_Y_sigma_n",I4.4,"_TO",I1,"_rk",I2.2,"_t",E8.2,".txt")') nx, torder, rk, tau
                  open(unit=iunit2, file=trim(basepath)//oname)
                  call stamp_logfile_header(iunit2, 'Observability Gramian', rk, tau, Tend, torder)
                  write(iunit2,*) 'DLRA: T    sigma_i    d(sigma-i)/sigma-1    d(sigma_i)/sigma_i    ||Sum(sigma_i)||_2'
               end if
               write(*,'(A16,A4,A4,A10,A6,A8,A16,A20)') 'DLRA:','  rk',' TO','dt','steps','Tend','|| X_DLRA ||_2', 'Elapsed time'
               nsteps = nint(Tend/tau)
               do irep = 1, nrep
                  ! run integrator
                  etime = 0.0_wp
                  call system_clock(count=clock_start)     ! Start Timer
                  call numerical_low_rank_splitting_lyapunov_integrator(X, LTI%prop, LTI%CT, Tend, tau, torder, info, &
                                                                     & exptA=exptA, iftrans=.false., ifverb=.false.)
                  call system_clock(count=clock_stop)      ! Stop Timer
                  etime = etime + real(clock_stop-clock_start)/real(clock_rate)
                  ! Compute LR basis spectrum
                  call dsval(X%S, vals)
                  if (if_save_logs) then
                     write(iunit2,'("sigma ",F8.4)',ADVANCE='NO') Ttot
                     do k = 1, rk
                        write(iunit2,'(E14.6)', ADVANCE='NO') vals(k)
                     end do
                     write (iunit2,'(A)', ADVANCE='NO') ' | '
                     do k = 1, rk
                        write(iunit2,'(E14.6)', ADVANCE='NO') abs(vals(k) - lagsvd(k))/lagsvd(1)
                     end do
                     write (iunit2,'(A)', ADVANCE='NO') ' | '
                     do k = 1, rk
                        write(iunit2,'(E14.6)', ADVANCE='NO') abs(vals(k) - lagsvd(k))/lagsvd(k)
                     end do
                     write (iunit2,'(A)', ADVANCE='NO') ' | '
                     lagsvd(1:rk) = lagsvd(1:rk) - vals
                     sfro = 0.0_wp
                     do k = 1, rk
                        sfro = sfro + lagsvd(k)**2
                     end do
                     sfro = sqrt(sfro)
                     write(iunit2,'(E14.6)'), sfro          
                  end if
                  lagsvd(1:rk) = vals

                  ! Reconstruct solution
                  call get_state(U_out(:,1:rk), Y%U)
                  X_out = matmul(U_out(:,1:rk), matmul(Y%S, transpose(U_out(:,1:rk))))

                  Ttot = Ttot + Tend
                  write(*,'(I4," ",A11,I4," TO",I1,F10.6,I6,F8.4,E16.8,F18.4," s")') irep, 'Yobs OUTPUT', &
                                    & rk, torder, tau, nsteps, Ttot, norm2(X_out), etime
                  if (if_save_logs) then
                     write(iunit1,'(I4," ",A11,I6,F8.4,E16.8,F18.4," s")') irep, 'Yobs OUTPUT', nsteps, Ttot, norm2(X_out), etime
                  end if
               end do
               if (if_save_logs) then
                  close(iunit1)
                  close(iunit2)
               end if
               if (if_save_npy) then
                  write(oname,'("GL_Ydata_TO",I1,"_rk",I2.2,"_t",I1,".npy")') torder, rk, j
                  call save_npy(trim(basepath)//oname, X_out, iostatus)
                  if (iostatus /= 0) then; write(*,*) "Error saving file ", trim(oname); STOP 2; end if
               end if
               deallocate(Y%U)
               deallocate(Y%S)
            end do
         end do
      end do

      return
   end subroutine run_DLRA_test

   !
   ! BALANCED TRUNCATION
   !

   subroutine run_BT_test(LTI, U0, S0, rk, tau, torder, Tmax, nrep, ifsave, ifload, ifverb, iflogs)
      ! LTI system
      type(lti_system),              intent(inout) :: LTI
      ! Initial condition
      type(state_vector),            intent(in) :: U0(:)
      real(kind=wp),                 intent(in) :: S0(:,:)
      integer,                       intent(in) :: rk
      real(kind=wp),                 intent(inout) :: tau
      integer,                       intent(inout) :: torder
      real(kind=wp),                 intent(in) :: Tmax
      integer,                       intent(in) :: nrep
      ! Optional
      logical, optional,             intent(in) :: ifsave
      logical                                   :: if_save_npy
      logical, optional,             intent(in) :: ifload
      logical                                   :: if_load_npy
      logical, optional,             intent(in) :: ifverb
      logical                                   :: verb
      logical, optional,             intent(in) :: iflogs
      logical                                   :: if_save_logs

      ! Internal variables
      type(LR_state),     allocatable           :: X     ! Controllability
      type(LR_state),     allocatable           :: Y     ! Observability
      type(state_vector), allocatable           :: Utmp(:)
      real(kind=wp)                             :: U0_in(2*nx, rkmax)
      real(kind=wp)                             :: U_out(2*nx,rkmax)
      real(kind=wp)                             :: X_out(2*nx,2*nx)
      real(kind=wp),      allocatable           :: vecs(:,:)
      real(kind=wp),      allocatable           :: vals(:)
      real(kind=wp)                             :: sfro
      real(kind=wp)                             :: Tend, Ttot, etime, etime_tot
      integer                                   :: i, j, k, irep, nsteps
      integer                                   :: info, iostatus
      real(kind=wp)                             :: lagsvd(rkmax)


      ! ROM
      real(kind=wp),      allocatable           :: Swrk(:,:)
      real(kind=wp),      allocatable           :: Ahat(:,:)
      real(kind=wp),      allocatable           :: Bhat(:,:)
      real(kind=wp),      allocatable           :: Chat(:,:)
      real(kind=wp),      allocatable           :: D(:,:)

      ! BT
      type(state_vector), allocatable           :: T(:)
      type(state_vector), allocatable           :: Tinv(:)
      real(kind=wp),      allocatable           :: S(:)
      real(kind=wp),      allocatable           :: U_load(:,:)
      real(kind=wp),      allocatable           :: S_load(:,:)

      ! SVD
      real(kind=wp)  :: U_svd(2*nx,2*nx)
      real(kind=wp)  :: S_svd(rkmax)
      real(kind=wp)  :: V_svd(rkmax,rkmax)

      character*128      :: oname
      character*128      :: onameU
      character*128      :: onameS
      logical                                   :: existU, existS
      integer            :: clock_rate, clock_start, clock_stop

      if_save_npy  = optval(ifsave, .false.)
      if_load_npy  = optval(ifload, .false.)
      verb         = optval(ifverb, .false.)
      if_save_logs = optval(iflogs, .false.)

      call system_clock(count_rate=clock_rate)

      Tend   = Tmax/nrep
      nsteps = nint(Tend/tau)

      write(*,*) ''
      write(*,*) '----------------------'
      write(*,*) '   CONTROLLABILITY'
      write(*,*) '----------------------'
      write(*,*) ''

      onameU = trim(basepath)//"GL_Xctl_U.npy"
      onameS = trim(basepath)//"GL_Xctl_S.npy"

      X = LR_state()
      if (if_load_npy) then
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
         etime_tot = 0.0_wp
         write(*,'(A16,A4,A4,A10,A6,A8,A16,A20)') 'DLRA:','  rk',' TO','dt','steps','Tend','|| X_DLRA ||_2', 'Elapsed time'
         do irep = 1, nrep
            ! run integrator
            etime = 0.0_wp
            call system_clock(count=clock_start)     ! Start Timer
            call numerical_low_rank_splitting_lyapunov_integrator(X, LTI%prop, LTI%B, Tend, tau, torder, info, &
                                                               & exptA=exptA, iftrans=.false., ifverb=.false.)
            call system_clock(count=clock_stop)      ! Stop Timer
            etime = etime + real(clock_stop-clock_start)/real(clock_rate)

            ! Reconstruct solution
            call get_state(U_out(:,1:rk), X%U)
            X_out = matmul(U_out(:,1:rk), matmul(X%S, transpose(U_out(:,1:rk))))

            Ttot = Ttot + Tend
            write(*,'(I4," ",A11,I4," TO",I1,F10.6,I6,F8.4,E16.8,F18.4," s")') irep, 'Xctl OUTPUT', &
                              & rk, torder, tau, nsteps, Ttot, norm2(X_out), etime
            etime_tot = etime_tot + etime
         end do
         if (verb) write(*,*) 'Total integration time (DLRA):', etime_tot, 's'
         if (if_save_npy) then
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
      if (if_load_npy) then
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
         etime_tot = 0.0_wp
         write(*,'(A16,A4,A4,A10,A6,A8,A16,A20)') 'DLRA:','  rk',' TO','dt','steps','Tend','|| X_DLRA ||_2', 'Elapsed time'
         do irep = 1, nrep
            ! run integrator
            etime = 0.0_wp
            call system_clock(count=clock_start)     ! Start Timer
            call numerical_low_rank_splitting_lyapunov_integrator(Y, LTI%prop, LTI%CT, Tend, tau, torder, info, &
                                                               & exptA=exptA, iftrans=.true., ifverb=.false.)
            call system_clock(count=clock_stop)      ! Stop Timer
            etime = etime + real(clock_stop-clock_start)/real(clock_rate)

            ! Reconstruct solution
            call get_state(U_out(:,1:rk), Y%U)
            X_out = matmul(U_out(:,1:rk), matmul(Y%S, transpose(U_out(:,1:rk))))

            Ttot = Ttot + Tend
            write(*,'(I4," ",A11,I4," TO",I1,F10.6,I6,F8.4,E16.8,F18.4," s")') irep, 'Yobs OUTPUT', &
                              & rk, torder, tau, nsteps, Ttot, norm2(X_out), etime
            etime_tot = etime_tot + etime
         end do
         if (verb) write(*,*) 'Total integration time (DLRA):', etime_tot, 's'
         if (if_save_npy) then
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

      if (if_save_npy) then
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

   end subroutine run_BT_test

   !
   ! Logfiles
   !

   subroutine stamp_logfile_header(iunit, problem, rk, tau, Tend, torder)
      integer,       intent(in) :: iunit
      character(*),  intent(in) :: problem
      integer,       intent(in) :: rk
      real(kind=wp), intent(in) :: tau
      real(kind=wp), intent(in) :: Tend
      integer,       intent(in) :: torder

      write(iunit,*) '-----------------------'
      write(iunit,*) '    GINZBURG LANDAU'
      write(iunit,*) '-----------------------'
      write(iunit,*) 'nu    = ', nu
      write(iunit,*) 'gamma = ', gamma
      write(iunit,*) 'mu_0  = ', mu_0
      write(iunit,*) 'c_mu  = ', c_mu
      write(iunit,*) 'mu_2  = ', mu_2
      write(iunit,*) '-----------------------'
      write(iunit,*) problem
      write(iunit,*) '-----------------------'
      write(iunit,*) 'nx    = ', nx
      write(iunit,*) 'rk_b  = ', rk_b
      write(iunit,*) 'x_b   = ', x_b
      write(iunit,*) 's_b   = ', s_b
      write(iunit,*) 'rk_c  = ', rk_c
      write(iunit,*) 'x_c   = ', x_c
      write(iunit,*) 's_c   = ', s_c
      write(iunit,*) '-----------------------'
      write(iunit,*) 'Time Integration: DLRA'
      write(iunit,*) '-----------------------'
      write(iunit,*) 'Tend   =', Tend
      write(iunit,*) 'torder =', torder
      write(iunit,*) 'tau    =', tau
      write(iunit,*) 'rk     =', rk
      write(iunit,*) '---------------------'
      write(iunit,*) '---------------------'

      return
   end subroutine stamp_logfile_header



end module Ginzburg_Landau_Utils