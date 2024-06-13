module Ginzburg_Landau_Tests
   ! Standard Library.
   use stdlib_optval, only : optval
   use stdlib_linalg, only : diag
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
   use Ginzburg_Landau_RK_Lyapunov
   use Ginzburg_Landau_Utils
   !use fortime
   implicit none

   ! IO
   integer,       parameter :: iunit1 = 1
   integer,       parameter :: iunit2 = 2
   character*128, parameter :: basepath = 'local/'
   integer,       parameter :: rkmax = 40
   integer,       parameter :: rk_X0 = 40

   private :: this_module
   public  :: iunit1, iunit2, basepath, rkmax, rk_X0
   public  :: run_DLRA_lyapunov_test
   public  :: run_BT_test
   public  :: run_kexpm_test
   public  :: run_DLRA_riccati_test
   public  :: run_lyap_convergence_test

   character*128, parameter :: this_module = 'Ginzburg_Landau_Tests'

contains

   subroutine run_DLRA_lyapunov_test(LTI, U0, S0, rkv, tauv, TOv, Tend, nrep, ifsave, ifverb, iflogs)
      ! LTI system
      type(lti_system),              intent(inout) :: LTI
      ! Initial condition
      type(state_vector),            intent(in) :: U0(:)
      real(wp),                      intent(in) :: S0(:,:)
      ! vector of dt values
      real(wp),                      intent(in) :: tauv(:)
      ! vector of rank values
      integer,                       intent(in) :: rkv(:)
      ! vector of torders
      integer,                       intent(in) :: TOv(:)
      real(wp),                      intent(in) :: Tend
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
      real(wp)                                  :: U_out(2*nx,rkmax)
      real(wp)                                  :: X_out(2*nx,2*nx)
      real(wp),           allocatable           :: vals(:)
      real(wp)                                  :: sfro
      real(wp)                                  :: tau, Ttot, etime, etime_tot
      integer                                   :: i, j, k, ito, rk, irep, nsteps
      integer                                   :: info, torder, iostatus
      real(wp)                                  :: lagsvd(rkmax)
      real(wp)                                  :: res(N**2)
      character*128      :: oname
      character*128      :: onameU
      character*128      :: onameS
      integer            :: clock_rate, clock_start, clock_stop

      if_save_npy  = optval(ifsave, .false.)
      verb         = optval(ifverb, .false.)
      if_save_logs = optval(iflogs, .false.)

      call system_clock(count_rate=clock_rate)

      write(*,*) ''
      write(*,*) '----------------------'
      write(*,*) '   CONTROLLABILITY'
      write(*,*) '----------------------'
      write(*,*) ''

      X = LR_state()
      do ito = 1, size(TOv)
         torder = TOv(ito)
         do i = 1, size(rkv)
            rk = rkv(i)
            if (allocated(vals)) deallocate(vals)
            allocate(vals(1:rk))
            do j = 1, size(tauv)
               tau = tauv(j)
               ! Initialize low-rank representation with rank rk
               if (verb) write(*,*) 'Initialize LR state, rk =', rk
               call X%initialize_LR_state(U0, S0, rk)
               ! Reset time
               Ttot = 0.0_wp
               lagsvd = 0.0_wp
               if (verb) write(*,*) 'Run DRLA'
               if (if_save_logs) then
                  write(oname,'("output_GL_X_norm__n",I4.4,"_TO",I1,"_rk",I2.2,"_t",E8.2,".txt")') nx, torder, rk, tau
                  open(unit=iunit1, file=trim(basepath)//oname)
                  call stamp_logfile_header(iunit1, 'Controllability Gramian', rk, tau, Tend, torder)
                  write(iunit1,'(A16,A4,A10,A18,A18,A20)') 'DLRA:','  rk',' Tend','|| X_DLRA ||_2/N','|| res ||_2/N', 'Elapsed time'
                  write(oname,'("output_GL_X_sigma_n",I4.4,"_TO",I1,"_rk",I2.2,"_t",E8.2,".txt")') nx, torder, rk, tau
                  open(unit=iunit2, file=trim(basepath)//oname)
                  call stamp_logfile_header(iunit2, 'Controllability Gramian', rk, tau, Tend, torder)
                  write(iunit2,*) 'DLRA: T    sigma_i    d(sigma-i)/sigma-1    d(sigma_i)/sigma_i    ||Sum(sigma_i)||_2'
               end if
               write(*,'(A16,A4,A4,A10,A6,A8,A18,A18,A20)') 'DLRA:','  rk',' TO','dt','steps','Tend', &
                                    & '|| X_DLRA ||_2/N','|| res ||_2/N', 'Elapsed time'
               nsteps = nint(Tend/tau)
               etime_tot = 0.0_wp
               do irep = 1, nrep
                  ! run integrator
                  etime = 0.0_wp
                  call system_clock(count=clock_start)     ! Start Timer
                  call numerical_low_rank_splitting_lyapunov_integrator(X, LTI%prop, LTI%B, Tend, tau, torder, info, &
                                                                        & exptA=exptA, iftrans=.false., ifverb=.false.)
                  call system_clock(count=clock_stop)      ! Stop Timer
                  etime = etime + real(clock_stop-clock_start)/real(clock_rate)
                  ! Compute LR basis spectrum
                  call sval(X%S(1:rk,1:rk), vals)
                  if (if_save_logs) then
                     write(iunit2,'("sigma ",F8.4)',ADVANCE='NO') Ttot
                     do k = 1, rk; write(iunit2,'(E14.6)', ADVANCE='NO') vals(k); end do
                     write (iunit2,'(A)', ADVANCE='NO') ' | '
                     do k = 1, rk; write(iunit2,'(E14.6)', ADVANCE='NO') abs(vals(k) - lagsvd(k))/lagsvd(1); end do
                     write (iunit2,'(A)', ADVANCE='NO') ' | '
                     do k = 1, rk; write(iunit2,'(E14.6)', ADVANCE='NO') abs(vals(k) - lagsvd(k))/lagsvd(k); end do
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
                  call reconstruct_solution(X_out, X)
                  Ttot = Ttot + Tend
                  call CALE(res, reshape(X_out, shape(res)), BBTW_flat, .false.)
                  write(*,'(I4," ",A11,I4," TO",I1,F10.6,I6,F8.4,E18.8,E18.8,F18.4," s")') irep, 'Xctl OUTPUT', &
                                    & rk, torder, tau, nsteps, Ttot, norm2(X_out)/N, norm2(res)/N, etime
                  if (if_save_logs) then
                     write(iunit1,'(I4," ",A11,I6,F8.4,E18.8,E18.8,E18.8,F18.4," s")') irep, 'Xctl OUTPUT', &
                                    & nsteps, Ttot, norm2(X_out)/N, norm2(res)/N, etime
                  end if
                  etime_tot = etime_tot + etime
               end do
               if (verb) write(*,*) 'Total integration time (DLRA):', etime_tot, 's'
               if (if_save_logs) then
                  write(iunit1,*) 'Total integration time (DLRA):', etime_tot, 's'; close(iunit1)
                  write(Iunit2,*) 'Total integration time (DLRA):', etime_tot, 's'; close(iunit2)
               end if
               if (if_save_npy) then
                  write(onameU,'("data_GLXY_XU_n",I4.4,"_TO",I1,"_rk",I2.2,"_t",E8.2,".npy")') nx, torder, rk, tau
                  write(onameS,'("data_GLXY_XS_n",I4.4,"_TO",I1,"_rk",I2.2,"_t",E8.2,".npy")') nx, torder, rk, tau
                  call save_npy(trim(basepath)//onameU, U_out(:,1:rk), iostatus)
                  if (iostatus /= 0) then; write(*,*) "Error saving file", trim(onameU); STOP 2; end if
                  call save_npy(trim(basepath)//onameS, X%S(1:rk,1:rk), iostatus)
                  if (iostatus /= 0) then; write(*,*) "Error saving file", trim(onameS); STOP 2; end if
               end if
               deallocate(X%U)
               deallocate(X%S)
            end do
         end do
      end do

      STOP 1
      write(*,*) ''
      write(*,*) '--------------------'
      write(*,*) '   OBSERVABILITY'
      write(*,*) '--------------------'
      write(*,*) ''

      Y = LR_state()
      do ito = 1, size(TOv)
         torder = TOv(ito)
         do i = 1, size(rkv)
            rk = rkv(i)
            if (allocated(vals)) deallocate(vals)
            allocate(vals(1:rk))
            do j = 1, size(tauv)
               tau = tauv(j)
               ! Initialize low-rank representation with rank rk
               if (verb) write(*,*) 'Initialize LR state, rk =', rk
               call Y%initialize_LR_state(U0, S0, rk)
               ! Reset time
               Ttot = 0.0_wp
               if (verb) write(*,*) 'Run DRLA'
               if (if_save_logs) then
                  write(oname,'("output_GL_Y_norm__n",I4.4,"_TO",I1,"_rk",I2.2,"_t",E8.2,".txt")') nx, torder, rk, tau
                  open(unit=iunit1, file=trim(basepath)//oname)
                  call stamp_logfile_header(iunit1, 'Observability Gramian', rk, tau, Tend, torder)
                  write(iunit1,'(A16,A4,A10,A18,A18,A20)') 'DLRA:','  rk',' Tend','|| X_DLRA ||_2/N','|| res ||_2/N', 'Elapsed time'
                  write(oname,'("output_GL_Y_sigma_n",I4.4,"_TO",I1,"_rk",I2.2,"_t",E8.2,".txt")') nx, torder, rk, tau
                  open(unit=iunit2, file=trim(basepath)//oname)
                  call stamp_logfile_header(iunit2, 'Observability Gramian', rk, tau, Tend, torder)
                  write(iunit2,*) 'DLRA: T    sigma_i    d(sigma-i)/sigma-1    d(sigma_i)/sigma_i    ||Sum(sigma_i)||_2'
               end if
               write(*,'(A16,A4,A4,A10,A6,A8,A18,A18,A20)') 'DLRA:','  rk',' TO','dt','steps','Tend', &
                                    & '|| X_DLRA ||_2/N','|| res ||_2/N', 'Elapsed time'
               nsteps = nint(Tend/tau)
               etime_tot = 0.0_wp
               do irep = 1, nrep
                  ! run integrator
                  etime = 0.0_wp
                  call system_clock(count=clock_start)     ! Start Timer
                  call numerical_low_rank_splitting_lyapunov_integrator(Y, LTI%prop, LTI%CT, Tend, tau, torder, info, & 
                                                                        & exptA=exptA, iftrans=.true., ifverb=.false.)
                  call system_clock(count=clock_stop)      ! Stop Timer
                  etime = etime + real(clock_stop-clock_start)/real(clock_rate)
                  ! Compute LR basis spectrum
                  call sval(Y%S(1:rk,1:rk), vals)
                  if (if_save_logs) then
                     write(iunit2,'("sigma ",F8.4)',ADVANCE='NO') Ttot
                     do k = 1, rk; write(iunit2,'(E14.6)', ADVANCE='NO') vals(k); end do
                     write (iunit2,'(A)', ADVANCE='NO') ' | '
                     do k = 1, rk; write(iunit2,'(E14.6)', ADVANCE='NO') abs(vals(k) - lagsvd(k))/lagsvd(1); end do
                     write (iunit2,'(A)', ADVANCE='NO') ' | '
                     do k = 1, rk; write(iunit2,'(E14.6)', ADVANCE='NO') abs(vals(k) - lagsvd(k))/lagsvd(k); end do
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
                  call reconstruct_solution(X_out, Y)
                  Ttot = Ttot + Tend
                  call CALE(res, reshape(X_out, shape(res)), CTCW_flat, .true.)
                  write(*,'(I4," ",A11,I4," TO",I1,F10.6,I6,F8.4,E18.8,E18.8,F18.4," s")') irep, 'Yobs OUTPUT', &
                                    & rk, torder, tau, nsteps, Ttot, norm2(X_out)/N, norm2(res)/N, etime
                  if (if_save_logs) then
                     write(iunit1,'(I4," ",A11,I6,F8.4,E18.8,E18.8,F18.4," s")') irep, 'Yobs OUTPUT', &
                                    & nsteps, Ttot, norm2(X_out)/N, norm2(res)/N, etime
                  end if
                  etime_tot = etime_tot + etime
               end do
               if (verb) write(*,*) 'Total integration time (DLRA):', etime_tot, 's'
               if (if_save_logs) then
                  write(iunit1,*) 'Total integration time (DLRA):', etime_tot, 's'; close(iunit1)
                  write(Iunit2,*) 'Total integration time (DLRA):', etime_tot, 's'; close(iunit2)
               end if
               if (if_save_npy) then
                  write(onameU,'("data_GLXY_YU_n",I4.4,"_TO",I1,"_rk",I2.2,"_t",E8.2,".npy")') nx, torder, rk, tau
                  write(onameS,'("data_GLXY_YS_n",I4.4,"_TO",I1,"_rk",I2.2,"_t",E8.2,".npy")') nx, torder, rk, tau
                  call save_npy(trim(basepath)//onameU, U_out(:,1:rk), iostatus)
                  if (iostatus /= 0) then; write(*,*) "Error saving file", trim(onameU); STOP 2; end if
                  call save_npy(trim(basepath)//onameS, Y%S(1:rk,1:rk), iostatus)
                  if (iostatus /= 0) then; write(*,*) "Error saving file", trim(onameS); STOP 2; end if
               end if
               deallocate(Y%U)
               deallocate(Y%S)
            end do
         end do
      end do

      return
   end subroutine run_DLRA_lyapunov_test

   !
   ! BALANCED TRUNCATION
   !

   subroutine run_BT_test(LTI, U0, S0, rk, tau, torder, Tmax, nrep, ifsave, ifload, ifverb, iflogs)
      ! LTI system
      type(lti_system),              intent(inout) :: LTI
      ! Initial condition
      type(state_vector),            intent(in) :: U0(:)
      real(wp),                      intent(in) :: S0(:,:)
      integer,                       intent(in) :: rk
      real(wp),                      intent(inout) :: tau
      integer,                       intent(inout) :: torder
      real(wp),                      intent(in) :: Tmax
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
      real(wp)                                  :: U0_in(2*nx, rkmax)
      real(wp)                                  :: U_out(2*nx,rkmax)
      real(wp)                                  :: X_out(2*nx,2*nx)
      real(wp),           allocatable           :: vecs(:,:)
      real(wp),           allocatable           :: vals(:)
      real(wp)                                  :: sfro
      real(wp)                                  :: Tend, Ttot, etime, etime_tot
      integer                                   :: i, j, k, irep, nsteps
      integer                                   :: info, iostatus
      real(wp)                                  :: lagsvd(rkmax)
      real(wp)                                  :: res(N**2)

      ! ROM
      real(wp),           allocatable           :: Swrk(:,:)
      real(wp),           allocatable           :: Ahat(:,:)
      real(wp),           allocatable           :: Bhat(:,:)
      real(wp),           allocatable           :: Chat(:,:)
      real(wp),           allocatable           :: D(:,:)

      ! BT
      type(state_vector), allocatable           :: T(:)
      type(state_vector), allocatable           :: Tinv(:)
      real(wp),           allocatable           :: S(:)
      real(wp),           allocatable           :: U_load(:,:)
      real(wp),           allocatable           :: S_load(:,:)

      ! SVD
      real(wp)       :: U_svd(2*nx,2*nx)
      real(wp)       :: S_svd(rkmax)
      real(wp)       :: V_svd(rkmax,rkmax)

      character*128      :: oname
      character*128      :: onameU
      character*128      :: onameS
      logical            :: existU, existS
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
         write(*,'(A16,A4,A4,A10,A6,A8,A18,A18,A20)') 'DLRA:','  rk',' TO','dt','steps','Tend', &
                           & '|| X_DLRA ||_2/N','|| res ||_2/N', 'Elapsed time'
         do irep = 1, nrep
            ! run integrator
            etime = 0.0_wp
            call system_clock(count=clock_start)     ! Start Timer
            call numerical_low_rank_splitting_lyapunov_integrator(X, LTI%prop, LTI%B, Tend, tau, torder, info, &
                                                                  & exptA=exptA, iftrans=.false., ifverb=.false.)
            call system_clock(count=clock_stop)      ! Stop Timer
            etime = etime + real(clock_stop-clock_start)/real(clock_rate)

            ! Reconstruct solution
            call reconstruct_solution(X_out, X)
            Ttot = Ttot + Tend
            call CALE(res, reshape(X_out, shape(res)), BBTW_flat, .false.)
            write(*,'(I4," ",A11,I4," TO",I1,F10.6,I6,F8.4,E18.8,E18.8,F18.4," s")') irep, 'Xctl OUTPUT', &
                              & rk, torder, tau, nsteps, Ttot, norm2(X_out)/N, norm2(res)/N, etime
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
         write(*,'(A16,A4,A4,A10,A6,A8,A18,A18,A20)') 'DLRA:','  rk',' TO','dt','steps','Tend', &
                                    & '|| X_DLRA ||_2/N','|| res ||_2/N', 'Elapsed time'
         do irep = 1, nrep
            ! run integrator
            etime = 0.0_wp
            call system_clock(count=clock_start)     ! Start Timer
            call numerical_low_rank_splitting_lyapunov_integrator(Y, LTI%prop, LTI%CT, Tend, tau, torder, info, &
                                                                  & exptA=exptA, iftrans=.true., ifverb=.false.)
            call system_clock(count=clock_stop)      ! Stop Timer
            etime = etime + real(clock_stop-clock_start)/real(clock_rate)

            ! Reconstruct solution
            call reconstruct_solution(X_out, Y)
            Ttot = Ttot + Tend
            call CALE(res, reshape(X_out, shape(res)), CTCW_flat, .true.)
            write(*,'(I4," ",A11,I4," TO",I1,F10.6,I6,F8.4,E18.8,E18.8,F18.4," s")') irep, 'Yobs OUTPUT', &
                              & rk, torder, tau, nsteps, Ttot, norm2(X_out)/N, norm2(res)/N, etime
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
      call sqrtm(X%S(1:rk,1:rk), Swrk(1:rk,1:rk), info)
      block
         class(abstract_vector_rdp), allocatable :: Xwrk(:)
         call linear_combination(Xwrk, X%U(1:rk), Swrk(1:rk,1:rk))
         call copy_basis(Utmp, Xwrk)
      end block
      !call linear_combination(Utmp, X%U, Swrk(1:rk,1:rk))
      call get_state(U0_in(:,1:rk), Utmp)
      ! compute SVD of updated X%U
      call svd(U0_in(:,1:rk), U_svd(:,1:2*nx), S_svd(1:rk), V_svd(1:rk,1:rk))
      call set_state(X%U(1:rk), matmul(U_svd(:,1:rk), diag(S_svd(1:rk))))

      ! compute sqrt of coefficient matrix Y%S and right-multiply it to Y%U
      Swrk = 0.0_wp
      call sqrtm(Y%S(1:rk,1:rk), Swrk(1:rk,1:rk), info)
      block
         class(abstract_vector_rdp), allocatable :: Xwrk(:)
         call linear_combination(Xwrk, Y%U(1:rk), Swrk(1:rk,1:rk))
         call copy_basis(Utmp, Xwrk)
      end block
      !call linear_combination(Utmp, Y%U, Swrk(1:rk,1:rk))
      call get_state(U0_in(:,1:rk), Utmp)
      ! compute SVD of updated Y%U
      call svd(U0_in(:,1:rk), U_svd(:,1:2*nx), S_svd(1:rk), V_svd(1:rk,1:rk)) 
      call set_state(Y%U(1:rk), matmul(U_svd(:,1:rk), diag(S_svd(1:rk))))   

      ! compute balancing transformation based on SVD of Gramians
      allocate(T(1:rk), source=U0(1)); allocate(Tinv(1:rk), source=U0(1)); allocate(S(1:rk))
      call Balancing_Transformation(T, S, Tinv, X%U(1:rk), Y%U(1:rk))
      
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

   subroutine run_DLRA_riccati_test(LTI, U0, S0, Qc, Rinv, rkv, tauv, Tend, nrep, ifsave, ifverb, iflogs)
      type(lti_system),              intent(inout) :: LTI
      !! Considered LTI system
      type(state_vector),            intent(in)    :: U0(:)
      real(wp),                      intent(in)    :: S0(:,:)
      !! Initial condition
      real(wp),                      intent(in)    :: Qc(:,:)
      !! Measurement weights.
      real(wp),                      intent(in)    :: Rinv(:,:)
      !! Inverse of the actuation weights.
      real(wp),                      intent(in)    :: tauv(:)
      !! vector of dt values
      integer,                       intent(in)    :: rkv(:)
      !! vector of rank values
      real(wp),                      intent(in)    :: Tend
      integer,                       intent(in)    :: nrep
      ! Optional
      logical, optional,             intent(in)    :: ifsave
      logical                                      :: if_save_npy
      logical, optional,             intent(in)    :: ifverb
      logical                                      :: verb
      logical, optional,             intent(in)    :: iflogs
      logical                                      :: if_save_logs
      
      ! Internal variables
      type(LR_state),     allocatable              :: X     ! Controllability
      type(LR_state),     allocatable              :: Y     ! Observability
      real(wp)                                     :: U_out(2*nx,rkmax)
      real(wp)                                     :: X_out(2*nx,2*nx)
      real(wp),           allocatable              :: vecs(:,:)
      real(wp),           allocatable              :: vals(:)
      real(wp)                                     :: sfro
      real(wp)                                     :: tau, Ttot, etime, etime_tot
      integer                                      :: i, j, k, rk, irep, nsteps
      integer                                      :: info, torder, iostatus
      real(wp)                                     :: res(N**2)
      character*128      :: oname
      character*128      :: onameU
      character*128      :: onameS
      integer            :: clock_rate, clock_start, clock_stop

      if_save_npy  = optval(ifsave, .false.)
      verb         = optval(ifverb, .false.)
      if_save_logs = optval(iflogs, .false.)

      call system_clock(count_rate=clock_rate)

      write(*,*) '----------------------'
      write(*,*) '   RICCATI EQUATION'
      write(*,*) '----------------------'

      X = LR_state()
      do torder = 1, 1
         do i = 1, size(rkv)
            rk = rkv(i)
            do j = 1, size(tauv)
               tau = tauv(j)
               ! Initialize low-rank representation with rank rk
               if (verb) write(*,*) 'Initialize LR state, rk =', rk
               call X%initialize_LR_state(U0, S0, rk)
               ! Reset time
               Ttot = 0.0_wp
               if (verb) write(*,*) 'Run DRLA'
               write(*,'(A16,A4,A4,A10,A6,A8,A18,A18,A20)') 'DLRA:','  rk',' TO','dt','steps','Tend', &
                              & '|| X_DLRA ||_2', '|| res ||_2','Elapsed time'
               nsteps = nint(Tend/tau)
               etime_tot = 0.0_wp
               do irep = 1, nrep
                  ! run integrator
                  etime = 0.0_wp
                  call system_clock(count=clock_start)     ! Start Timer
                  call numerical_low_rank_splitting_riccati_integrator(X, LTI%prop, LTI%B, LTI%CT, Qc, Rinv, &
                                                                     & Tend, tau, torder, info, &
                                                                     & exptA=exptA, iftrans=.false., ifverb=verb)
                  call system_clock(count=clock_stop)      ! Stop Timer
                  etime = etime + real(clock_stop-clock_start)/real(clock_rate)

                  ! Reconstruct solution
                  call reconstruct_solution(X_out, X)
                  Ttot = Ttot + Tend
                  call CARE(res, reshape(X_out, shape(res)), reshape(CTQcCW_mat, shape(res)), BRinvBTW_mat, .false.)
                  write(*,'(I4," ",A11,I4," TO",I1,F10.6,I6,F8.4,E18.8,E18.8,F18.4," s")') irep, 'Xricc OUTPUT', &
                                    & rk, torder, tau, nsteps, Ttot, norm2(X_out)/N, norm2(res)/N, etime
                  etime_tot = etime_tot + etime
               end do
               if (verb) write(*,*) 'Total integration time (DLRA):', etime_tot, 's'
               if (if_save_npy) then
                  write(onameU,'("data_GL_Riccati_XU_n",I4.4,"_TO",I1,"_rk",I2.2,"_t",E8.2,".npy")') nx, torder, rk, tau
                  write(onameS,'("data_GL_Riccati_XS_n",I4.4,"_TO",I1,"_rk",I2.2,"_t",E8.2,".npy")') nx, torder, rk, tau
                  call save_npy(trim(basepath)//onameU, U_out(:,1:rk), iostatus)
                  if (iostatus /= 0) then; write(*,*) "Error saving file", trim(onameU); STOP 2; end if
                  call save_npy(trim(basepath)//onameS, X%S(1:rk,1:rk), iostatus)
                  if (iostatus /= 0) then; write(*,*) "Error saving file", trim(onameS); STOP 2; end if
               end if
               deallocate(X%U)
               deallocate(X%S)
            end do
         end do
      end do

   end subroutine run_DLRA_riccati_test


   subroutine run_kexpm_test(A, prop, U0, tauv, torder, N)
      class(abstract_linop_rdp),     intent(inout) :: A
      !! Linear operator: A
      class(abstract_linop_rdp),     intent(inout) :: prop
      !! Linear operator: exponential propagator
      class(abstract_vector_rdp),    intent(in)    :: U0
      !! Abstract vector as a source
      real(wp),                      intent(in)    :: tauv(:)
      !! vector of dt values
      integer,                       intent(inout) :: torder
      !! torder
      integer,                       intent(in)    :: n
      !! Number of repeats

      ! internal variables
      class(abstract_vector_rdp),    allocatable   :: U, V_kryl, V_rk   ! scratch bases
      integer                                      :: i, j, info
      real(wp)                                     :: tau
      real(wp)                                     :: tv_kryl(N), tv_rk(N)
      real(wp)                                     :: etime_kryl, etime_rk, stddev_kryl, stddev_rk
      integer            :: clock_rate, clock_start, clock_stop
      !type(timer)                                  :: tmr

      call system_clock(count_rate=clock_rate)
      
      allocate(U, source=U0); call U%zero()
      allocate(V_kryl, source=U0); call V_kryl%zero()
      allocate(V_rk,   source=U0); call V_rk%zero()

      tv_kryl = 0.0_wp
      tv_rk = 0.0_wp

      write(*,*) 'Comparison over N = ', N, 'runs.'
      write(*,'(A8," | ",2(10X,A10,10X,5X))') 'tau','KRYLOV','R-K'
      write(*,'(A8," | ",3(A10,1X),3X,3(A10,1X))') ' ','TOT','AVG','STDDEV','TOT','AVG','STDDEV'
      write(*,*) '------------------------------------------------------------------------------'

      do i = 1, size(tauv)
         tau = tauv(i)
         ! Reset time
         etime_kryl = 0.0_wp
         etime_rk   = 0.0_wp
         
         !call tmr%timer_start()
         !call U%rand()
         !call k_exptA(V_kryl, A, U, tau, info, .false.)
         !call tmr%timer_stop(nloops=N,message="test" , print=.true., color='red')
         !call tmr%timer_start()
         !call U%rand()
         !call exptA(V_rk, prop, U, tau, info, .false.)
         !call tmr%timer_stop(nloops=N, message="test" , print=.true., color='blue')

         do j = 1, N
            ! generate random vecor
            call U%rand()
            ! run Krylov based exponential propagator
            call system_clock(count=clock_start)     ! Start Timer
            call k_exptA(V_kryl, A, U, tau, info, .false.)
            call system_clock(count=clock_stop)      ! Stop Timer
            tv_kryl(j) = real(clock_stop-clock_start)/real(clock_rate)
            etime_kryl = etime_kryl + tv_kryl(j)
            ! run RK integrator
            call system_clock(count=clock_start)     ! Start Timer
            call exptA(V_rk, prop, U, tau, info, .false.)
            call system_clock(count=clock_stop)      ! Stop Timer
            tv_rk(j) = real(clock_stop-clock_start)/real(clock_rate)
            etime_rk = etime_rk + tv_rk(j)
            ! Check solution
            call V_kryl%axpby(1.0_wp, V_rk, -1.0_wp)
            if (V_kryl%norm()/(2*nx) > 10*atol_dp     ) then
               write(*,*) "Iteration", j, ": Solutions do not match!"
               write(*,* ) " tol", 10*atol_dp     , "delta = ", V_kryl%norm()/(2*nx)
            end if
         end do
         tv_kryl = tv_kryl - etime_kryl/N
         tv_rk = tv_rk - etime_rk/N
         stddev_kryl = 0.0_wp
         stddev_rk   = 0.0_wp
         do j = 1, N
            stddev_kryl = stddev_kryl + tv_kryl(j)**2
            stddev_rk   = stddev_rk + tv_rk(j)**2
         end do
         write(*,'(F8.6," | ",3(F10.6,1X),3X,3(F10.6,1X))') tau, &
                                 & etime_kryl, etime_kryl/N, sqrt(stddev_kryl/(N-1)), &
                                 & etime_rk,   etime_rk/N,   sqrt(stddev_rk/(N-1))
         
      end do

      return
   end subroutine run_kexpm_test

   subroutine run_lyap_convergence_test(LTI, U0, S0, Tend, tauv, rkv, TOv, ifsave, ifverb)
      ! LTI system
      type(lti_system),              intent(inout) :: LTI
      ! Initial condition
      type(state_vector),            intent(inout)    :: U0(:)
      real(wp),                      intent(inout)    :: S0(:,:)
      real(wp),                      intent(in)    :: Tend
      ! vector of dt values
      real(wp),                      intent(in)    :: tauv(:)
      ! vector of rank values
      integer,                       intent(in)    :: rkv(:)
      ! vector of torders
      integer,                       intent(in)    :: TOv(:)
      ! Optional
      logical, optional,             intent(in)    :: ifsave
      logical                                      :: if_save_npy
      logical, optional,             intent(in)    :: ifverb
      logical                                      :: verb

      ! Internals
      type(LR_state),                allocatable   :: X_state
      type(rk_lyapunov),             allocatable   :: RK_propagator
      type(state_matrix)                           :: X_mat(2)
      real(wp),                      allocatable   :: X_RKlib(:,:,:)
      real(wp)                                     :: X_mat_ref(N,N)
      real(wp)                                     :: U0_mat(N, rkmax)
      integer                                      :: info
      integer                                      :: i, j, ito, rk, nsteps, torder
      integer                                      :: irep, nrep
      real(wp)                                     :: etime, tau
      ! OUTPUT
      real(wp)                                     :: U_out(N,rkmax)
      real(wp)                                     :: X_out(N,N)
      real(wp)                                     :: Bmat(N,2)
      real(wp)                                     :: tmp(N,2)
      real(wp),           allocatable              :: U_load(:,:)
      real(wp)                                     :: X_mat_flat(N**2)
      real(wp)                                     :: res_flat(N**2)
      character*128      :: oname
      character*128      :: onameU
      character*128      :: onameS
      integer            :: iostatus
      ! timer
      integer            :: clock_rate, clock_start, clock_stop
      logical            :: existfile, load_data

      if_save_npy  = optval(ifsave, .false.)
      verb         = optval(ifverb, .false.)

      load_data = .false.

      call system_clock(count_rate=clock_rate)

      ! initialize exponential propagator
      RK_propagator = RK_lyapunov(Tend)

      oname = trim(basepath)//"Xctl_RK_W/data_BS_X_W.npy"
      inquire(file=oname, exist=existfile)
      if (existfile) then
         call load_npy(trim(oname), U_load, iostatus)
         print *, 'Load Bartels-Stuart solution ', trim(oname)
         if (iostatus /= 0) then; write(*,*) "Error loading file", trim(oname); STOP 2; end if
      else
         write(*,*) 'Cannot find ', trim(oname); STOP 12
      end if
      call CALE(res_flat, reshape(U_load, shape(res_flat)), BBTW_flat, .false.)
      print *, '   || res ||_2/N = ', norm2(res_flat)/N

      nrep = 1
      allocate(X_RKlib(N, N, nrep))
      if (load_data) then
         ! Load initial condition
         write(oname,'("data_GL_lyapconv_X0_RK_n",I4.4,".npy")') nx
         inquire(file=trim(basepath)//oname, exist=existfile)
         if (existfile) then
            call load_npy(trim(basepath)//trim(oname), U_load, iostatus)
            print *, 'Load initial condition: ', trim(basepath)//trim(oname)
            if (iostatus /= 0) then; write(*,*) "Error loading file", trim(basepath)//trim(oname); STOP 2; end if
         else
            write(*,*) 'Cannot find ', trim(basepath)//trim(oname); STOP 12
         end if
         X_out = U_load(1:N, 1:N)
         ! Load reference RK solution
         write(oname,'("data_GL_lyapconv_X_RK_n",I4.4,"_r",I3.3,".npy")') nx, nrep
         inquire(file=trim(basepath)//trim(oname), exist=existfile)
         if (existfile) then
            print *, 'Load RK solution ', trim(basepath)//trim(oname)
            call load_npy(trim(basepath)//trim(oname), U_load, iostatus)
            if (iostatus /= 0) then; write(*,*) "Error loading file", trim(basepath)//trim(oname); STOP 2; end if
         else
            write(*,*) 'Cannot find ', trim(basepath)//trim(oname); STOP 12
         end if
         ! Set reference RK solution
         X_mat_ref = U_load(1:N, 1:N)
      else
         ! Set random initial condition
         call get_state(U0_mat(:,1:rk_X0), U0(1:rk))
         X_out = matmul( U0_mat, matmul( S0, transpose(U0_mat ) ) )
         if (if_save_npy) then
            ! Save forcing RK
            write(oname,'("data_GL_lyapconv_BBTW_RK_n",I4.4,".npy")') nx
            write(*,*) 'Save ', trim(basepath)//trim(oname)
            call save_npy(trim(basepath)//oname, reshape(BBTW_flat, (/ N,N /)), iostatus)
            if (iostatus /= 0) then; write(*,*) "Error saving file", trim(basepath)//trim(oname); STOP 2; end if
            ! Save initial condition
            write(oname,'("data_GL_lyapconv_X0_RK_n",I4.4,".npy")') nx
            write(*,*) 'Save ', trim(basepath)//trim(oname)
            call save_npy(trim(basepath)//oname, X_out, iostatus)
            if (iostatus /= 0) then; write(*,*) "Error saving file", trim(basepath)//trim(oname); STOP 2; end if
            ! Save forcing DLRA
            Bmat = 0.0_wp
            call get_state(Bmat, LTI%B(1:rk_b))
            X_out = matmul(Bmat, transpose(Bmat))
            write(oname,'("data_GL_lyapconv_BBTW_DLRA_n",I4.4,".npy")') nx
            write(*,*) 'Save ', trim(basepath)//trim(oname)
            call save_npy(trim(basepath)//oname, X_out, iostatus)
            if (iostatus /= 0) then; write(*,*) "Error saving file", trim(basepath)//trim(oname); STOP 2; end if
         end if
         X_out = matmul( U0_mat, matmul( S0, transpose(U0_mat ) ) )
         ! Set initial condition for RK
         call set_state(X_mat(1:1), X_out)
         write(*,'(A10,A26,A26,A26,A20)') 'RKlib:','Tend','|| X_RK ||_2/N', '|| res ||_2/N','Elapsed time'
         write(*,*) '         ------------------------------------------------------------------------'
         do irep = 1, nrep
            call system_clock(count=clock_start)     ! Start Timer
            ! integrate
            call RK_propagator%matvec(X_mat(1), X_mat(2))
            call system_clock(count=clock_stop)      ! Stop Timer
            ! recover output
            call get_state(X_RKlib(:,:,irep), X_mat(2:2))
            call CALE(res_flat, reshape(X_RKlib(:,:,irep), shape(res_flat)), BBTW_flat, .false.)
            ! replace input
            call set_state(X_mat(1:1), X_RKlib(:,:,irep))
            write(*,'(I10,F26.4,E26.8,E26.8,F18.4," s")') irep, irep*Tend, norm2(X_RKlib(:,:,irep))/N, norm2(res_flat)/N, &
                           & real(clock_stop-clock_start)/real(clock_rate)
            if (if_save_npy) then
               write(oname,'("data_GL_lyapconv_X_RK_n",I4.4,"_r",I3.3,".npy")') nx, irep
               write(*,*) 'Save ', trim(basepath)//trim(oname)
               call save_npy(trim(basepath)//oname, X_RKlib(:,:,irep), iostatus)
               if (iostatus /= 0) then; write(*,*) "Error saving file", trim(basepath)//trim(oname); STOP 2; end if
            end if
         enddo
         ! Set reference RK solution
         X_mat_ref = X_RKlib(:,:,nrep)
      end if
      
      write(*,*) ''
      write(*,'(A16,A4,A4,A10,A6,A8,A26,A20)') 'DLRA:','  rk',' TO','dt','steps','Tend', &
                                 & '|| X_DLRA - X_RK ||_2/N', 'Elapsed time'
      X_state = LR_state()
      do ito = 1, size(TOv)
         torder = TOv(ito)
         do i = 1, size(rkv)
            rk = rkv(i)
            do j = 1, size(tauv)
               tau = tauv(j)
               ! Initialize low-rank representation with rank rk
               if (verb) write(*,*) 'Initialize LR state, rk =', rk
               call X_state%initialize_LR_state(U0, S0, rk)
               if (verb) write(*,*) 'Run DRLA'
               nsteps = nint(Tend/tau)
               ! run integrator
               call system_clock(count=clock_start)     ! Start Timer
               call numerical_low_rank_splitting_lyapunov_integrator(X_state, LTI%prop, LTI%B, Tend, tau, torder, info, &
                                                                     & exptA=exptA, iftrans=.false., ifverb=.false.)
               call system_clock(count=clock_stop)      ! Stop Timer
               etime = real(clock_stop-clock_start)/real(clock_rate)
               ! Reconstruct solution
               call reconstruct_solution(X_out, X_state)
               write(*,'(I4," ",A11,I4," TO",I1,F10.6,I6,F8.4,E26.8,F18.4," s")') 1, 'Xctl OUTPUT', &
                                 & rk, torder, tau, nsteps, Tend, norm2(X_out - X_mat_ref)/N, etime
               if (verb) write(*,*) 'Total integration time (DLRA):', etime, 's'
               if (j == size(tauv) .and. if_save_npy) then
                  write(onameU,'("data_GL_lyapconv_XU_n",I4.4,"_TO",I1,"_rk",I3.3,"_t",E8.2,".npy")') nx, torder, rk, tau
                  write(onameS,'("data_GL_lyapconv_XS_n",I4.4,"_TO",I1,"_rk",I3.3,"_t",E8.2,".npy")') nx, torder, rk, tau
                  write(*,*) 'Save ', trim(basepath)//trim(onameU)
                  call save_npy(trim(basepath)//onameU, U_out(:,1:rk), iostatus)
                  if (iostatus /= 0) then; write(*,*) "Error saving file", trim(onameU); STOP 2; end if
                  write(*,*) 'Save ', trim(basepath)//trim(onameS)
                  call save_npy(trim(basepath)//onameS, X_state%S(1:rk,1:rk), iostatus)
                  if (iostatus /= 0) then; write(*,*) "Error saving file", trim(onameS); STOP 2; end if
               end if
               deallocate(X_state%U)
               deallocate(X_state%S)
            end do
         end do
      end do
      
   end subroutine run_lyap_convergence_test

end module Ginzburg_Landau_Tests