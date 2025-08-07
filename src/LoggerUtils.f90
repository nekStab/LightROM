module LightROM_LoggerUtils
   ! stdlib
   use stdlib_strings, only: padl
   use stdlib_optval, only : optval
   ! LightKrylov
   use LightKrylov, only : dp, wp => dp
   use LightKrylov_Constants, only: io_rank
   use LightKrylov_Logger, only: log_message, log_information, log_warning, log_debug, check_info, stop_error
   ! LightROM
   use LightROM_AbstractLTIsystems, only: abstract_sym_low_rank_state_rdp
   use LightROM_Utils, only: dlra_opts
   
   implicit none 

   private :: this_module
   character(len=*), parameter :: this_module = 'LR_LoggerUtils'
   integer, parameter :: iline = 4  ! number of svals per line
   logical :: if_overwrite = .true.
   integer :: rename_counter = 0

   public :: write_logfile_headers, reset_logfiles, stamp_logfiles
   public :: log_settings, log_svals

contains

   subroutine write_logfile_headers(basename)
      character(*), intent(in) :: basename
      ! parameters
      character(len=*), parameter :: this_procedure = 'write_logfile_headers'
      ! internal
      integer :: i
      character(len=128) :: fname, suffix, info
      if (io_rank() .and. if_overwrite) then
         do i = 1, 2
            select case (i)
               case (1)
                  suffix = '_abs.dat'
                  info   = 'singular value'
               case (2)
                  suffix = '_rel.dat'
                  info   = 'singular value relative difference'
            end select
            write(fname,'(A,A)') trim(basename), trim(suffix)
            open (1234, file=fname, status='replace', action='write')
            write (1234, '(A8,A8,2(A15,1X),A4,4X,A)') 'icall', 'istep', 'time', 'lag', 'rk', trim(info)
            close (1234)
         end do
      end if
      if_overwrite = .false.
      return
   end subroutine write_logfile_headers

   subroutine stamp_logfiles(basename, X, lag, svals, dsvals, icall)
      character(*), intent(in) :: basename
      class(abstract_sym_low_rank_state_rdp),  intent(in) :: X
      real(dp), intent(in) :: lag
      real(dp), dimension(:), intent(in) :: svals
      real(dp), dimension(:), intent(in) :: dsvals
      integer, intent(in) :: icall
      ! parameters
      character(len=*), parameter :: this_procedure = 'stamp_logfiles'
      ! internal
      integer :: i
      character(len=128) :: fname, suffix
      if (io_rank()) then
         do i = 1, 2
            select case (i)
               case (1)
                  suffix = '_abs.dat'
                  write(fname,'(A,A)') trim(basename), trim(suffix)
                  ! SVD absolute
                  open (1234, file=fname, status='old', action='write', position='append')
                  write (1234, '(I8,1X,I7,2(1X,F15.9),I4)', ADVANCE='NO') icall, X%tot_step, X%tot_time, lag, X%rk
                  write (1234, '(*(1X,F15.9))') svals
                  close (1234)
               case (2)
                  suffix = '_rel.dat'
                  write(fname,'(A,A)') trim(basename), trim(suffix)
                  ! dSVD relative
                  open (1234, file=fname, status='old', action='write', position='append')
                  write (1234, '(I8,1X,I7,2(1X,F15.9),I4)', ADVANCE='NO') icall, X%tot_step, X%tot_time, lag, X%rk
                  write (1234, '(*(1X,F15.9))') dsvals
                  close (1234)
            end select 
         end do
      end if
      return
   end subroutine stamp_logfiles

   subroutine reset_logfiles(basename, if_rename)
      character(*), intent(in) :: basename
      logical, optional, intent(in) :: if_rename
      ! parameters
      character(len=*), parameter :: this_procedure = 'reset_logfiles'
      ! internal
      integer :: i
      logical :: rename_logfiles, exist_origin, exist_target
      character(len=128) :: fname, fname_new, suffix, msg
      if_overwrite = .true.
      rename_logfiles = optval(if_rename, .true.)
      if (rename_logfiles) then
         rename_counter = rename_counter + 1
         do i = 1, 2
            select case (i)
               case (1)
                  suffix = '_abs.dat'
               case (2)
                  suffix = '_rel.dat'
            end select
            write(fname,    '(A,A)')      trim(basename),                 trim(suffix)
            write(fname_new,'(A,I3.3,A)') trim(basename), rename_counter, trim(suffix)
            inquire(file=fname,     exist=exist_origin)
            inquire(file=fname_new, exist=exist_target)
            if (exist_origin) then
               if (exist_target) then
                  msg = trim(fname_new)//' exists and will be overwritten.'
                  call log_message(msg, this_module, this_procedure)
               end if
               msg = 'Renaming '//trim(fname)//' --> '//trim(fname_new)
               call log_message(msg, this_module, this_procedure)
               call rename(fname, fname_new)
            end if
         end do
      else
         msg = 'Logfiles not renamed. Files may be overwritten.'
         call log_warning(msg, this_module, this_procedure)
      end if
   end subroutine reset_logfiles

   subroutine log_settings(X, Tend, tau, nsteps, opts)
      class(abstract_sym_low_rank_state_rdp),  intent(in) :: X
      real(dp), intent(in) :: Tend
      real(dp), intent(in) :: tau
      integer, intent(in) :: nsteps
      type(dlra_opts), intent(in) :: opts
      ! parameters
      character(len=*), parameter :: this_procedure = 'log_settings'
      ! internals
      character(len=128) :: msg, ctype

      call log_message('###### solver settings ######', this_module, this_procedure)
      write(msg,'(A15," : ", F15.8)') padl('t0',15), X%tot_time
      call log_message(msg, this_module, this_procedure)
      write(msg,'(A15," : ", F15.8)') padl('tf',15), X%tot_time + Tend
      call log_message(msg, this_module, this_procedure)
      write(msg,'(A15," : ", F15.8)') padl('dt',15), tau
      call log_message(msg, this_module, this_procedure)
      write(msg,'(A15," : ", I0)')    padl('nsteps',15),  nsteps
      call log_message(msg, this_module, this_procedure)
      write(msg,'(A15," : ", I0)')    padl('t-order',15), opts%mode
      call log_message(msg, this_module, this_procedure)
      write(msg,'(A15," : ", L)')     padl('adaptive rank',15), opts%if_rank_adaptive
      call log_message(msg, this_module, this_procedure)
      if (opts%if_rank_adaptive) then
         write(msg,'(A15," : ", I0)') padl('rk_init',15), X%rk
         call log_message(msg, this_module, this_procedure)
         write(msg,'(A15," : ", I0)') padl('rk_max',15), size(X%U)
         call log_message(msg, this_module, this_procedure)
         write(msg,'(A15," : sigma_{r+1} < ", E15.8)') padl('adapt. tol.',15), opts%tol
         call log_message(msg, this_module, this_procedure)
      else
         write(msg,'(A15," : ", I0)') padl('rk',15), X%rk
         call log_message(msg, this_module, this_procedure)
      end if
      if (opts%relative_inc) then
         ctype = 'relative'
      else
         ctype = 'absolute'
      end if
      write(msg,'(A15," : ",A,A)')    padl('convergence',15), trim(ctype), ' increment of the solution 2-norm'
      call log_message(msg, this_module, this_procedure)
      write(msg,'(A15," : ", E15.8)') padl('tol',15), opts%inc_tol
      call log_message(msg, this_module, this_procedure)
      if (opts%chkctrl_time) then
         write(msg,'("  Output every ",F8.4," time units (",I0," steps)")') opts%chktime, nint(opts%chktime/tau)
         call log_message(msg, this_module, this_procedure)
      else
         write(msg,'("  Output every ",I0," steps (",F8.4," time units)")') opts%chkstep, opts%chkstep*tau
         call log_message(msg, this_module, this_procedure)
      end if
      call log_message('###### solver settings ######', this_module, this_procedure)
   end subroutine log_settings

   subroutine log_step(X, istep, nsteps)
      class(abstract_sym_low_rank_state_rdp),  intent(in) :: X
      integer, intent(in) :: istep
      integer, intent(in) :: nsteps
      ! parameters
      character(len=*), parameter :: this_procedure = 'log_step'
      ! internal
      integer :: irk, ifmt
      character(len=128) :: msg, fmt
      ifmt = max(5,ceiling(log10(real(nsteps))))
      write(fmt,'(A,2(I0,A))') '("Step ",I', ifmt, ',"/",I', ifmt, ',": T= ",F10.4,", Ttot= ",F10.4)'
      write(msg,fmt) istep, nsteps, X%time, X%tot_time
      call log_information(msg, this_module, 'DLRA_main')
   end subroutine log_step

   subroutine log_svals(basename, X, lag, svals, svals_lag, icall, istep, nsteps)
      character(*), intent(in) :: basename
      class(abstract_sym_low_rank_state_rdp),  intent(in) :: X
      real(dp), intent(in) :: lag
      real(dp), dimension(:), intent(in) :: svals
      real(dp), dimension(:), intent(in) :: svals_lag
      integer, intent(in) :: icall
      integer, intent(in) :: istep
      integer, intent(in) :: nsteps
      ! parameters
      character(len=*), parameter :: this_procedure = 'log_svals'
      ! internal
      integer :: i, j, is, ie, irk, ifmt, irkfmt
      real(wp), dimension(:), allocatable :: dsvals
      character(len=128) :: msg, fmt

      irk = min(size(svals), size(svals_lag))
      allocate(dsvals(irk)); dsvals = 0.0_dp
      ifmt = max(5,ceiling(log10(real(nsteps))))
      irkfmt = max(3,ceiling(log10(real(size(X%U)))))
      write(fmt,'(A,5(I0,A))') '("Step ",I', ifmt, ',"/",I', ifmt, ',": T= ",F10.4,1X,I', irkfmt, '" : ",A,"[",I', irkfmt, ',"-",I', irkfmt, ',"]",*(E12.5))'
      do i = 1, irk
         dsvals(i) = abs(svals(i)-svals_lag(i))/svals(i)
      end do
      do i = 1, ceiling(float(X%rk)/iline)
         is = (i-1)*iline+1; ie = min(X%rk,i*iline)
         write(msg,fmt) istep, nsteps, X%tot_time, X%rk, " SVD abs", is, ie, ( svals(j), j = is, ie )
         call log_information(msg, this_module, 'DLRA_main')
      end do
      do i = 1, ceiling(float(irk)/iline)
         is = (i-1)*iline+1; ie = min(irk,i*iline)
         write(msg,fmt) istep, nsteps, X%tot_time, X%rk, "dSVD rel", is, ie, ( dsvals(j) , j = is, ie )
         call log_information(msg, this_module, 'DLRA_main')
      end do
      ! Add information to the SVD logfile as well
      call stamp_logfiles(basename, X, lag, svals, dsvals, icall)
   end subroutine log_svals

end module LightROM_LoggerUtils
