module Ginzburg_Landau_Utils
   ! Standard Library.
   use stdlib_stats_distribution_normal, only: normal => rvs_normal
   use stdlib_io_npy, only: save_npy
   use stdlib_strings, only: replace_all
   use stdlib_optval, only : optval
   use stdlib_linalg, only : eye, diag, svd
   !use fortime
   ! LightKrylov for linear algebra.
   use LightKrylov
   use LightKrylov_Constants, only : zero_rdp
   use LightKrylov_AbstractVectors
   use LightKrylov_Utils, only : assert_shape
   ! LightROM
   use LightROM_AbstractLTIsystems
   use LightROM_Utils
   ! Lyapunov Solver
   use LightROM_LyapunovSolvers
   use LightROM_LyapunovUtils
   ! Riccati Solver
   use LightROM_RiccatiSolvers
   use LightROM_RiccatiUtils
   ! Ginzburg Landau
   use Ginzburg_Landau_Base
   use Ginzburg_Landau_Operators
   use Ginzburg_Landau_RKlib

   implicit none

   character(len=*), parameter, private :: this_module = 'Ginzburg_Landau_Utils'

   ! initialize mesh
   public :: initialize_parameters
   ! utilities for state_vectors
   public :: set_state, get_state, init_rand, reconstruct_solution
   ! initial conditions
   public :: generate_random_initial_condition
   ! misc
   public :: CALE, CARE
   ! printing helpers
   public :: print_header, print_rklib_output, print_dlra_output, print_svdvals
   ! saving helpers
   public :: save_LR_state_npy, make_filename

   

   interface reconstruct_solution
      module procedure reconstruct_solution_X
      module procedure reconstruct_solution_US
   end interface

contains

   !--------------------------------------------------------------
   !-----     CONSTRUCT THE MESH AND PHYSICAL PARAMETERS     -----
   !--------------------------------------------------------------

   subroutine initialize_parameters()
      implicit none
      character(len=*), parameter :: this_procedure = 'initialize_parameters'
      ! Mesh array.
      real(dp), allocatable :: x(:)
      real(dp)              :: x2(1:2*nx)
      real(dp), allocatable :: mat(:,:), matW(:,:)
      integer               :: i

      ! Construct mesh.
      x = linspace(-L/2, L/2, nx+2)
      dx = x(2)-x(1)

      ! Construct mu(x)
      mu(:) = (mu_0 - c_mu**2) + (mu_2 / 2.0_dp) * x(2:nx+1)**2

      ! Define integration weights
      weight          = dx
      weight_mat      = eye(N)*dx
      weight_flat     = dx

      ! Construct B & C
      ! B = [ [ Br, -Bi ], [ Bi, Br ] ]
      ! B = [ [ Cr, -Ci ], [ Ci, Cr ] ]
      ! where Bi = Ci = 0

      ! actuator is a Guassian centered just upstream of branch I
      ! column 1
      x2       = 0.0_dp
      x2(1:nx) = x(2:nx+1)
      B(1)%state = exp(-((x2 - x_b)/s_b)**2)
      ! column 2
      x2            = 0.0_dp
      x2(nx+1:2*nx) = x(2:nx+1)
      B(2)%state = exp(-((x2 - x_b)/s_b)**2)

      ! the sensor is a Gaussian centered at branch II
      ! column 1
      x2       = 0.0_dp
      x2(1:nx) = x(2:nx+1)
      CT(1)%state = exp(-((x2 - x_c)/s_c)**2)
      ! column 2
      x2            = 0.0_dp
      x2(nx+1:2*nx) = x(2:nx+1)
      CT(2)%state = exp(-((x2 - x_c)/s_c)**2)

      ! RK lyap & riccati
      Qc   = eye(rk_c)
      Rinv = eye(rk_b)

      allocate(mat(N, rk_b), matW(N, rk_b))
      call get_state(mat(:,:rk_b), B(:rk_b), this_procedure)
      matW = matmul(mat, weight_mat(:rk_b,:rk_b)) ! incorporate weights
      BBTW = matmul(mat, transpose(matW))
      BRinvBTW  = matmul(mat, matmul(Rinv, transpose(matW)))
      deallocate(mat, matW)

      allocate(mat(N, rk_c), matW(N, rk_c))
      call get_state(mat(:,:rk_c), CT(:rk_c), this_procedure)
      matW = matmul(mat, weight_mat(:rk_c,:rk_c)) ! incorporate weights
      CTCW = matmul(mat, transpose(matW))
      CTQcCW =  matmul(mat, matmul(Qc, transpose(matW)))
      deallocate(mat, matW)

      print '(A)', ' ----------------------------------------'
      print '(A)', '    LINEAR GINZBURG LANDAU PARAMETERS'
      print '(A)', ' ----------------------------------------'
      print '(4X,A,F10.6," + ",F10.6," i")', 'nu    = ', nu
      print '(4X,A,F10.6," + ",F10.6," i")', 'gamma = ', gamma
      print '(4X,A,F10.6)', 'mu_0  = ', mu_0
      print '(4X,A,F10.6)', 'c_mu  = ', c_mu
      print '(4X,A,F10.6)', 'mu_2  = ', mu_2
      print '(4X,A)', '-----------------------'
      print '(4X,A)', ' Inhomogeneities'
      print '(4X,A)', '-----------------------'
      print '(4X,A,I10,A)',   'rk_b  = ', rk_b, '     ! forcing rank'
      print '(4X,A,F10.6,A)', 'x_b   = ', x_b,  '     ! forcing location'
      print '(4X,A,F10.6,A)', 's_b   = ', s_b,  '     ! std dev of gaussian distribution'
      print '(4X,A,I10,A)',   'rk_c  = ', rk_c, '     ! sensor rank'
      print '(4X,A,F10.6,A)', 'x_c   = ', x_c,  '     ! sensor location'
      print '(4X,A,F10.6,A)', 's_c   = ', s_c,  '     ! std dev of gaussian distribution'
      print '(A)', ' ----------------------------------------'

   end subroutine initialize_parameters

   !--------------------------------------------------------------------
   !-----     UTILITIES FOR STATE_VECTOR AND STATE MATRIX TYPES    -----
   !--------------------------------------------------------------------

   subroutine get_state(mat_out, state_in, procedure)
      !! Utility function to transfer data from a state vector to a real array
      real(dp),                   intent(out) :: mat_out(:,:)
      class(abstract_vector_rdp), intent(in)  :: state_in(:)
      character(len=*),           intent(in)  :: procedure
      ! internal variables
      character(len=*), parameter :: this_procedure = 'get_state'
      integer :: k, kdim
      mat_out = 0.0_dp
      select type (state_in)
      type is (state_vector)
         kdim = size(state_in)
         call assert_shape(mat_out, [ N, kdim ], 'mat_out', this_module, this_procedure//': '//trim(procedure))
         do k = 1, kdim
            mat_out(:,k) = state_in(k)%state
         end do
      type is (state_matrix)
         call assert_shape(mat_out, [ N, N ], 'mat_out', this_module, this_procedure//': '//trim(procedure))
         mat_out = reshape(state_in(1)%state, [ N, N ])
      class default
         call type_error('state', 'state_vector or state_matrix', 'IN', this_module, this_procedure)
      end select
   end subroutine get_state

   subroutine set_state(state_out, mat_in, procedure)
      !! Utility function to transfer data from a real array to a state vector
      class(abstract_vector_rdp), intent(out) :: state_out(:)
      real(dp),                   intent(in)  :: mat_in(:,:)
      character(len=*),           intent(in)  :: procedure
      ! internal variables
      character(len=*), parameter :: this_procedure = 'set_state'
      integer       :: k, kdim
      select type (state_out)
      type is (state_vector)
         kdim = size(state_out)
         call assert_shape(mat_in, [ N, kdim ], 'mat_in', this_module, this_procedure//': '//trim(procedure))
         call zero_basis(state_out)
         do k = 1, kdim
            state_out(k)%state = mat_in(:,k)
         end do
      type is (state_matrix)
         call assert_shape(mat_in, [ N, N ], 'mat_in', this_module, this_procedure//': '//trim(procedure))
         call zero_basis(state_out)
         state_out(1)%state = reshape(mat_in, shape(state_out(1)%state))
      class default
         call type_error('state', 'state_vector or state_matrix', 'IN', this_module, this_procedure)
      end select
   end subroutine set_state

   subroutine init_rand(state, ifnorm)
      !! Utility function to initialize a state vector with random data
      class(abstract_vector_rdp), intent(inout)  :: state(:)
      logical, optional,          intent(in)     :: ifnorm
      ! internal variables
      character(len=*), parameter :: this_procedure = 'init_rand'
      integer :: k, kdim
      logical :: normalize
      normalize = optval(ifnorm,.true.)
      select type (state)
      type is (state_vector)
         kdim = size(state)
         do k = 1, kdim
            call state(k)%rand(ifnorm = normalize)
         end do
      type is (state_matrix)
         kdim = size(state)
         do k = 1, kdim
            call state(k)%rand(ifnorm = normalize)
         end do
      class default
         call type_error('state', 'state_vector or state_matrix', 'INOUT', this_module, this_procedure)
      end select
   end subroutine init_rand

   subroutine reconstruct_solution_X(X, LR_X)
      real(dp),          intent(out) :: X(:,:)
      type(LR_state),    intent(in)  :: LR_X
      
      ! internals
      character(len=*), parameter :: this_procedure = 'reconstruct_solution_X'
      real(dp) :: Uwrk(N, LR_X%rk)

      call assert_shape(X, [ N, N ], 'X', this_module, this_procedure)
      call get_state(Uwrk, LR_X%U(1:LR_X%rk), this_procedure)
      X = matmul(matmul(Uwrk, matmul(LR_X%S(1:LR_X%rk,1:LR_X%rk), transpose(Uwrk))), weight_mat)

   end subroutine reconstruct_solution_X

   subroutine reconstruct_solution_US(X, U, S)
      real(dp),           intent(out) :: X(:,:)
      type(state_vector), intent(in)  :: U(:)
      real(dp),           intent(in)  :: S(:,:)
      
      ! internals
      character(len=*), parameter :: this_procedure = 'reconstruct_solution_US'
      integer  :: rk
      real(dp) :: Uwrk(N, size(U))

      rk = size(U)
      call assert_shape(X, [ N, N ], 'X', this_module, this_procedure)
      call assert_shape(S, [ rk, rk ], 'S', this_module, this_procedure)
      call get_state(Uwrk, U, this_procedure)
      X = matmul(matmul(Uwrk, matmul(S, transpose(Uwrk))), weight_mat)

   end subroutine reconstruct_solution_US

   !------------------------------------
   !-----     INITIAL CONDIIONS    -----
   !------------------------------------

   subroutine generate_random_initial_condition(U, S, rk)
      class(state_vector),   intent(out) :: U(:)
      real(dp),              intent(out) :: S(:,:)
      integer,               intent(in)  :: rk
      ! internals
      character(len=*), parameter :: this_procedure = 'generate_random_initial_condition'
      integer                            :: i
      real(dp)                           :: mean, std
      character(len=128) :: msg

      ! sanity check
      if (size(U) < rk) then
         write(msg,'(A,I0)') 'Input krylov basis size incompatible with requested rank ', rk
         call stop_error(msg, this_module, this_procedure)
      end if

      call zero_basis(U)
      call initialize_random_orthonormal_basis(U(:rk))
      call assert_shape(S, [ rk,rk ], 'S', this_module, this_procedure)
      S = zero_rdp
      mean = zero_rdp
      std  = one_rdp
      do i = 1, rk
         S(i,i) = normal(mean,std)
      end do
      write(msg,'(A,I0,A,I0,A)') 'size(U) = [ ', size(U),' ]: filling the first ', rk, ' columns with orthonormal noise.'
      call logger%log_information(msg, this_module, this_procedure)
   end subroutine

   !-------------------------
   !-----      MISC     -----
   !-------------------------

   function CALE(X, adjoint) result(res)
      ! solution
      real(dp)          :: X(N,N)
      ! adjoint
      logical, optional :: adjoint
      ! residual
      real(dp)          :: res(N,N)

      ! internals
      real(dp), dimension(N,N) :: AX, XAH

      AX = 0.0_dp; XAH = 0.0_dp
      call GL_mat(AX,  X,             adjoint = adjoint, transpose = .false.)
      call GL_mat(XAH, transpose(X),  adjoint = adjoint, transpose = .true. )

      ! construct Lyapunov equation
      if (adjoint) then
         res = AX + XAH + CTCW
      else
         res = AX + XAH + BBTW
      end if

   end function CALE

   function CARE(X, CTQcCW, BRinvBTW) result(res)
      ! solution
      real(dp)          :: X(N,N)
      ! inhomogeneity
      real(dp)          :: CTQcCW(N,N)
      ! inhomogeneity
      real(dp)          :: BRinvBTW(N,N)
      ! residual
      real(dp)          :: res(N,N)
      
      ! internals
      real(dp), dimension(N,N) :: AHX, XA

      AHX = 0.0_dp; XA = 0.0_dp
      call GL_mat(AHX, X,            adjoint = .true., transpose = .false.)
      call GL_mat(XA, transpose(X),  adjoint = .true., transpose = .true. )

      ! construct Lyapunov equation
      res = AHX + XA + CTQcCW - matmul(X, matmul(BRinvBTW, X))

   end function CARE

   !-------------------------------------
   !-----      Printing helpers     -----
   !-------------------------------------

   subroutine print_header(case, eq)
      character(len=*), intent(in) :: case
      character(len=*), intent(in) :: eq
      ! internal
      character(len=4) :: ref
      ref = merge('X_BS', 'X_SD', eq=='lyap')
      select case (trim(case))
      case ('RKLIB')
         write(*,'(A7,A10,A19,A19,A19,A12)') &
              ' RKlib:', 'Tend', '| X_RK |/N', '| X_RK - '//trim(ref)//' |/N', '| res_RK |/N', 'etime'
         write(*,*) repeat('-', 85)
      case ('DLRA_FIXED')
         print '(A16,A8,A4,A10,A8,A10,4(A19),A12)', &
         'DLRA:', '  rk', ' TO', 'dt', 'steps', 'Tend', &
         '| X_D |/N', '| X_D - X_RK |/N', '| X_D - '//trim(ref)//' |/N', '| res_D |/N', 'etime'
         write(*,*) repeat('-', 144)
      case ('DLRA_ADAPT')
         print '(A16,A8,A4,A10,A8,A10,4(A19),A12)', &
         'DLRA:', 'rk_end', ' TO', 'dt', 'steps', 'Tend', &
         '| X_D |/N', '| X_D - X_RK |/N', '| X_D - '//trim(ref)//' |/N', '| res_D |/N', 'etime'
         write(*,*) repeat('-', 144)
      case default
         error stop 'Unknown header case in print_header'
      end select
   end subroutine print_header

   subroutine print_rklib_output(eq, irep, Tstep, X_RK, Xref, etime, note, adjoint)
      character(len=*),  intent(in) :: eq
      integer,           intent(in) :: irep
      real(dp),          intent(in) :: Tstep, etime
      real(dp),          intent(in) :: X_RK(:,:,:)
      real(dp),          intent(in) :: Xref(:,:)
      character(len=*),  intent(in) :: note
      logical, optional, intent(in) :: adjoint
      ! internal
      integer :: N
      real(dp) :: nrm_x, nrm_diff, nrm_res

      N = size(Xref,1)

      nrm_x    = norm2(X_RK(:,:,irep)) / N
      nrm_diff = norm2(X_RK(:,:,irep) - Xref) / N
      if (trim(eq) == 'lyap') then
         nrm_res  = norm2(CALE(X_RK(:,:,irep), adjoint)) / N
      else
         nrm_res  = norm2(CARE(X_RK(:,:,irep), CTQcCW, BRinvBTW)) / N
      end if

      write(*,'(I7,F10.4,3(1X,E18.6),F10.4," s",A)') &
           irep, irep*Tstep, nrm_x, nrm_diff, nrm_res, etime, trim(note)

   end subroutine print_rklib_output

   subroutine print_dlra_output(eq, rk, torder, tau, nsteps, tend, X_out, Xref_RK, Xref, etime, note, adjoint)
      character(len=*),  intent(in) :: eq
      integer,           intent(in) :: rk, torder, nsteps
      real(dp),          intent(in) :: tau, tend, etime
      real(dp),          intent(in) :: X_out(:,:), Xref_RK(:,:), Xref(:,:)
      character(len=*),  intent(in) :: note
      logical, optional, intent(in) :: adjoint
      ! internal
      integer :: N
      real(dp) :: nrm_x, nrm_rk, nrm_bs, nrm_res

      N = size(X_out,1)

      nrm_x   = norm2(X_out) / N
      nrm_rk  = norm2(X_out - Xref_RK) / N
      nrm_bs  = norm2(X_out - Xref) / N
      if (trim(eq) == 'lyap') then
         nrm_res = norm2(CALE(X_out, adjoint)) / N
      else
         nrm_res = norm2(CARE(X_out, CTQcCW, BRinvBTW)) / N
      end if

      write(*,'(I4," ",A4,1X,A6,I8," TO",I1,F10.6,I8,F10.4,4(E19.8),F10.2," s")') &
         1, note, 'OUTPUT', rk, torder, tau, nsteps, tend, &
         nrm_x, nrm_rk, nrm_bs, nrm_res, etime

   end subroutine print_dlra_output

   subroutine print_svdvals(A, label, nprint, irow)
      real(dp),          intent(in) :: A(:,:)
      character(len=*),  intent(in) :: label
      integer,           intent(in) :: nprint, irow
      ! internal
      real(dp), allocatable :: svals(:)
      integer :: k, j, is, ie, nblocks
      real(dp), parameter :: one_rdp = 1.0_dp
   
      svals = svdvals(A)
   
      nblocks = ceiling(nprint * one_rdp / irow)
   
      do k = 1, nblocks
         is = (k-1)*irow + 1
         ie = min(k*irow, size(svals))
         print '(1X,A,I2,A,I2,*(1X,F16.12))', 'SVD('//trim(label)//') ', is, '-', ie, ( svals(j), j = is, ie )
      end do
      print *, ''
   
   end subroutine print_svdvals

   !-----------------------------------
   !-----      Saving helpers     -----
   !-----------------------------------

   subroutine save_LR_state_npy(filename, X, weight_mat)
      
      implicit none
   
      character(len=*), intent(in) :: filename
      type(LR_state),   intent(in) :: X
      real(dp),         intent(in) :: weight_mat(:,:)
      ! internal
      character(len=*), parameter :: this_procedure = 'save_LR_state_npy'
      real(dp), allocatable :: Uwrk(:,:)
      character(len=:), allocatable :: base
   
      ! strip extension
      base = trim(filename)
      if (index(base, '.npy') > 0) base = base(:index(base,'.npy')-1)
   
      ! Save low-rank components
      allocate(Uwrk(N, X%rk), source=zero_rdp)
      call get_state(Uwrk, X%U(:X%rk), this_procedure)
      call save_npy(base//'_U.npy', Uwrk)
      call save_npy(base//'_S.npy', X%S(1:X%rk,1:X%rk))   
      ! Save weight matrix
      call save_npy(base//'_W.npy', weight_mat)
   
   end subroutine save_LR_state_npy

   pure function make_filename(fldr, case, note, rk, TO, tau, Tend, tol) result(name)
      implicit none

      character(len=*), intent(in)           :: fldr
      character(len=*), intent(in)           :: case
      character(len=*), intent(in)           :: note
      integer,          intent(in)           :: rk, TO
      real(dp),         intent(in)           :: tau, Tend
      real(dp),         intent(in), optional :: tol

      character(len=256) :: name
      character(len=32)  :: tolstr, taustr

      ! ---- tau string ----------------------------------------------
      taustr = ''
      write(tolstr,'(ES10.2)') tau
      taustr = adjustl(trim(taustr))
      taustr = replace_all(taustr, 'E', 'e')
      taustr = replace_all(taustr, '+', '')

      ! ---- tolerance string ----------------------------------------------
      tolstr = ''
      if (present(tol)) then
         write(tolstr,'(ES10.2)') tol
         tolstr = adjustl(trim(tolstr))
         tolstr = replace_all(tolstr, 'E', 'e')
         tolstr = replace_all(tolstr, '+', '')
      end if

      ! ---- assemble filename ---------------------------------------------
      if (present(tol)) then
         write(name,'(A,A,"_rk",I3.3,"_TO",I1,"_tau",A,"_Tend",I3.3,"_tol",A,"_",A,".npy")') &
            trim(fldr), trim(case), rk, TO, trim(taustr), int(Tend), trim(tolstr), trim(note)
      else
         write(name,'(A,A,"_rk",I3.3,"_TO",I1,"_tau",A,"_Tend",I3.3,"_",A,".npy")') &
            trim(fldr), trim(case), rk, TO, trim(taustr), int(Tend), trim(note)
      end if

   end function make_filename


end module Ginzburg_Landau_Utils