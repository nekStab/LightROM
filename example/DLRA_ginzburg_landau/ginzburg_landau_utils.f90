module Ginzburg_Landau_Utils
   ! Standard Library.
   use stdlib_stats_distribution_normal, only: normal => rvs_normal
   use stdlib_strings, only: padr
   use stdlib_io_npy, only: save_npy, load_npy
   use stdlib_strings, only: replace_all
   use stdlib_optval, only : optval
   use stdlib_linalg, only : eye, diag
   !use fortime
   ! LightKrylov for linear algebra.
   use LightKrylov
   use LightKrylov_Constants, only : zero_rdp, one_rdp
   use LightKrylov_AbstractVectors
   use LightKrylov_Utils, only : assert_shape
   ! LightROM
   use LightROM_AbstractLTIsystems
   use LightROM_Utils
   use LightROM_Timing
   ! DLRA Integrators
   !use LightROM_LyapunovUtils
   !use LightROM_RiccatiUtils
   use LightROM_DLRAIntegrators
   ! Ginzburg Landau
   use Ginzburg_Landau_Base
   use Ginzburg_Landau_Operators
   use Ginzburg_Landau_RKlib

   implicit none

   character(len=*), parameter, private :: this_module = 'Ginzburg_Landau_Utils'

   integer, parameter, public :: rk_X0_lyapunov = 10
   integer, parameter, public :: rk_X0_riccati  = 10
   integer, parameter, public :: rkmax = 80
   integer, parameter, public :: irow = 8

   ! initialize mesh
   public :: initialize_parameters
   ! utilities for state_vectors
   public :: set_state, get_state, init_rand, reconstruct_solution
   ! initial conditions
   public :: generate_random_initial_condition
   ! misc
   public :: CALE, CARE
   ! readers
   public :: exist_RK_file 
   public :: exist_X_file, load_X_from_file
   ! printing helpers
   public :: print_test_header, print_test_info
   public :: print_header, print_rklib_output, get_metadata, print_dlra_output, print_svdvals
   ! saving helpers
   public :: save_RK_state_npy, save_LR_state_npy, save_metadata
   public :: make_labels, make_filename_RK, make_filename, make_filename_meta

   interface reconstruct_solution
      module procedure reconstruct_solution_X
      module procedure reconstruct_solution_US
   end interface

   interface print_svdvals
      module procedure print_svdvals_from_matrix
      module procedure print_svdvals_direct
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
      x2            = 0.0_dp
      x2(1:nx)      = x(2:nx+1)
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

      ! RK lyapunov and riccati
      Qc   = eye(rk_c)        ! state energy
      Rinv = eye(rk_b)        ! inverse control cost
      Vinv = eye(rk_c)        ! variance of sensor noise
      Qe   = eye(rk_b)        ! variance of actuator noise

      ! apply weights
      allocate(mat(N, rk_b), matW(N, rk_b))
      call get_state(mat(:,:rk_b), B(:rk_b), this_procedure)
      matW = matmul(mat, weight_mat(:rk_b,:rk_b)) ! incorporate weights
      BBTW     = matmul(mat, transpose(matW))
      BRinvBTW = matmul(mat, matmul(Rinv, transpose(matW)))
      BQeBTW   = matmul(mat, matmul(Qe, transpose(matW))) !
      deallocate(mat, matW)
      allocate(mat(N, rk_c), matW(N, rk_c))
      call get_state(mat(:,:rk_c), CT(:rk_c), this_procedure)
      matW = matmul(mat, weight_mat(:rk_c,:rk_c)) ! incorporate weights
      CTCW     = matmul(mat, transpose(matW))
      CTQcCW   = matmul(mat, matmul(Qc, transpose(matW)))
      CTVinvCW = matmul(mat, matmul(Vinv, transpose(matW)))
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

   logical function exist_X_file(fbase, load_W_) result(exist_file)
      character(len=*), intent(in) :: fbase
      logical, optional, intent(in) :: load_W_
      ! internal
      logical :: load_W
      exist_file = .false.
      load_W = optval(load_W_, .false.)
      inquire(file=trim(fbase)//'_U.npy', exist=exist_file)
      inquire(file=trim(fbase)//'_S.npy', exist=exist_file)
      if (load_W) inquire(file=trim(fbase)//'_W.npy', exist=exist_file)
   end function exist_X_file
   
   logical function exist_RK_file(fbase) result(exist_file)
      character(len=*), intent(in) :: fbase
      exist_file = .false.
      inquire(file=trim(fbase)//'_X.npy', exist=exist_file)
      inquire(file=trim(fbase)//'_meta.npy', exist=exist_file)
   end function exist_RK_file

   subroutine load_X_from_file(X, meta, fbase, U0, load_W_)
      type(LR_state), intent(inout) :: X
      real(dp), allocatable, intent(out) :: meta(:)
      character(len=*), intent(in) :: fbase
      type(state_vector), intent(in) :: U0(:) ! for type reference
      logical, optional, intent(in) :: load_W_
      ! internal
      real(dp), allocatable :: W_load(:,:)
      real(dp), allocatable :: M_load(:,:)
      logical :: load_W, exist_file
      if (allocated(X%U)) deallocate(X%U)
      if (allocated(X%S)) deallocate(X%S)

      load_W = optval(load_W_, .false.)
      ! U
      call load_npy(trim(fbase)//'_U.npy', M_load)
      X%rk = size(M_load, 2)
      allocate(X%U(X%rk), source=U0(1))
      if (load_W) then
         call load_npy(trim(fbase)//'_W.npy', W_load)
         call set_state(X%U, matmul(sqrt(W_load),M_load), 'load_X_from_file')
      else
         call set_state(X%U, M_load, 'load_X_from_file')
      end if
      ! S
      call load_npy(trim(fbase)//'_S.npy', M_load)
      allocate(X%S(X%rk,X%rk), source=M_load)
      inquire(file=trim(fbase)//'_meta.npy', exist=exist_file)
      if (exist_file) call load_npy(trim(fbase)//'_meta.npy', meta)

   end subroutine load_X_from_file

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

   function CARE(X, adjoint) result(res)
      ! solution
      real(dp)          :: X(N,N)
      ! adjoint
      logical, optional :: adjoint
      ! residual
      real(dp)          :: res(N,N)
      
      ! internals
      real(dp), dimension(N,N) :: LX, XL
      logical :: flip

      LX = 0.0_dp; XL = 0.0_dp
      flip = .not. adjoint
      call GL_mat(LX, X,            adjoint = flip, transpose = .false.)
      call GL_mat(XL, transpose(X), adjoint = flip, transpose = .true. )

      ! construct Riccati equation
      if (adjoint) then
         res = LX + XL + BQeBTW - matmul(X, matmul(CTVinvCW, X))
      else
         res = LX + XL + CTQcCW - matmul(X, matmul(BRinvBTW, X))
      end if

   end function CARE

   !-------------------------------------
   !-----      Printing helpers     -----
   !-------------------------------------

   subroutine print_test_info(if_lyapunov, if_adj, main_run, run_fixed_rank_test, run_rank_adaptive_test, run_eigenvalue_test, run_LQG_test)
      implicit none 
      logical, intent(in) :: if_lyapunov
      logical, intent(in) :: if_adj
      logical, intent(in) :: main_run
      logical, intent(in) :: run_fixed_rank_test
      logical, intent(in) :: run_rank_adaptive_test
      logical, intent(in) :: run_eigenvalue_test
      logical, intent(in) :: run_LQG_test

      print *, 'Cases to be run:'
      print '(2X,A45,L4)', padr('if_lyapunov:',50), if_lyapunov
      print '(2X,A45,L4)', padr('if_adj:',50), if_adj
      print '(2X,A45,L4)', padr('main_run:',50), main_run
      print '(2X,A45,L4)', padr('run_fixed_rank_test:',50), run_fixed_rank_test
      print '(2X,A45,L4)', padr('run_rank_adaptive_test:',50), run_rank_adaptive_test
      print '(2X,A45,L4)', padr('run_eigenvalue_test:',50), run_eigenvalue_test
      print '(2X,A45,L4)', padr('run_LQG_test:',50), run_LQG_test
   end subroutine print_test_info

   subroutine print_test_header(if_lyapunov, if_adj, eq, refid, rk_X0, Xref)
      implicit none
      logical, intent(in) :: if_lyapunov
      logical, intent(in) :: if_adj
      character(len=4), intent(out) :: eq
      character(len=2), intent(out) :: refid
      integer, intent(out) :: rk_X0
      real(dp), intent(out) :: Xref(N,N)

      ! internal
      integer :: nprint
      character(len=128) :: oname
      real(dp), allocatable :: U_load(:,:), svals(:)

      print *, ''
      call initialize_parameters()
      print *, ''

      print *, '#########################################################################'
      print *, '#                                                                       #'
      print *, '#               DYNAMIC LOW-RANK APPROXIMATION  -  DLRA                 #'
      print *, '#                                                                       #'
      print *, '#########################################################################'
      print *, ''
      print *, '             THE NON-PARALLEL LINEAR GINZBURG-LANDAU EQUATION:'
      print *, ''
      print *, '                 A = mu(x) * I + nu * D_x + gamma * D2_x'
      print *, ''
      print *, '                   with mu(x) = mu_0 * x + mu_2 * x^2'
      print *, ''
      if (if_lyapunov) then
         print *, '                     Algebraic Lyapunov equation:'
         if (if_adj) then
            print *, '                     0 = A.T @ X + X @ A + C.T @ C'
         else
            print *, '                     0 = A @ X + X @ A.T + B @ B.T'
         end if
         print *, ''               
         print *, '                   Differential Lyapunov equation:'
         if (if_adj) then
            print *, '                  \dot{X} = A.T @ X + X @ A + C.T @ C'
         else
            print *, '                  \dot{X} = A @ X + X @ A.T + B @ B.T'
         end if
         print *, ''
         print *, ''
         print '(13X,A,I4,"x",I4)', 'Complex problem size:          ', nx, nx
         print '(13X,A,I4,"x",I4)', 'Equivalent real problem size:  ', N, N
         print *, ''
         print *, '            Initial condition: rank(X0)  =', rk_X0_lyapunov
         print *, '            Inhomogeneity:     rank(B)   =', rk_B
         print *, '            Inhomogeneity:     rank(C.T) =', rk_C
         if (if_adj) then
            svals = svdvals(CTCW)
            print '(1X,A)', 'Inhomogeneity: CTCW'
            print '(1X,A,*(F16.12,X))', 'SVD(1:3)     = ', svals(1:3)
            print '(1X,A,F16.12)',      '|  CTCW  |/N = ', norm2(CTCW)/N
            oname = './example/DLRA_ginzburg_landau/CGL_Lyapunov_Observability_Yref_BS_W.npy'
         else
            svals = svdvals(BBTW)
            print '(1X,A)', 'Inhomogeneity: BBTW'
            print '(1X,A,*(F16.12,X))', 'SVD(1:3)     = ', svals(1:3)
            print '(1X,A,F16.12)',      '|  BBTW  |/N = ', norm2(BBTW)/N
            oname = './example/DLRA_ginzburg_landau/CGL_Lyapunov_Controllability_Xref_BS_W.npy'
         end if
         print *, ''
         print *, 'Check residual computation with Bartels-Stuart solution:'
         eq = 'lyap'
         refid = 'BS'
         rk_X0 = rk_X0_lyapunov
      else
         print *, '                     Algebraic Riccati equation:'
         if (if_adj) then
            print *, '     0 = A @ X + X @ A.T - X @ C.T @ V^{-1} @ C @ X + B @ Qe @ B.T'
         else
            print *, '     0 = A.T @ X + X @ A - X @ B @ R^{-1} @ B.T @ X + C.T @ Qc @ C'
         end if
         print *, ''               
         print *, '                   Differential Riccati equation:'
         if (if_adj) then
            print *, '   \dot{X} = A @ X + X @ A.T - X @ C,T @ V^{-1} @ C @ X + B @ Qe @ B.T'
         else
            print *, '   \dot{X} = A.T @ X + X @ A - X @ B @ R^{-1} @ B.T @ X + C.T @ Qc @ C'
         end if
         print *, ''
         print '(13X,A,I4,"x",I4)', 'Complex problem size:                       ', nx, nx
         print '(13X,A,I4,"x",I4)', 'Equivalent real problem size:               ', N, N
         print *, ''
         print *, '            Initial condition: rank(X0)               =', rk_X0_riccati
         if (if_adj) then
            print *, '            Nonlinearity:      rank(C.T @ V^{-1} @ C) =', rk_C
            print *, '            Inhomogeneity:     rank(B @ Qe @ B.T)     =', rk_B
            svals = svdvals(BQeBTW)
            print '(1X,A)', 'Inhomogeneity: BQeBTW'
            print '(1X,A,*(F16.12,X))', 'SVD(1:3)         = ', svals(1:3)
            print '(1X,A,F16.12)',      '|  BQeBTW  |/N   = ', norm2(BQeBTW)/N
            print *, ''
            svals = svdvals(CTVinvCW)
            print '(1X,A)', 'Nonlinearity: CTVinvCW'
            print '(1X,A,*(F16.12,X))', 'SVD(1:3)         = ', svals(1:3)
            print '(1X,A,F16.12)',      '|  CTVinvCW  |/N = ', norm2(CTVinvCW)/N
            print *, ''
            oname = './example/DLRA_ginzburg_landau/CGL_Riccati_Pref_Schur_Adjoint_W.npy'
         else
            print *, '            Nonlinearity:      rank(B @ R^{-1} @ B.T) =', rk_B
            print *, '            Inhomogeneity:     rank(C.T @ Qc @ C)     =', rk_C
            svals = svdvals(CTQcCW)
            print '(1X,A)', 'Inhomogeneity: CTQcCW'
            print '(1X,A,*(F16.12,X))', 'SVD(1:3)         = ', svals(1:3)
            print '(1X,A,F16.12)',      '|  CTQcCW  |/N   = ', norm2(CTQcCW)/N
            print *, ''
            svals = svdvals(BRinvBTW)
            print '(1X,A)', 'Nonlinearity: BRinvBTW'
            print '(1X,A,*(F16.12,X))', 'SVD(1:3)         = ', svals(1:3)
            print '(1X,A,F16.12)',      '|  BRinvBTW  |/N = ', norm2(BRinvBTW)/N
            print *, ''
            oname = './example/DLRA_ginzburg_landau/CGL_Riccati_Pref_Schur_Direct_W.npy'
         end if
         eq = 'ricc'
         refid = 'SD'
         rk_X0 = rk_X0_riccati
      end if
      print *, ''
      print *, '#########################################################################'
      print *, ''
      call load_npy(oname, U_load)
      Xref = U_load
      print *, ''
      if (if_lyapunov) then
         print '(A,A,A,F16.12)', '  |  X_', refid, '  |/N = ', norm2(Xref)/N
         print '(A,A,A,F16.12)', '  | res_', refid, ' |/N = ', norm2(CALE(Xref, if_adj))/N
      else
         print '(A,A,A,F16.12)', '  |  X_', refid, '  |/N = ', norm2(Xref)/N
         print '(A,A,A,F16.12)', '  | res_', refid, ' |/N = ', norm2(CARE(Xref, if_adj))/N
      end if
      print *, ''
      ! compute svd
      nprint = 60
      call print_svdvals(Xref, 'X_'//refid, nprint)
   end subroutine print_test_header

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
         write(*,*) repeat('-', 86)
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
      logical,           intent(in) :: adjoint
      ! internal
      integer :: N
      real(dp) :: nrm_x, nrm_diff, nrm_res

      N = size(Xref,1)

      nrm_x    = norm2(X_RK(:,:,irep)) / N
      nrm_diff = norm2(X_RK(:,:,irep) - Xref) / N
      if (trim(eq) == 'lyap') then
         nrm_res  = norm2(CALE(X_RK(:,:,irep), adjoint)) / N
      else
         nrm_res  = norm2(CARE(X_RK(:,:,irep), adjoint)) / N
      end if

      write(*,'(I7,F10.4,3(1X,E18.6),F10.4," s",A)') &
           irep, irep*Tstep, nrm_x, nrm_diff, nrm_res, etime, trim(note)

   end subroutine print_rklib_output

   subroutine get_metadata(meta, eq, rk, torder, tau, nsteps, tend, X_out, Xref_RK, Xref, etime, adjoint)
      real(dp), allocatable, intent(out) :: meta(:)
      character(len=*),  intent(in) :: eq
      integer,           intent(in) :: rk, torder, nsteps
      real(dp),          intent(in) :: tau, tend, etime
      real(dp),          intent(in) :: X_out(:,:), Xref_RK(:,:), Xref(:,:)
      logical,           intent(in) :: adjoint
      ! internal
      integer :: N
      real(dp) :: nrm_x, nrm_rk, nrm_bs, nrm_res

      if (allocated(meta)) deallocate(meta)
      N = size(X_out,1)

      nrm_x   = norm2(X_out) / N
      nrm_rk  = norm2(X_out - Xref_RK) / N
      nrm_bs  = norm2(X_out - Xref) / N
      if (trim(eq) == 'lyap') then
         nrm_res = norm2(CALE(X_out, adjoint)) / N
      else
         nrm_res = norm2(CARE(X_out, adjoint)) / N
      end if
      meta = [ etime, 1.0_dp*rk, 1.0_dp*torder, 1.0_dp*nsteps, tau, Tend, nrm_X, nrm_rk, nrm_bs, nrm_res ]
   end subroutine get_metadata

   subroutine print_dlra_output(eq, note, meta)
      character(len=*),  intent(in) :: eq
      character(len=*),  intent(in) :: note
      real(dp),          intent(in) :: meta(:)
      
      write(*,'(I4," ",A4,1X,A6,I8," TO",I1,F10.6,I8,F10.4,4(E19.8),F10.2," s")') &
         int(meta(3)), note, 'OUTPUT', int(meta(2)), int(meta(3)), meta(5), int(meta(4)), meta(6), &
         meta(7), meta(8), meta(9), meta(10), meta(1)

   end subroutine print_dlra_output

   subroutine print_svdvals_direct(svals, label, nprint)
      real(dp),          intent(in) :: svals(:)
      character(len=*),  intent(in) :: label
      integer,           intent(in) :: nprint
      ! internal
      integer :: k, j, is, ie, nblocks
   
      nblocks = ceiling(nprint * one_rdp / irow)
   
      do k = 1, nblocks
         is = (k-1)*irow + 1
         ie = min(k*irow, size(svals))
         print '(1X,A,I2,A,I2,*(1X,F16.12))', 'SVD('//label//') ', is, '-', ie, ( svals(j), j = is, ie )
      end do
      print *, ''
   
   end subroutine print_svdvals_direct

   subroutine print_svdvals_from_matrix(A, label, nprint)
      real(dp),          intent(in) :: A(:,:)
      character(len=*),  intent(in) :: label
      integer,           intent(in) :: nprint
      ! internal
      real(dp), allocatable :: svals(:)
      integer :: k, j, is, ie, nblocks
   
      svals = svdvals(A)
   
      nblocks = ceiling(nprint * one_rdp / irow)
   
      do k = 1, nblocks
         is = (k-1)*irow + 1
         ie = min(k*irow, size(svals))
         print '(1X,A,I2,A,I2,*(1X,F16.12))', 'SVD('//label//') ', is, '-', ie, ( svals(j), j = is, ie )
      end do
      print *, ''
   
   end subroutine print_svdvals_from_matrix

   !-----------------------------------
   !-----      Saving helpers     -----
   !-----------------------------------

   subroutine make_labels(Tstr, taustr, tolstr, Tend, tau, tol)
      implicit none
      character(len=32), intent(out)  :: Tstr, taustr, tolstr
      real(dp), intent(in)            :: Tend, tau, tol

      character(len=256) :: name
      
      write(Tstr,'(I3.3)') int(Tend)
      write(taustr,'(ES10.0)') tau
      taustr = adjustl(trim(taustr))
      taustr = replace_all(taustr, 'E', 'e')
      taustr = replace_all(taustr, '.', '')
      write(tolstr,'(ES10.0)') tol
      tolstr = adjustl(trim(tolstr))
      tolstr = replace_all(tolstr, 'E', 'e')
      tolstr = replace_all(tolstr, '.', '')
   end subroutine make_labels

   pure function make_filename_RK(fldr, fbase, eq, Tend, adjoint) result(name)
      implicit none
      character(len=*), intent(in)           :: fldr
      character(len=*), intent(in)           :: fbase
      character(len=4), intent(in)           :: eq
      real(dp), intent(in)                   :: Tend
      logical, intent(in)                    :: adjoint

      character(len=256) :: name
      character(len=32)  :: Tstr, mode

      write(Tstr,'(I3.3)') int(Tend)
      if (adjoint) then
         mode = '_adjoint'
      else
         mode = '_direct'
      end if
      write(name,'(A,A,A,A,"_Tend",A)') trim(fldr), trim(fbase), eq, trim(mode), Tstr
   end function make_filename_RK

   subroutine save_RK_state_npy(bname, X_mat, meta)
      implicit none
      character(len=*), intent(in) :: bname
      real(dp),         intent(in) :: X_mat(:,:)
      real(dp),         intent(in) :: meta(:)
      ! internal
      character(len=*), parameter :: this_procedure = 'save_RK_state_npy'
   
      ! Save low-rank components
      call save_npy(trim(bname)//'_X.npy', X_mat)
      ! Save weight matrix
      call save_npy(trim(bname)//'_meta.npy', meta)
   
   end subroutine save_RK_state_npy

   function make_filename(fldr, case, eq, note, rk, TO, tau, Tend, tol) result(name)
      implicit none
      character(len=*), intent(in)           :: fldr
      character(len=*), intent(in)           :: case
      character(len=*), intent(in)           :: eq
      character(len=*), intent(in)           :: note
      integer,          intent(in)           :: rk, TO
      real(dp),         intent(in)           :: tau, Tend
      real(dp),         intent(in), optional :: tol

      character(len=256) :: name
      character(len=32)  :: Tstr, tolstr, taustr

      if (present(tol)) then
         call make_labels(Tstr, taustr, tolstr, Tend, tau, tol)
         if (len(trim(note)) == 0) then
            write(name,'(A,A,"_",A,"_Tend",A,"_TO",I1,"_tau",A,"_tol",A)') &
               trim(fldr), trim(case), trim(eq), trim(Tstr), TO, trim(taustr), trim(tolstr)
         else
            write(name,'(A,A,"_",A,"_Tend",A,"_TO",I1,"_tau",A,"_tol",A,"_",A)') &
               trim(fldr), trim(case), trim(eq), trim(Tstr), TO, trim(taustr), trim(tolstr), trim(note)
         end if
      else
         call make_labels(Tstr, taustr, tolstr, Tend, tau, 0.0_dp)
         if (len(trim(note)) == 0) then
            write(name,'(A,A,"_",A,"_Tend",A,"_rk",I3.3,"_TO",I1,"_tau",A)') &
               trim(fldr), trim(case), trim(eq), trim(Tstr), rk, TO, trim(taustr)
         else
            write(name,'(A,A,"_",A,"_Tend",A,"_rk",I3.3,"_TO",I1,"_tau",A,"_",A)') &
               trim(fldr), trim(case), trim(eq), trim(Tstr), rk, TO, trim(taustr), trim(note)
         end if
      end if
   end function make_filename

   subroutine save_LR_state_npy(bname, X, weight_mat, meta)
      implicit none
      character(len=*), intent(in) :: bname
      type(LR_state),   intent(in) :: X
      real(dp),         intent(in) :: weight_mat(:,:)
      real(dp),         intent(in) :: meta(:)
      ! internal
      character(len=*), parameter :: this_procedure = 'save_LR_state_npy'
      real(dp), allocatable :: Uwrk(:,:)
      integer :: rk
   
      ! Save low-rank components
      rk = max(1, X%rk)
      allocate(Uwrk(N, rk), source=zero_rdp)
      call get_state(Uwrk, X%U(:rk), this_procedure)
      call save_npy(trim(bname)//'_U.npy', Uwrk)
      call save_npy(trim(bname)//'_S.npy', X%S(1:rk,1:rk))   
      ! Save weight matrix
      call save_npy(trim(bname)//'_W.npy', weight_mat)
      call save_npy(trim(bname)//'_meta.npy', meta)
   
   end subroutine save_LR_state_npy
   
   subroutine write_int(unit, key, value)
      integer, intent(in) :: unit, value
      character(*), intent(in) :: key
      write(unit,'(A," = ",I0)') trim(key), value
   end subroutine
  
   subroutine write_real(unit, key, value)
      integer, intent(in) :: unit
      real(dp), intent(in) :: value
      character(*), intent(in) :: key
      write(unit,'(A," = ",ES16.9)') trim(key), value
   end subroutine

   pure function make_filename_meta(fbase) result(name)
      implicit none
      character(len=*), intent(in)           :: fbase
      character(len=256) :: name
      name = trim(fbase)//'_meta.txt'
   end function make_filename_meta

   subroutine save_metadata(fbase, case, eq, rk, torder, tau, nsteps, Tend, adjoint)
      implicit none
      character(len=*), intent(in) :: fbase
      character(len=*), intent(in) :: case
      character(len=*), intent(in) :: eq
      integer,          intent(in) :: rk, torder, nsteps
      real(dp),         intent(in) :: tau, tend
      logical,          intent(in) :: adjoint
      ! internal
      character(len=256) :: fname
      character(len=20) :: str
      integer :: i, unit, n_called
      real(dp) :: etime, etmin, etmax, etimp
      integer :: lcount, rcount, gcount
      character(len=128), allocatable :: names(:)

      ! write to file
      fname = make_filename_meta(fbase)
      open(newunit=unit, file=fname, status='replace', action='write')
      
      write(unit, '(A20,A)')     padr('case :',20), case
      write(unit, '(A20,A)')     padr('equation :',20), eq
      write(unit, '(A20,I16)')   padr('torder :',20), torder
      str = merge('rk(adapt)', 'rk(fixed)', case == 'DLRA_ADAPT')
      write(unit, '(A20,I16)')   padr(str,20), rk
      write(unit, '(A20,F16.9)') padr('tau :',20), tau
      write(unit, '(A20,F16.9)') padr('Tend :',20), Tend
      write(unit, '(A20,L16)')   padr('adjoint :',20), adjoint
      ! get timer that were called since last reset
      call lk_timer%get_called(n_called, names)
      write(unit, '(A50,2X,A12,*(1X,A13))') padr('TIMERS:',50), 'count', 'etime', 'min', 'max', 'paused', 'avg'
      ! get timer info
      do i = 1, n_called
         call lk_timer%get_data(names(i), restart=.false., &
                     & etime=etime, etmin=etmin, etmax=etmax, etimp=etimp, &
                     & lcount=lcount, rcount=rcount, gcount=gcount)
         
         write(unit, '(A5,A45,2X,I12,*(1X,F13.6)))') 'LK % ', padr(trim(names(i)),45), lcount, etime, etmin, etmax, etimp, etime/lcount
      end do
      call global_lightROM_timer%get_called(n_called, names)
      ! get timer info
      do i = 1, n_called
         call global_lightROM_timer%get_data(names(i), restart=.false., &
                     & etime=etime, etmin=etmin, etmax=etmax, etimp=etimp, &
                     & lcount=lcount, rcount=rcount, gcount=gcount)
         
         write(unit, '(A5,A45,2X,I12,*(1X,F13.6)))') 'LR % ', padr(trim(names(i)),45), lcount, etime, etmin, etmax, etimp, etime/lcount
      end do
      close(unit)
   end subroutine save_metadata

end module Ginzburg_Landau_Utils