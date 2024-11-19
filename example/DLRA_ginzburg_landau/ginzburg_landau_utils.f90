module Ginzburg_Landau_Utils
   ! Standard Library.
   use stdlib_optval, only : optval
   use stdlib_linalg, only : eye, diag, svd
   !use fortime
   ! LightKrylov for linear algebra.
   use LightKrylov
   use LightKrylov, only : wp => dp
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
   use Ginzburg_Landau_RK_Lyapunov

   implicit none

   private :: this_module
   ! initialize mesh
   public  :: initialize_parameters
   ! utilities for state_vectors
   public  :: set_state, get_state, init_rand, reconstruct_solution
   ! initial conditions
   public  :: generate_random_initial_condition
   ! misc
   public  :: CALE, CARE

   character(len=*), parameter :: this_module = 'Ginzburg_Landau_Utils'

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
      ! Mesh array.
      real(wp), allocatable :: x(:)
      real(wp)              :: x2(1:2*nx)
      real(wp), allocatable :: mat(:,:), matW(:,:)
      integer               :: i

      ! Construct mesh.
      x = linspace(-L/2, L/2, nx+2)
      dx = x(2)-x(1)

      ! Construct mu(x)
      mu(:) = (mu_0 - c_mu**2) + (mu_2 / 2.0_wp) * x(2:nx+1)**2

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
      x2       = 0.0_wp
      x2(1:nx) = x(2:nx+1)
      B(1)%state = exp(-((x2 - x_b)/s_b)**2)!*sqrt(weight)
      ! column 2
      x2            = 0.0_wp
      x2(nx+1:2*nx) = x(2:nx+1)
      B(2)%state = exp(-((x2 - x_b)/s_b)**2)!*sqrt(weight)

      ! the sensor is a Gaussian centered at branch II
      ! column 1
      x2       = 0.0_wp
      x2(1:nx) = x(2:nx+1)
      CT(1)%state = exp(-((x2 - x_c)/s_c)**2) !/sqrt(weight)
      ! column 2
      x2            = 0.0_wp
      x2(nx+1:2*nx) = x(2:nx+1)
      CT(2)%state = exp(-((x2 - x_c)/s_c)**2) !/sqrt(weight)

      ! RK lyap & riccati
      Qc   = eye(rk_c)
      Rinv = eye(rk_b)

      allocate(mat(N, rk_b), matW(N, rk_b))
      call get_state(mat(:,1:rk_b), B(1:rk_b), 'initialize_parameters')
      matW = matmul(mat, weight_mat(:rk_b,:rk_b)) ! incorporate weights
      BBTW = matmul(mat, transpose(matW))
      BRinvBTW_mat  = matmul(mat, matmul(Rinv, transpose(matW)))
      deallocate(mat, matW)

      allocate(mat(N, rk_c), matW(N, rk_c))
      call get_state(mat(:,1:rk_c), CT(1:rk_c), 'initialize_parameters')
      matW = matmul(mat, weight_mat(:rk_c,:rk_c)) ! incorporate weights
      CTCW = matmul(mat, transpose(matW))
      CTQcCW_mat =  matmul(mat, matmul(Qc, transpose(matW)))
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

      return
   end subroutine initialize_parameters

   !--------------------------------------------------------------------
   !-----     UTILITIES FOR STATE_VECTOR AND STATE MATRIX TYPES    -----
   !--------------------------------------------------------------------

   subroutine get_state(mat_out, state_in, procedure)
      !! Utility function to transfer data from a state vector to a real array
      real(wp),                   intent(out) :: mat_out(:,:)
      class(abstract_vector_rdp), intent(in)  :: state_in(:)
      character(len=*),           intent(in)  :: procedure
      ! internal variables
      integer :: k, kdim
      mat_out = 0.0_wp
      select type (state_in)
      type is (state_vector)
         kdim = size(state_in)
         call assert_shape(mat_out, [ N, kdim ], 'mat_out', this_module, 'get_state:'//trim(procedure))
         do k = 1, kdim
            mat_out(:,k) = state_in(k)%state
         end do
      type is (state_matrix)
         call assert_shape(mat_out, [ N, N ], 'mat_out', this_module, 'get_state:'//trim(procedure))
         mat_out = reshape(state_in(1)%state, [ N, N ])
      end select
      return
   end subroutine get_state

   subroutine set_state(state_out, mat_in, procedure)
      !! Utility function to transfer data from a real array to a state vector
      class(abstract_vector_rdp), intent(out) :: state_out(:)
      real(wp),                   intent(in)  :: mat_in(:,:)
      character(len=*),           intent(in)  :: procedure
      ! internal variables
      integer       :: k, kdim
      select type (state_out)
      type is (state_vector)
         kdim = size(state_out)
         call assert_shape(mat_in, [ N, kdim ], 'mat_in', this_module, 'set_state:'//trim(procedure))
         call zero_basis(state_out)
         do k = 1, kdim
            state_out(k)%state = mat_in(:,k)
         end do
      type is (state_matrix)
         call assert_shape(mat_in, [ N, N ], 'mat_in', this_module, 'set_state:'//trim(procedure))
         call zero_basis(state_out)
         state_out(1)%state = reshape(mat_in, shape(state_out(1)%state))
      end select
      return
   end subroutine set_state

   subroutine init_rand(state, ifnorm)
      !! Utility function to initialize a state vector with random data
      class(abstract_vector_rdp), intent(inout)  :: state(:)
      logical, optional,          intent(in)     :: ifnorm
      ! internal variables
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
      end select
      return
   end subroutine init_rand

   subroutine reconstruct_solution_X(X, LR_X)
      real(wp),          intent(out) :: X(:,:)
      type(LR_state),    intent(in)  :: LR_X
      
      ! internals
      real(wp) :: Uwrk(N, LR_X%rk)

      call assert_shape(X, [ N, N ], 'X', this_module, 'reconstruct_solution_X')
      call get_state(Uwrk, LR_X%U(1:LR_X%rk), 'reconstruct_solution_X')
      X = matmul(matmul(Uwrk, matmul(LR_X%S(1:LR_X%rk,1:LR_X%rk), transpose(Uwrk))), weight_mat)

      return
   end subroutine reconstruct_solution_X

   subroutine reconstruct_solution_US(X, U, S)
      real(wp),           intent(out) :: X(:,:)
      type(state_vector), intent(in)  :: U(:)
      real(wp),           intent(in)  :: S(:,:)
      
      ! internals
      integer  :: rk
      real(wp) :: Uwrk(N, size(U))

      rk = size(U)
      call assert_shape(X, [ N, N ], 'X', this_module, 'reconstruct_solution_US')
      call assert_shape(S, [ rk, rk ], 'S', this_module, 'reconstruct_solution_US')
      call get_state(Uwrk, U, 'reconstruct_solution_US')
      X = matmul(matmul(Uwrk, matmul(S, transpose(Uwrk))), weight_mat)

      return
   end subroutine reconstruct_solution_US

   !------------------------------------
   !-----     INITIAL CONDIIONS    -----
   !------------------------------------

   subroutine generate_random_initial_condition(U, S, rk)
      class(state_vector),   intent(out) :: U(:)
      real(wp),              intent(out) :: S(:,:)
      integer,               intent(in)  :: rk
      ! internals
      class(state_vector),   allocatable :: Utmp(:)
      ! SVD
      real(wp)                           :: U_svd(rk,rk)
      real(wp)                           :: S_svd(rk)
      real(wp)                           :: V_svd(rk,rk)
      integer                            :: i, info
      character(len=128) :: msg

      if (size(U) < rk) then
         write(msg,'(A,I0)') 'Input krylov basis size incompatible with requested rank ', rk
         call stop_error(msg, module=this_module, procedure='generate_random_initial_condition')
         STOP 1
      else
         call zero_basis(U)
         do i = 1,rk
            call U(i)%rand(.false.)
         end do
      end if
      call assert_shape(S, [ rk,rk ], 'S', this_module, 'generate_random_initial_condition')
      S = 0.0_wp
      
      ! perform QR
      allocate(Utmp(rk), source=U(:rk))
      call qr(Utmp, S, info)
      call check_info(info, 'qr', module=this_module, procedure='generate_random_initial_condition')
      ! perform SVD
      call svd(S(:rk,:rk), S_svd, U_svd, V_svd)
      S(:rk,:rk) = diag(S_svd)
      block
         class(abstract_vector_rdp), allocatable :: Xwrk(:)
         call linear_combination(Xwrk, Utmp, U_svd)
         call copy(U, Xwrk)
      end block
      write(msg,'(A,I0,A,I0,A)') 'size(U) = [ ', size(U),' ]: filling the first ', rk, ' columns with noise.'
      call logger%log_information(msg, module=this_module, procedure='generate_random_initial_condition')
      return
   end subroutine

   !-------------------------
   !-----      MISC     -----
   !-------------------------

   function CALE(X, adjoint) result(res)
      
      ! solution
      real(wp)          :: X(N,N)
      ! adjoint
      logical, optional :: adjoint
      ! residual
      real(wp)          :: res(N,N)

      ! internals
      real(wp), dimension(N,N) :: AX, XAH

      AX = 0.0_wp; XAH = 0.0_wp
      call GL_mat(AX,  X,             adjoint = adjoint, transpose = .false.)
      call GL_mat(XAH, transpose(X),  adjoint = adjoint, transpose = .true. )

      ! construct Lyapunov equation
      if (adjoint) then
         res = AX + XAH + CTCW
      else
         res = AX + XAH + BBTW
      end if

   end function CALE

   function CARE(X, CTQcCW, BRinvBTW, adjoint) result(res)
      ! solution
      real(wp)          :: X(N,N)
      ! inhomogeneity
      real(wp)          :: CTQcCW(N,N)
      ! inhomogeneity
      real(wp)          :: BRinvBTW(N,N)
      ! adjoint
      logical, optional :: adjoint
      ! residual
      real(wp)          :: res(N,N)
      
      ! internals
      real(wp), dimension(N,N) :: AX, XAH

      AX = 0.0_wp; XAH = 0.0_wp
      call GL_mat(AX,  X,             adjoint = adjoint, transpose = .false.)
      call GL_mat(XAH, transpose(X),  adjoint = adjoint, transpose = .true. )

      ! construct Lyapunov equation
      res = AX + XAH + CTCW + matmul(X, matmul(BRinvBTW, X))

   end function CARE
end module Ginzburg_Landau_Utils