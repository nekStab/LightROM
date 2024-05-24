module Ginzburg_Landau_RK_Lyapunov
   ! LightKrylov for linear algebra.
   use LightKrylov
   use LightKrylov, only : wp => dp
   use LightKrylov_utils, only : assert_shape
   ! LightROM
   use LightROM_AbstractLTIsystems
   use LightROM_Utils ! zero_basis for now
   ! Ginzburg Landau
   use Ginzburg_Landau_Base
   use Ginzburg_Landau_Operators
   ! RKLIB module for time integration.
   use rklib_module
   ! Standard Library.
   use stdlib_optval, only : optval
   implicit none

   private
   public :: CALE, CARE, GL_mat
   public :: get_state_mat, set_state_mat, init_rand

   !-------------------------------------------
   !-----     LIGHTKRYLOV VECTOR TYPE     -----
   !-------------------------------------------

   type, extends(abstract_vector_rdp), public :: state_matrix
      real(wp) :: state(N**2) = 0.0_wp
   contains
      private
      procedure, pass(self), public :: zero => matrix_zero
      procedure, pass(self), public :: dot => matrix_dot
      procedure, pass(self), public :: scal => matrix_scal
      procedure, pass(self), public :: axpby => matrix_axpby
      procedure, pass(self), public :: rand => matrix_rand
   end type state_matrix

   !-------------------------------
   !-----     RK LYAPUNOV     -----
   !-------------------------------

   type, extends(abstract_linop_rdp), public :: RK_lyapunov
      real(wp) :: tau ! Integration time.
   contains
      private
      procedure, pass(self), public :: matvec => direct_solver_lyap
      procedure, pass(self), public :: rmatvec => adjoint_solver_lyap
   end type RK_lyapunov

contains

   subroutine CALE(res_flat, x_flat, Q_flat, adjoint)
      ! residual
      real(wp),                 intent(out) :: res_flat(:)
      ! solution
      real(wp),                 intent(in)  :: x_flat(:)
      ! inhomogeneity
      real(wp),                 intent(in)  :: Q_flat(:)
      !> Adjoint
      logical, optional :: adjoint
      logical           :: adj

      ! internals
      real(wp),   dimension(N**2) :: x_tmp, AX_flat, XAH_flat

      !> Deal with optional argument
      adj  = optval(adjoint,.false.)

      res_flat = 0.0_wp; AX_flat = 0.0_wp; XAH_flat = 0.0_wp; x_tmp = 0.0_wp
      call GL_mat( AX_flat, x_flat, adjoint = adj, transpose = .false.)
      x_tmp    = reshape(transpose(reshape(x_flat,   (/ N,N /))), shape(x_flat))
      call GL_mat(XAH_flat, x_tmp,  adjoint = adj, transpose = .true. )
      ! construct Lyapunov equation
      res_flat = AX_flat + XAH_flat + Q_flat

   end subroutine CALE

   subroutine CARE(res_flat, x_flat, CTQcC_flat, BRinvBT_mat, adjoint)
      ! residual
      real(wp),                 intent(out) :: res_flat(:)
      ! solution
      real(wp),                 intent(in)  :: x_flat(:)
      ! inhomogeneity
      real(wp),                 intent(in)  :: CTQcC_flat(:)
      ! inhomogeneity
      real(wp),                 intent(in)  :: BRinvBT_mat(:,:)
      !> Adjoint
      logical, optional :: adjoint
      logical           :: adj

      ! internals
      real(wp),   dimension(N**2) :: x_tmp, AX_flat, XAH_flat, NL_flat
      real(wp),   dimension(N,N)  :: x_mat

      !> Deal with optional argument
      adj  = optval(adjoint,.false.)

      res_flat = 0.0_wp; AX_flat = 0.0_wp; XAH_flat = 0.0_wp; x_tmp = 0.0_wp
      call GL_mat( AX_flat, x_flat, adjoint = adj, transpose = .false.)
      x_mat = reshape(x_flat, (/ N,N /))
      x_tmp = reshape(transpose(x_mat), shape(x_flat))
      call GL_mat(XAH_flat, x_tmp,  adjoint = adj, transpose = .true. )
      NL_flat = reshape(matmul(x_mat, matmul(BRinvBTW_mat, x_mat)), shape(NL_flat))
      ! construct Lyapunov equation
      res_flat = AX_flat + XAH_flat + CTQcC_flat + NL_flat

   end subroutine CARE

   !-----     TYPE-BOUND PROCEDURE FOR MATRICES     -----

   subroutine matrix_zero(self)
      class(state_matrix), intent(inout) :: self
      self%state = 0.0_wp
      return
   end subroutine matrix_zero

   real(wp) function matrix_dot(self, vec) result(alpha)
      class(state_matrix),        intent(in) :: self
      class(abstract_vector_rdp), intent(in) :: vec
      select type(vec)
      type is(state_matrix)
         alpha = dot_product(self%state, weight_mat*vec%state)
      end select
      return
   end function matrix_dot

   subroutine matrix_scal(self, alpha)
      class(state_matrix), intent(inout) :: self
      real(wp),            intent(in)    :: alpha
      self%state = self%state * alpha
      return
   end subroutine matrix_scal  

   subroutine matrix_axpby(self, alpha, vec, beta)
      class(state_matrix),        intent(inout) :: self
      class(abstract_vector_rdp), intent(in)    :: vec
      real(wp)         , intent(in)    :: alpha, beta
      select type(vec)
      type is(state_matrix)
         self%state = alpha*self%state + beta*vec%state
      end select
      return
   end subroutine matrix_axpby

   subroutine matrix_rand(self, ifnorm)
      class(state_matrix), intent(inout) :: self
      logical, optional,   intent(in)    :: ifnorm
      ! internals
      logical :: normalize
      real(wp) :: alpha
      normalize = optval(ifnorm, .true.)
      call random_number(self%state)
      if (normalize) then
         alpha = self%norm()
         call self%scal(1.0/alpha)
      endif
      return
   end subroutine matrix_rand

   subroutine GL_mat(flat_mat_out, flat_mat_in, adjoint, transpose)
      
      !> State vector.
      real(wp), dimension(:), intent(in)  :: flat_mat_in
      !> Time-derivative.
      real(wp), dimension(:), intent(out) :: flat_mat_out
      !> Adjoint
      logical, optional :: adjoint
      logical           :: adj
      logical, optional :: transpose
      logical           :: trans
      
      !> Internal variables.
      integer :: j
      real(wp), dimension(N,N) :: mat, dmat
      
      !> Deal with optional argument
      adj   = optval(adjoint,.false.)
      trans = optval(transpose,.false.)

      !> Sets the internal variables.
      mat  = reshape(flat_mat_in(1:N**2),(/N, N/))
      dmat = 0.0_wp
      
      if (adj) then
         if (trans) then
            do j = 1,N
               call adjoint_GL(mat(:,j), dmat(j,:))
            end do
         else
            do j = 1,N
               call adjoint_GL(mat(:,j), dmat(:,j))
            end do
         end if
      else
         if (trans) then
            do j = 1,N
               call direct_GL(mat(:,j), dmat(j,:))
            end do
         else
            do j = 1,N
               call direct_GL(mat(:,j), dmat(:,j))
            end do
         end if
      endif

      !> Reshape for output
      flat_mat_out = reshape(dmat, shape(flat_mat_in))
      
      return
   end subroutine GL_mat

   !--------------------------------------
   !-----     WRAPPERS FOR RKLIB     -----
   !--------------------------------------

   subroutine rhs_lyap(me, t, x_flat, f_flat)
      ! Time-integrator.
      class(rk_class), intent(inout)             :: me
      ! Current time.
      real(wp),        intent(in)                :: t
      ! State vector.
      real(wp),        dimension(:), intent(in)  :: x_flat
      ! Time-derivative.
      real(wp),        dimension(:), intent(out) :: f_flat

      ! internals
      real(wp),        dimension(N**2) :: x_tmp, AX_flat, XAH_flat

      f_flat = 0.0_wp; AX_flat = 0.0_wp; XAH_flat = 0.0_wp; x_tmp = 0.0_wp
      ! A @ X
      call GL_mat( AX_flat, x_flat, adjoint = .false., transpose = .false.)
      ! build X.T
      x_tmp    = reshape(transpose(reshape(x_flat,   (/ N,N /))), shape(x_flat))
      ! build ( A @ X.T ).T = X @ A.T
      call GL_mat(XAH_flat, x_tmp,  adjoint = .false., transpose = .true.)
      ! construct Lyapunov equation
      f_flat = AX_flat + XAH_flat + BBTW_flat

      return
   end subroutine rhs_lyap

   subroutine adjoint_rhs_lyap(me, t, x_flat, f_flat)
      ! Time-integrator.
      class(rk_class), intent(inout)             :: me
      ! Current time.
      real(wp),        intent(in)                :: t
      ! State vector.
      real(wp),        dimension(:), intent(in)  :: x_flat
      ! Time-derivative.
      real(wp),        dimension(:), intent(out) :: f_flat

      ! internals
      real(wp),        dimension(N**2) :: x_tmp, AHX_flat, XA_flat

      f_flat = 0.0_wp; AHX_flat = 0.0_wp; XA_flat = 0.0_wp; x_tmp = 0.0_wp
      ! A.T @ X
      call GL_mat(AHX_flat, x_flat, adjoint = .true., transpose = .false.)
      ! build X.T
      x_tmp   = reshape(transpose(reshape(x_flat,  (/ N,N /))), shape(x_flat))
      ! build ( A.T @ X.T ).T = X @ A
      call GL_mat( XA_flat, x_tmp,  adjoint = .true., transpose = .true.)
      ! construct Lyapunov equation
      f_flat = AHX_flat + XA_flat + CTCW_flat

      return
   end subroutine adjoint_rhs_lyap

   !------------------------------------------------------------------------
   !-----     TYPE-BOUND PROCEDURES FOR THE EXPONENTIAL PROPAGATOR     -----
   !------------------------------------------------------------------------

   subroutine direct_solver_lyap(self, vec_in, vec_out)
      ! Linear Operator.
      class(rk_lyapunov),          intent(in)  :: self
      ! Input vector.
      class(abstract_vector_rdp),  intent(in)  :: vec_in
      ! Output vector.
      class(abstract_vector_rdp),  intent(out) :: vec_out

      ! Time-integrator.
      type(rks54_class) :: prop
      real(kind=wp)     :: dt = 1.0_wp

      select type(vec_in)
      type is(state_matrix)
         select type(vec_out)
         type is(state_matrix)
            ! Initialize propagator.
            call prop%initialize(n=N**2, f=rhs_lyap)
            ! Integrate forward in time.
            call prop%integrate(0.0_wp, vec_in%state, dt, self%tau, vec_out%state)
         end select
      end select
      return
   end subroutine direct_solver_lyap

   subroutine adjoint_solver_lyap(self, vec_in, vec_out)
      ! Linear Operator.
      class(rk_lyapunov),          intent(in)  :: self
      ! Input vector.
      class(abstract_vector_rdp),  intent(in)  :: vec_in
      ! Output vector.
      class(abstract_vector_rdp),  intent(out) :: vec_out

      ! Time-integrator.
      type(rks54_class) :: prop
      real(kind=wp)     :: dt = 1.0_wp

      select type(vec_in)
      type is(state_matrix)
         select type(vec_out)
         type is(state_matrix)
            ! Initialize propagator.
            call prop%initialize(n=N**2, f=adjoint_rhs_lyap)
            ! Integrate forward in time.
            call prop%integrate(0.0_wp, vec_in%state, dt, self%tau, vec_out%state)
         end select
      end select
      return
   end subroutine adjoint_solver_lyap

   !-----------------------------------------------
   !-----     UTILITIES FOR STATE_MATRICES    -----
   !-----------------------------------------------

   subroutine get_state_mat(mat_out, state_in)
      !! Utility function to transfer data from a state vector to a real array
      real(wp),                   intent(out) :: mat_out(:,:)
      class(abstract_vector_rdp), intent(in)  :: state_in(:)
      ! internal variables
      mat_out = 0.0_wp
      select type (state_in)
      type is (state_matrix)
         call assert_shape(mat_out, (/ N, N /), 'get_state -> state_matrix', 'mat_out')
         mat_out = reshape(state_in(1)%state, (/ N, N /))
      end select
      return
   end subroutine get_state_mat

   subroutine set_state_mat(state_out, mat_in)
      !! Utility function to transfer data from a real array to a state vector
      class(abstract_vector_rdp), intent(out) :: state_out(:)
      real(wp),                   intent(in)  :: mat_in(:,:)
      ! internal variables
      select type (state_out)
      type is (state_matrix)
         call assert_shape(mat_in, (/ N, N /), 'set_state -> state_matrix', 'mat_in')
         call zero_basis(state_out)
         state_out(1)%state = reshape(mat_in, shape(state_out(1)%state))
      end select
      return
   end subroutine set_state_mat

   subroutine init_rand_mat(state, ifnorm)
      !! Utility function to initialize a state vector with random data
      class(abstract_vector_rdp), intent(inout)  :: state(:)
      logical, optional,          intent(in)     :: ifnorm
      ! internal variables
      integer :: k, kdim
      logical :: normalize
      normalize = optval(ifnorm,.true.)
      select type (state)
      type is (state_matrix)
         kdim = size(state)
         do k = 1, kdim
            call state(k)%rand(ifnorm = normalize)
         end do
      end select
      return
   end subroutine init_rand_mat

end module Ginzburg_Landau_RK_Lyapunov