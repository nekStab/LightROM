module Ginzburg_Landau_RK_Lyapunov
   ! Standard Library.
   use stdlib_stats_distribution_normal, only: normal => rvs_normal
   use stdlib_optval, only : optval
   ! RKLIB module for time integration.
   use rklib_module
   ! LightKrylov for linear algebra.
   use LightKrylov
   use LightKrylov, only : wp => dp
   ! Ginzburg Landau
   use Ginzburg_Landau_Base
   use Ginzburg_Landau_Operators
   implicit none

   private :: this_module
   character(len=*), parameter :: this_module = 'Ginzburg_Landau_RK_Lyapunov'
   
   public  :: GL_mat

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
      procedure, pass(self), public :: get_size => matrix_get_size
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
         alpha = dot_product(self%state, weight_flat*vec%state)
      end select
      return
   end function matrix_dot

   integer function matrix_get_size(self) result(N)
     class(state_matrix), intent(in) :: self
     N = N**2
     return
   end function matrix_get_size

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
      real(wp), dimension(N**2) :: mean, std
      normalize = optval(ifnorm, .true.)
      mean = 0.0_wp
      std  = 1.0_wp
      self%state = normal(mean,std)
      if (normalize) then
         alpha = self%norm()
         call self%scal(1.0/alpha)
      endif
      return
   end subroutine matrix_rand

   subroutine GL_mat(mat_out, mat_in, adjoint, transpose)
      
      ! State vector.
      real(wp), dimension(:,:), intent(in)  :: mat_in
      ! Time-derivative.
      real(wp), dimension(:,:), intent(out) :: mat_out
      ! Adjoint
      logical, optional :: adjoint
      logical           :: adj
      logical, optional :: transpose
      logical           :: trans
      
      ! Internal variables.
      integer :: j
      
      ! Deal with optional argument
      adj   = optval(  adjoint,.false.)
      trans = optval(transpose,.false.)
      
      mat_out = 0.0_wp
      if (adj) then
         if (trans) then
            do j = 1,N
               call adjoint_GL(mat_in(:,j), mat_out(j,:))
            end do
         else
            do j = 1,N
               call adjoint_GL(mat_in(:,j), mat_out(:,j))
            end do
         end if
      else
         if (trans) then
            do j = 1,N
               call direct_GL(mat_in(:,j), mat_out(j,:))
            end do
         else
            do j = 1,N
               call direct_GL(mat_in(:,j), mat_out(:,j))
            end do
         end if
      endif
      
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
      real(wp),        dimension(N,N)  :: AX, XAH

      f_flat = 0.0_wp
      AX = 0.0_wp; XAH = 0.0_wp;
      ! A @ X
      call GL_mat( AX, reshape(x_flat, [N,N]),             adjoint = .false., transpose = .false.)
      ! build ( A @ X.T ).T = X @ A.T
      call GL_mat(XAH, transpose(reshape(x_flat, [N,N])),  adjoint = .false., transpose = .true.)
      ! construct Lyapunov equation
      f_flat = reshape(AX + XAH + BBTW, [ N**2 ])

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
      real(wp),        dimension(N,N) ::  AHX, XA

      f_flat = 0.0_wp
      AHX = 0.0_wp; XA = 0.0_wp;
      ! A.T @ X
      call GL_mat(AHX, reshape(x_flat, [N,N]),            adjoint = .true., transpose = .false.)
      ! build ( A.T @ X.T ).T = X @ A
      call GL_mat( XA, transpose(reshape(x_flat, [N,N])), adjoint = .true., transpose = .true.)
      ! construct Lyapunov equation
      f_flat = reshape(AHX + XA +  CTCW, [ N**2 ])

      return
   end subroutine adjoint_rhs_lyap

   !------------------------------------------------------------------------
   !-----     TYPE-BOUND PROCEDURES FOR THE EXPONENTIAL PROPAGATOR     -----
   !------------------------------------------------------------------------

   subroutine direct_solver_lyap(self, vec_in, vec_out)
      ! Linear Operator.
      class(rk_lyapunov),          intent(inout)  :: self
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
      class(rk_lyapunov),          intent(inout)  :: self
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

end module Ginzburg_Landau_RK_Lyapunov