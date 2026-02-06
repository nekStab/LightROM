module Ginzburg_Landau_RKlib
   ! Standard Library.
   use stdlib_stats_distribution_normal, only: normal => rvs_normal
   use stdlib_optval, only : optval
   ! RKLIB module for time integration.
   use rklib_module
   ! LightKrylov for linear algebra.
   use LightKrylov
   ! Ginzburg Landau
   use Ginzburg_Landau_Base
   use Ginzburg_Landau_Operators
   implicit none

   private :: this_module
   character(len=*), parameter :: this_module = 'Ginzburg_Landau_RKlib'
   
   public  :: GL_mat

   !-------------------------------------------
   !-----     LIGHTKRYLOV VECTOR TYPE     -----
   !-------------------------------------------

   type, extends(abstract_vector_rdp), public :: state_matrix
      real(dp) :: state(N**2) = 0.0_dp
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
      real(dp) :: tau ! Integration time.
   contains
      private
      procedure, pass(self), public :: matvec => direct_solver_lyap
      procedure, pass(self), public :: rmatvec => adjoint_solver_lyap
   end type RK_lyapunov

   !------------------------------
   !-----     RK RICCATI     -----
   !------------------------------

   type, extends(abstract_linop_rdp), public :: RK_riccati
      real(dp) :: tau ! Integration time.
   contains
      private
      procedure, pass(self), public :: matvec => direct_solver_ricc
      procedure, pass(self), public :: rmatvec => adjoint_solver_ricc
   end type RK_riccati

contains

   !-----     TYPE-BOUND PROCEDURE FOR MATRICES     -----

   subroutine matrix_zero(self)
      class(state_matrix), intent(inout) :: self
      self%state = 0.0_dp
   end subroutine matrix_zero

   real(dp) function matrix_dot(self, vec) result(alpha)
      class(state_matrix),        intent(in) :: self
      class(abstract_vector_rdp), intent(in) :: vec
      select type(vec)
      type is(state_matrix)
         alpha = dot_product(self%state, weight_flat*vec%state)
      class default
         call stop_error('vec must be a state_matrix', this_module, 'matrix_dot')
      end select
   end function matrix_dot

   integer function matrix_get_size(self) result(N)
     class(state_matrix), intent(in) :: self
     N = N**2
   end function matrix_get_size

   subroutine matrix_scal(self, alpha)
      class(state_matrix), intent(inout) :: self
      real(dp),            intent(in)    :: alpha
      self%state = self%state * alpha
   end subroutine matrix_scal  

   subroutine matrix_axpby(alpha, vec, beta, self)
      class(state_matrix),        intent(inout) :: self
      class(abstract_vector_rdp), intent(in)    :: vec
      real(dp)         , intent(in)    :: alpha, beta
      select type(vec)
      type is(state_matrix)
         self%state = beta*self%state + alpha*vec%state
      class default
         call stop_error('vec must be a state_matrix', this_module, 'matrix_axpby')
      end select
   end subroutine matrix_axpby

   subroutine matrix_rand(self, ifnorm)
      class(state_matrix), intent(inout) :: self
      logical, optional,   intent(in)    :: ifnorm
      ! internals
      logical :: normalize
      real(dp) :: alpha
      real(dp), dimension(N**2) :: mean, std
      normalize = optval(ifnorm, .true.)
      mean = 0.0_dp
      std  = 1.0_dp
      self%state = normal(mean,std)
      if (normalize) then
         alpha = self%norm()
         call self%scal(1.0/alpha)
      endif
   end subroutine matrix_rand

   subroutine GL_mat(mat_out, mat_in, adjoint, transpose)
      
      ! State vector.
      real(dp), dimension(:,:), intent(in)  :: mat_in
      ! Time-derivative.
      real(dp), dimension(:,:), intent(out) :: mat_out
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
      
      mat_out = 0.0_dp
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

   end subroutine GL_mat

   !--------------------------------------
   !-----     WRAPPERS FOR RKLIB     -----
   !--------------------------------------

   subroutine rhs_lyap(me, t, x_flat, f_flat)
      ! Time-integrator.
      class(rk_class), intent(inout)             :: me
      ! Current time.
      real(dp),        intent(in)                :: t
      ! State vector.
      real(dp),        dimension(:), intent(in)  :: x_flat
      ! Time-derivative.
      real(dp),        dimension(:), intent(out) :: f_flat

      ! internals
      real(dp),        dimension(N,N) :: X, AX, XAT

      f_flat = 0.0_dp
      AX = 0.0_dp; XAT = 0.0_dp;
      X = reshape(x_flat, [N,N])
      ! A @ X
      call GL_mat(AX, X,              adjoint = .false., transpose = .false.)
      ! build ( A @ X.T ).T = X @ A.T
      call GL_mat(XAT, transpose(X),  adjoint = .false., transpose = .true.)
      ! construct Lyapunov equation
      f_flat = reshape(AX + XAT + BBTW, [ N**2 ])

   end subroutine rhs_lyap

   subroutine adjoint_rhs_lyap(me, t, x_flat, f_flat)
      ! Time-integrator.
      class(rk_class), intent(inout)             :: me
      ! Current time.
      real(dp),        intent(in)                :: t
      ! State vector.
      real(dp),        dimension(:), intent(in)  :: x_flat
      ! Time-derivative.
      real(dp),        dimension(:), intent(out) :: f_flat

      ! internals
      real(dp),        dimension(N,N) :: X, ATX, XA

      f_flat = 0.0_dp
      ATX = 0.0_dp; XA = 0.0_dp;
      X = reshape(x_flat, [N,N])
      ! A.T @ X
      call GL_mat(ATX, X,            adjoint = .true., transpose = .false.)
      ! build ( A.T @ X.T ).T = X @ A
      call GL_mat( XA, transpose(X), adjoint = .true., transpose = .true.)
      ! construct Lyapunov equation
      f_flat = reshape(ATX + XA + CTCW, [ N**2 ])

   end subroutine adjoint_rhs_lyap

   subroutine rhs_ricc(me, t, x_flat, f_flat)
      ! Time-integrator.
      class(rk_class), intent(inout)             :: me
      ! Current time.
      real(dp),        intent(in)                :: t
      ! State vector.
      real(dp),        dimension(:), intent(in)  :: x_flat
      ! Time-derivative.
      real(dp),        dimension(:), intent(out) :: f_flat

      ! internals
      real(dp),        dimension(N,N)  :: X, ATX, XA

      f_flat = 0.0_dp
      ATX = 0.0_dp; XA = 0.0_dp;
      X = reshape(x_flat, [N,N])
      ! A.T @ X
      call GL_mat(ATX, X,             adjoint = .true., transpose = .false.)
      ! build ( A.T @ X.T ).T = X @ A
      call GL_mat( XA, transpose(X),  adjoint = .true., transpose = .true.)
      ! construct Lyapunov equation
      f_flat = reshape(ATX + XA + CTQcCW - matmul(X, matmul(BRinvBTW, X)), [ N**2 ])

   end subroutine rhs_ricc

   subroutine adjoint_rhs_ricc(me, t, x_flat, f_flat)
      ! Time-integrator.
      class(rk_class), intent(inout)             :: me
      ! Current time.
      real(dp),        intent(in)                :: t
      ! State vector.
      real(dp),        dimension(:), intent(in)  :: x_flat
      ! Time-derivative.
      real(dp),        dimension(:), intent(out) :: f_flat

      ! internals
      real(dp),        dimension(N,N)  :: X, AX, XAT

      f_flat = 0.0_dp
      AX = 0.0_dp; XAT = 0.0_dp;
      X = reshape(x_flat, [N,N])
      ! A @ X
      call GL_mat(AX, X,              adjoint = .false., transpose = .false.)
      ! build ( A @ X.T ).T = X @ A.T
      call GL_mat(XAT, transpose(X),  adjoint = .false., transpose = .true.)
      ! construct Lyapunov equation
      f_flat = reshape(AX + XAT + BQeBTW - matmul(X, matmul(CTVinvCW, X)), [ N**2 ])

   end subroutine adjoint_rhs_ricc

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
      character(len=*), parameter :: this_procedure = 'direct_solver_lyap'
      type(rks54_class) :: prop
      real(kind=dp)     :: dt = 1.0_dp

      select type(vec_in)
      type is(state_matrix)
         select type(vec_out)
         type is(state_matrix)
            ! Initialize propagator.
            call prop%initialize(n=N**2, f=rhs_lyap)
            ! Integrate forward in time.
            call prop%integrate(0.0_dp, vec_in%state, dt, self%tau, vec_out%state)
         class default
            call type_error('vec_out', 'state_vector', 'OUT', this_module, this_procedure)
         end select
      class default
         call type_error('vec_in', 'state_vector', 'IN', this_module, this_procedure)
      end select
   end subroutine direct_solver_lyap

   subroutine adjoint_solver_lyap(self, vec_in, vec_out)
      ! Linear Operator.
      class(rk_lyapunov),          intent(inout)  :: self
      ! Input vector.
      class(abstract_vector_rdp),  intent(in)  :: vec_in
      ! Output vector.
      class(abstract_vector_rdp),  intent(out) :: vec_out

      ! Time-integrator.
      character(len=*), parameter :: this_procedure = 'adjoint_solver_lyap'
      type(rks54_class) :: prop
      real(kind=dp)     :: dt = 1.0_dp

      select type(vec_in)
      type is(state_matrix)
         select type(vec_out)
         type is(state_matrix)
            ! Initialize propagator.
            call prop%initialize(n=N**2, f=adjoint_rhs_lyap)
            ! Integrate forward in time.
            call prop%integrate(0.0_dp, vec_in%state, dt, self%tau, vec_out%state)
         class default
            call type_error('vec_out', 'state_vector', 'OUT', this_module, this_procedure)
         end select
      class default
         call type_error('vec_in', 'state_vector', 'IN', this_module, this_procedure)
      end select
   end subroutine adjoint_solver_lyap

   subroutine direct_solver_ricc(self, vec_in, vec_out)
      ! Linear Operator.
      class(rk_riccati) ,          intent(inout)  :: self
      ! Input vector.
      class(abstract_vector_rdp),  intent(in)  :: vec_in
      ! Output vector.
      class(abstract_vector_rdp),  intent(out) :: vec_out

      ! Time-integrator.
      character(len=*), parameter :: this_procedure = 'direct_solver_ricc'
      type(rks54_class) :: prop
      real(kind=dp)     :: dt = 1.0_dp

      select type(vec_in)
      type is(state_matrix)
         select type(vec_out)
         type is(state_matrix)
            ! Initialize propagator.
            call prop%initialize(n=N**2, f=rhs_ricc)
            ! Integrate forward in time.
            call prop%integrate(0.0_dp, vec_in%state, dt, self%tau, vec_out%state)
         class default
            call type_error('vec_out', 'state_vector', 'OUT', this_module, this_procedure)
         end select
      class default
         call type_error('vec_in', 'state_vector', 'IN', this_module, this_procedure)
      end select
      return
   end subroutine direct_solver_ricc

   subroutine adjoint_solver_ricc(self, vec_in, vec_out)
      ! Linear Operator.
      class(rk_riccati) ,          intent(inout)  :: self
      ! Input vector.
      class(abstract_vector_rdp),  intent(in)  :: vec_in
      ! Output vector.
      class(abstract_vector_rdp),  intent(out) :: vec_out

      ! Time-integrator.
      character(len=*), parameter :: this_procedure = 'adjoint_solver_ricc'
      type(rks54_class) :: prop
      real(kind=dp)     :: dt = 1.0_dp

      select type(vec_in)
      type is(state_matrix)
         select type(vec_out)
         type is(state_matrix)
            ! Initialize propagator.
            call prop%initialize(n=N**2, f=adjoint_rhs_ricc)
            ! Integrate forward in time.
            call prop%integrate(0.0_dp, vec_in%state, dt, self%tau, vec_out%state)
         class default
            call type_error('vec_out', 'state_vector', 'OUT', this_module, this_procedure)
         end select
      class default
         call type_error('vec_in', 'state_vector', 'IN', this_module, this_procedure)
      end select
      return
   end subroutine adjoint_solver_ricc

end module Ginzburg_Landau_RKlib