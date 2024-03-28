module cartpole_riccati
   use LightKrylov
   use rklib_module
   use cartpole
   use stdlib_linalg, only: eye
   use stdlib_optval, only: optval

   implicit none

   !-----------------------------------------------------------
   !-----     LIGHTKRYLOV VECTOR TYPE FOR MATRIX STATE    -----
   !-----------------------------------------------------------

   type, extends(abstract_vector), public :: state_matrix
      real(kind=wp) :: state(n**2) = 0.0_wp
   contains
      private
      procedure, pass(self), public :: zero
      procedure, pass(self), public :: dot
      procedure, pass(self), public :: scal
      procedure, pass(self), public :: axpby
      procedure, pass(self), public :: rand
   end type state_matrix

   type, extends(abstract_linop), public :: rklib_riccati_mat
      real(kind=wp) :: tau ! Integration time.
   contains
      private
      procedure, pass(self), public :: matvec  => direct_solver_mat
      procedure, pass(self), public :: rmatvec => direct_solver_mat                ! dummy, not used
   end type rklib_riccati_mat

contains

   !---------------------------------------------------------
   !-----     TYPE-BOUND PROCEDURE FOR MATRIX STATE     -----
   !---------------------------------------------------------

   subroutine zero(self)
      class(state_matrix), intent(inout) :: self
      self%state = 0.0_wp
      return
   end subroutine zero

   real(kind=wp) function dot(self, vec) result(alpha)
      class(state_matrix)   , intent(in) :: self
      class(abstract_vector), intent(in) :: vec
      select type(vec)
      type is (state_matrix)
         alpha = dot_product(self%state, vec%state)
      end select
      return
   end function dot

   subroutine scal(self, alpha)
      class(state_matrix), intent(inout) :: self
      real(kind=wp)      , intent(in)    :: alpha
      self%state = self%state * alpha
      return
   end subroutine scal

   subroutine axpby(self, alpha, vec, beta)
      class(state_matrix)   , intent(inout) :: self
      class(abstract_vector), intent(in)    :: vec
      real(kind=wp)         , intent(in)    :: alpha, beta
      select type(vec)
      type is (state_matrix)
         self%state = alpha*self%state + beta*vec%state
      end select
      return
   end subroutine axpby

   subroutine rand(self, ifnorm)
      class(state_matrix), intent(inout) :: self
      logical, optional,   intent(in)    :: ifnorm
      ! internals
      logical :: normalize
      real(kind=wp) :: alpha
      normalize = optval(ifnorm,.true.)
      call random_number(self%state)
      if (normalize) then
         alpha = self%norm()
         call self%scal(1.0/alpha)
      endif
      return
   end subroutine rand

   !--------------------------------------------
   !-----    Riccati equation for RKlib    -----
   !--------------------------------------------

   subroutine rhs_riccati(me, t, x, f)
      !> Time-integrator.
      class(rk_class), intent(inout)             :: me
      !> Current time.
      real(kind=wp)  , intent(in)                :: t
      !> State vector.
      real(kind=wp)  , dimension(:), intent(in)  :: x
      !> Time-derivative.
      real(kind=wp)  , dimension(:), intent(out) :: f

      !> Internal variables.
      integer :: i, j, k
      real(kind=wp), dimension(N,N) :: xm, vm

      xm = reshape(x, shape(xm))

      !> Sets the internal variables.
      vm  = 0.0_wp

      vm = matmul(A, xm) + matmul(xm, transpose(A)) + CTQC - matmul(xm, matmul(BRinvBT, xm))

      !> Combine the parts and copy result to the output array.
      f(1:N**2) = reshape(vm, shape(f))
      
      return
   end subroutine rhs_riccati

   !------------------------------------------------------------------------
   !-----     TYPE-BOUND PROCEDURES FOR THE EXPONENTIAL PROPAGATOR     -----
   !------------------------------------------------------------------------

   subroutine direct_solver_mat(self, vec_in, vec_out)
      !> Linear Operator.
      class(rklib_riccati_mat), intent(in)  :: self
      !> Input vector.
      class(abstract_vector) , intent(in)  :: vec_in
      !> Output vector.
      class(abstract_vector) , intent(out) :: vec_out
      !> Time-integrator.
      type(rks54_class) :: prop
      real(kind=wp)     :: dt = 0.1_wp

      select type(vec_in)
      type is (state_matrix)
         select type(vec_out)
         type is(state_matrix)
            !> Initialize propagator.
            call prop%initialize(n=N**2, f=rhs_riccati)
            !> Integrate forward in time.
            call prop%integrate(0.0_wp, vec_in%state, dt, self%tau, vec_out%state)
         end select
      end select

      return
   end subroutine direct_solver_mat

end module cartpole_riccati