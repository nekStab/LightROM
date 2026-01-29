module Ginzburg_Landau_Operators
   ! Standard Library.
   use stdlib_optval, only : optval
   ! RKLIB module for time integration.
   use rklib_module
   ! LightKrylov for linear algebra.
   use LightKrylov
   use LightKrylov_utils, only : assert_shape
   ! LightROM
   use LightROM_AbstractLTIsystems ! abstract_lti_system
   ! Ginzburg Landau
   use Ginzburg_Landau_Base
   implicit none

   private :: this_module
   character(len=*), parameter :: this_module = 'Ginzburg_Landau_Operators'
   public  :: exptA, direct_GL, adjoint_GL

   !-----------------------------------------------
   !-----     LIGHTKRYLOV LTI SYSTEM TYPE     -----
   !-----------------------------------------------

   type, extends(abstract_lti_system_rdp), public :: lti_system
   contains
      private
      procedure, pass(self), public :: initialize_lti_system
   end type lti_system

   !--------------------------------------
   !-----     LINEAR GL OPERATOR     -----
   !--------------------------------------

   type, extends(abstract_linop_rdp), public :: GL_operator
   contains
      private
      procedure, pass(self), public :: matvec  => direct_matvec_GL
      procedure, pass(self), public :: rmatvec => adjoint_matvec_GL
   end type GL_operator

   !------------------------------------------
   !-----     EXPONENTIAL PROPAGATOR     -----
   !------------------------------------------

   type, extends(abstract_linop_rdp), public :: exponential_prop
      real(dp) :: tau ! Integration time.
   contains
      private
      procedure, pass(self), public :: matvec => direct_solver
      procedure, pass(self), public :: rmatvec => adjoint_solver
   end type exponential_prop

contains

   !========================================================================
   !========================================================================
   !=====                                                              =====
   !=====     PHYSICAL MODEL : LINEARIZED GINZBURG-LANDAU EQUATION     =====
   !=====                                                              =====
   !========================================================================
   !========================================================================

   !---------------------------------------------------------
   !-----      LINEARIZED GINZBURG-LANDAU EQUATIONS     -----
   !---------------------------------------------------------

   subroutine direct_GL(vec_in, vec_out)

      !> State vector.
      real(dp), dimension(:), intent(in)  :: vec_in
      !> Time-derivative.
      real(dp), dimension(:), intent(out) :: vec_out

      !> Internal variables.
      integer                 :: i
      real(dp), dimension(nx) :: u, v, du, dv
      real(dp)                :: d2u, d2v, cu, cv

      u = vec_in(1:nx)     
      v = vec_in(nx+1:2*nx)

      !---------------------------------------------------
      !-----     Linear Ginzburg Landau Equation     -----
      !---------------------------------------------------

      cu = u(2) / (2*dx) ; cv = v(2) / (2*dx)
      du(1) = -(real(nu)*cu - aimag(nu)*cv) ! Convective term.
      dv(1) = -(aimag(nu)*cu + real(nu)*cv) ! Convective term.

      d2u = (u(2) - 2*u(1)) / dx**2 ; d2v = (v(2) - 2*v(1)) / dx**2
      du(1) = du(1) + real(gamma)*d2u - aimag(gamma)*d2v ! Diffusion term.
      dv(1) = dv(1) + aimag(gamma)*d2u + real(gamma)*d2v ! Diffusion term.

      du(1) = du(1) + mu(1)*u(1) ! Non-parallel term.
      dv(1) = dv(1) + mu(1)*v(1) ! Non-parallel term.

      ! Interior nodes.
      do i = 2, nx-1
         ! Convective term.
         cu = (u(i+1) - u(i-1)) / (2*dx)
         cv = (v(i+1) - v(i-1)) / (2*dx)
         du(i) = -(real(nu)*cu - aimag(nu)*cv)
         dv(i) = -(aimag(nu)*cu + real(nu)*cv)

         ! Diffusion term.
         d2u = (u(i+1) - 2*u(i) + u(i-1)) / dx**2
         d2v = (v(i+1) - 2*v(i) + v(i-1)) / dx**2
         du(i) = du(i) + real(gamma)*d2u - aimag(gamma)*d2v
         dv(i) = dv(i) + aimag(gamma)*d2u + real(gamma)*d2v

         ! Non-parallel term.
         du(i) = du(i) + mu(i)*u(i)
         dv(i) = dv(i) + mu(i)*v(i)
      enddo

      ! Right most boundary points.
      cu = -u(nx-1) / (2*dx) ; cv = -v(nx-1) / (2*dx)
      du(nx) = -(real(nu)*cu - aimag(nu)*cv) ! Convective term.
      dv(nx) = -(aimag(nu)*cu + real(nu)*cv) ! Convective term.

      d2u = (-2*u(nx) + u(nx-1)) / dx**2 ; d2v = (-2*v(nx) + v(nx-1)) / dx**2
      du(nx) = du(nx) + real(gamma)*d2u - aimag(gamma)*d2v ! Diffusion term.
      dv(nx) = dv(nx) + aimag(gamma)*d2u + real(gamma)*d2v ! Diffusion term.

      du(nx) = du(nx) + mu(nx)*u(nx) ! Non-parallel term.
      dv(nx) = dv(nx) + mu(nx)*v(nx) ! Non-parallel term.

      vec_out(1:nx)      = du
      vec_out(nx+1:2*nx) = dv

   end subroutine direct_GL

   !-----------------------------------------------------------
   !-----     Adjoint linear Ginzburg-Landau equation     -----
   !-----------------------------------------------------------

   subroutine adjoint_GL(vec_in, vec_out)
      !> State vector.
      real(dp), dimension(:), intent(in)  :: vec_in
      !> Time-derivative.
      real(dp), dimension(:), intent(out) :: vec_out

      ! Internal variables.
      integer :: i
      real(dp), dimension(nx) :: u, du
      real(dp), dimension(nx) :: v, dv
      real(dp)                :: d2u, d2v, cu, cv

      ! Sets the internal variables.
      u = vec_in(1:nx)     
      v = vec_in(nx+1:2*nx)

      !---------------------------------------------------
      !-----     Linear Ginzburg Landau Equation     -----
      !---------------------------------------------------

      ! Left most boundary points.
      cu = u(2) / (2*dx) ; cv = v(2) / (2*dx)
      du(1) = (real(nu)*cu + aimag(nu)*cv) ! Convective term.
      dv(1) = (-aimag(nu)*cu + real(nu)*cv) ! Convective term.

      d2u = (u(2) - 2*u(1)) / dx**2 ; d2v = (v(2) - 2*v(1)) / dx**2
      du(1) = du(1) + real(gamma)*d2u + aimag(gamma)*d2v ! Diffusion term.
      dv(1) = dv(1) - aimag(gamma)*d2u + real(gamma)*d2v ! Diffusion term.

      du(1) = du(1) + mu(1)*u(1) ! Non-parallel term.
      dv(1) = dv(1) + mu(1)*v(1) ! Non-parallel term.

      ! Interior nodes.
      do i = 2, nx-1
         ! Convective term.
         cu = (u(i+1) - u(i-1)) / (2*dx)
         cv = (v(i+1) - v(i-1)) / (2*dx)
         du(i) = (real(nu)*cu + aimag(nu)*cv)
         dv(i) = (-aimag(nu)*cu + real(nu)*cv)

         ! Diffusion term.
         d2u = (u(i+1) - 2*u(i) + u(i-1)) / dx**2
         d2v = (v(i+1) - 2*v(i) + v(i-1)) / dx**2
         du(i) = du(i) + real(gamma)*d2u + aimag(gamma)*d2v
         dv(i) = dv(i) - aimag(gamma)*d2u + real(gamma)*d2v

         ! Non-parallel term.
         du(i) = du(i) + mu(i)*u(i)
         dv(i) = dv(i) + mu(i)*v(i)
      enddo

      ! Right most boundary points.
      cu = -u(nx-1) / (2*dx) ; cv = -v(nx-1) / (2*dx)
      du(nx) = (real(nu)*cu + aimag(nu)*cv) ! Convective term.
      dv(nx) = (-aimag(nu)*cu + real(nu)*cv) ! Convective term.

      d2u = (-2*u(nx) + u(nx-1)) / dx**2 ; d2v = (-2*v(nx) + v(nx-1)) / dx**2
      du(nx) = du(nx) + real(gamma)*d2u + aimag(gamma)*d2v ! Diffusion term.
      dv(nx) = dv(nx) - aimag(gamma)*d2u + real(gamma)*d2v ! Diffusion term.

      du(nx) = du(nx) + mu(nx)*u(nx) ! Non-parallel term.
      dv(nx) = dv(nx) + mu(nx)*v(nx) ! Non-parallel term.

      ! Copy results to the output array.
      vec_out(1:nx)      = du
      vec_out(nx+1:2*nx) = dv

   end subroutine adjoint_GL

   !--------------------------------------
   !-----     WRAPPERS FOR RKLIB     -----
   !--------------------------------------

   subroutine rhs(me, t, x, f)
      ! Time-integrator.
      class(rk_class),               intent(inout) :: me
      ! Current time.
      real(dp),                      intent(in)    :: t
      ! State vector.
      real(dp),        dimension(:), intent(in)    :: x
      ! Time-derivative.
      real(dp),        dimension(:), intent(out)   :: f

      f = 0.0_dp
      call direct_GL(x, f)

   end subroutine rhs

   subroutine adjoint_rhs(me, t, x, f)
      ! Time-integrator.
      class(rk_class),               intent(inout) :: me
      ! Current time.
      real(dp),                      intent(in)    :: t
      ! State vector.
      real(dp),        dimension(:), intent(in)    :: x
      ! Time-derivative.
      real(dp),        dimension(:), intent(out)   :: f

      f = 0.0_dp
      call adjoint_GL(x, f)

   end subroutine adjoint_rhs

   !-------------------------------------------------------------
   !-----     TYPE-BOUND PROCEDURES FOR THE GL OPERATOR     -----
   !-------------------------------------------------------------

   subroutine direct_matvec_GL(self, vec_in, vec_out)
      !> Linear Operator.
      class(GL_operator),          intent(inout)  :: self
      !> Input vector.
      class(abstract_vector_rdp), intent(in)  :: vec_in
      !> Output vector.
      class(abstract_vector_rdp), intent(out) :: vec_out
      ! internal
      character(len=*), parameter :: this_procedure = 'direct_matvec_GL'
      select type(vec_in)
      type is (state_vector)
         select type(vec_out)
         type is (state_vector)
            call direct_GL(vec_in%state, vec_out%state)
         class default
            call type_error('vec_out', 'state_vector', 'OUT', this_module, this_procedure)
         end select
      class default
         call type_error('vec_in', 'state_vector', 'IN', this_module, this_procedure)
      end select
   end subroutine direct_matvec_GL

   subroutine adjoint_matvec_GL(self, vec_in, vec_out)
      !> Linear Operator.
      class(GL_operator),         intent(inout)  :: self
      !> Input vector.
      class(abstract_vector_rdp), intent(in)  :: vec_in
      !> Output vector.
      class(abstract_vector_rdp), intent(out) :: vec_out
      character(len=*), parameter :: this_procedure = 'adjoint_matvec_GL'
      select type(vec_in)
      type is (state_vector)
         select type(vec_out)
         type is (state_vector)
            call adjoint_GL(vec_in%state, vec_out%state)
         class default
            call type_error('vec_out', 'state_vector', 'OUT', this_module, this_procedure)
         end select
      class default
         call type_error('vec_in', 'state_vector', 'IN', this_module, this_procedure)
      end select
   end subroutine adjoint_matvec_GL

   !------------------------------------------------------------------------
   !-----     TYPE-BOUND PROCEDURES FOR THE EXPONENTIAL PROPAGATOR     -----
   !------------------------------------------------------------------------

   subroutine direct_solver(self, vec_in, vec_out)
      ! Linear Operator.
      class(exponential_prop),     intent(inout)  :: self
      ! Input vector.
      class(abstract_vector_rdp),  intent(in)  :: vec_in
      ! Output vector.
      class(abstract_vector_rdp),  intent(out) :: vec_out

      ! Time-integrator.
      character(len=*), parameter :: this_procedure = 'direct_solver'
      type(rks54_class) :: prop
      real(dp)          :: dt = 1.0_dp

      select type(vec_in)
      type is(state_vector)
         select type(vec_out)
         type is(state_vector)

            ! Initialize propagator.
            call prop%initialize(n=2*nx, f=rhs)
            ! Integrate forward in time.
            call prop%integrate(0.0_dp, vec_in%state, dt, self%tau, vec_out%state)

         class default
            call type_error('vec_out', 'state_vector', 'OUT', this_module, this_procedure)
         end select
      class default
         call type_error('vec_in', 'state_vector', 'IN', this_module, this_procedure)
      end select
   end subroutine direct_solver

   subroutine adjoint_solver(self, vec_in, vec_out)
      ! Linear Operator.
      class(exponential_prop),     intent(inout)  :: self
      ! Input vector.
      class(abstract_vector_rdp),  intent(in)  :: vec_in
      ! Output vector.
      class(abstract_vector_rdp),  intent(out) :: vec_out

      ! Time-integrator.
      character(len=*), parameter :: this_procedure = 'adjoint_solver'
      type(rks54_class) :: prop
      real(dp)          :: dt = 1.0_dp

      select type(vec_in)
      type is(state_vector)
         select type(vec_out)
         type is(state_vector)

            ! Initialize propagator.
            call prop%initialize(n=2*nx, f=adjoint_rhs)
            ! Integrate forward in time.
            call prop%integrate(0.0_dp, vec_in%state, dt, self%tau, vec_out%state)

         class default
            call type_error('vec_out', 'state_vector', 'OUT', this_module, this_procedure)
         end select
      class default
         call type_error('vec_in', 'state_vector', 'IN', this_module, this_procedure)
      end select
   end subroutine adjoint_solver

   !--------------------------------------
   !-----     EXP(tA) SUBROUTINE     -----
   !--------------------------------------

   subroutine exptA(vec_out, A, vec_in, tau, info, trans)
      !! Subroutine for the exponential propagator that conforms with the abstract interface
      !! defined in expmlib.f90
      class(abstract_vector_rdp),  intent(out)   :: vec_out
      !! Output vector
      class(abstract_linop_rdp),   intent(inout) :: A
      !! Linear operator
      class(abstract_vector_rdp),  intent(in)    :: vec_in
      !! Input vector.
      real(dp),                    intent(in)    :: tau
      !! Integration horizon
      integer,                     intent(out)   :: info
      !! Information flag
      logical, optional,           intent(in)    :: trans
      logical                                    :: transpose
      !! Direct or Adjoint?
      ! internal
      character(len=*), parameter :: this_procedure = 'exptA'

      ! optional argument
      transpose = optval(trans, .false.)

      ! time integrator
      select type (vec_in)
      type is (state_vector)
         select type (vec_out)
         type is (state_vector)
            select type (A)
            type is (exponential_prop)
               ! set integration time
               A%tau = tau
               if (transpose) then
                  call A%rmatvec(vec_in, vec_out)
               else
                  call A%matvec(vec_in, vec_out)
               end if 
            class default
               call type_error('A', 'exponential_prop', this_module, this_procedure)
            end select
         class default
            call type_error('vec_out', 'state_vector', 'OUT', this_module, this_procedure)
         end select
      class default
         call type_error('vec_in', 'state_vector', 'IN', this_module, this_procedure)
      end select
   end subroutine exptA

   !--------------------------------------------------------
   !-----     TYPE BOUND PROCEDURES FOR LTI SYSTEMS    -----
   !--------------------------------------------------------

   subroutine initialize_lti_system(self, A, prop, B_in, CT_in, D)
      class(lti_system),           intent(inout) :: self
      class(abstract_linop_rdp),   intent(in)    :: A
      class(abstract_linop_rdp),   intent(in)    :: prop
      class(abstract_vector_rdp),  intent(in)    :: B_in(:)
      class(abstract_vector_rdp),  intent(in)    :: CT_in(:)
      real(dp),          optional, intent(in)    :: D(:,:)

      character(len=*), parameter :: this_procedure = 'initialize_lti_system'
      ! Operator
      allocate(self%A, source=A)
      ! Exp prop
      allocate(self%prop, source=prop)
      ! Input
      allocate(self%B(rk_b), source=B_in(:rk_b))
      ! Output
      allocate(self%CT(rk_c), source=CT_in(:rk_c))
      ! Throughput
      allocate(self%D(rk_c,rk_b))
      if (present(D)) then
         call assert_shape(D, [ rk_c, rk_b ], 'D', this_module, this_procedure)
         self%D = D
      else
         self%D = 0.0_dp
      end if
   end subroutine initialize_lti_system

end module Ginzburg_Landau_Operators