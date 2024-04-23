module Ginzburg_Landau_RKlib
   ! LightKrylov for linear algebra.
   use LightKrylov
   use LightKrylov_utils, only : assert_shape
   ! LightROM
   use LightROM_AbstractLTIsystems
   ! Ginzburg Landau
   use Ginzburg_Landau_Base
   ! RKLIB module for time integration.
   use rklib_module
   ! Standard Library.
   use stdlib_optval, only : optval
   implicit none

   private
   public :: exptA

   !-----------------------------------------------
   !-----     LIGHTKRYLOV LTI SYSTEM TYPE     -----
   !-----------------------------------------------

   type, extends(abstract_lti_system), public :: lti_system
   contains
      private
      procedure, pass(self), public :: initialize_lti_system
   end type lti_system

   !------------------------------------------
   !-----     EXPONENTIAL PROPAGATOR     -----
   !------------------------------------------

   type, extends(abstract_linop), public :: exponential_prop
      real(kind=wp) :: tau ! Integration time.
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

   subroutine rhs(me, t, x, f)
      ! Time-integrator.
      class(rk_class), intent(inout)             :: me
      ! Current time.
      real(kind=wp)  , intent(in)                :: t
      ! State vector.
      real(kind=wp)  , dimension(:), intent(in)  :: x
      ! Time-derivative.
      real(kind=wp)  , dimension(:), intent(out) :: f

      ! Internal variables.
      integer :: i, j, k
      real(kind=wp), dimension(nx) :: u, du
      real(kind=wp), dimension(nx) :: v, dv
      real(kind=wp)                :: d2u, d2v, cu, cv

      ! Sets the internal variables.
      f = 0.0_wp
      u = x(1:nx)      ; du = f(1:nx)
      v = x(nx+1:2*nx) ; dv = f(nx+1:2*nx)

      !---------------------------------------------------
      !-----     Linear Ginzburg Landau Equation     -----
      !---------------------------------------------------

      ! Left most boundary points.
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

      ! Copy results to the output array.
      f(1:nx) = du ; f(nx+1:2*nx) = dv

      return
   end subroutine rhs

   !-----------------------------------------------------------
   !-----     Adjoint linear Ginzburg-Landau equation     -----
   !-----------------------------------------------------------

   subroutine adjoint_rhs(me, t, x, f)
      ! Time-integrator.
      class(rk_class), intent(inout)             :: me
      ! Current time.
      real(kind=wp)  , intent(in)                :: t
      ! State vector.
      real(kind=wp)  , dimension(:), intent(in)  :: x
      ! Time-derivative.
      real(kind=wp)  , dimension(:), intent(out) :: f

      ! Internal variables.
      integer :: i, j, k
      real(kind=wp), dimension(nx) :: u, du
      real(kind=wp), dimension(nx) :: v, dv
      real(kind=wp)                :: d2u, d2v, cu, cv

      ! Sets the internal variables.
      f = 0.0_wp
      u = x(1:nx)      ; du = f(1:nx)
      v = x(nx+1:2*nx) ; dv = f(nx+1:2*nx)

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
      f(1:nx) = du ; f(nx+1:2*nx) = dv

      return
   end subroutine adjoint_rhs

   !------------------------------------------------------------------------
   !-----     TYPE-BOUND PROCEDURES FOR THE EXPONENTIAL PROPAGATOR     -----
   !------------------------------------------------------------------------

   subroutine direct_solver(self, vec_in, vec_out)
      ! Linear Operator.
      class(exponential_prop), intent(in)  :: self
      ! Input vector.
      class(abstract_vector) , intent(in)  :: vec_in
      ! Output vector.
      class(abstract_vector) , intent(out) :: vec_out

      ! Time-integrator.
      type(rks54_class) :: prop
      real(kind=wp)     :: dt = 1.0_wp

      select type(vec_in)
      type is(state_vector)
         select type(vec_out)
         type is(state_vector)

            ! Initialize propagator.
            call prop%initialize(n=2*nx, f=rhs)
            ! Integrate forward in time.
            call prop%integrate(0.0_wp, vec_in%state, dt, self%tau, vec_out%state)

         end select
      end select
      return
   end subroutine direct_solver

   subroutine adjoint_solver(self, vec_in, vec_out)
      ! Linear Operator.
      class(exponential_prop), intent(in)  :: self
      ! Input vector.
      class(abstract_vector) , intent(in)  :: vec_in
      ! Output vector.
      class(abstract_vector) , intent(out) :: vec_out

      ! Time-integrator.
      type(rks54_class) :: prop
      real(kind=wp)     :: dt = 1.0_wp

      select type(vec_in)
      type is(state_vector)
         select type(vec_out)
         type is(state_vector)

            ! Initialize propagator.
            call prop%initialize(n=2*nx, f=adjoint_rhs)
            ! Integrate forward in time.
            call prop%integrate(0.0_wp, vec_in%state, dt, self%tau, vec_out%state)

         end select
      end select
      return
   end subroutine adjoint_solver

   !--------------------------------------
   !-----     EXP(tA) SUBROUTINE     -----
   !--------------------------------------

   subroutine exptA(vec_out, A, vec_in, tau, info, trans)
      !! Subroutine for the exponential propagator that conforms with the abstract interface
      !! defined in expmlib.f90
      class(abstract_vector),  intent(out)   :: vec_out
      !! Output vector
      class(abstract_linop),   intent(inout) :: A
      !! Linear operator
      class(abstract_vector),  intent(in)    :: vec_in
      !! Input vector.
      real(kind=wp),           intent(in)    :: tau
      !! Integration horizon
      integer,                 intent(out)   :: info
      !! Information flag
      logical, optional,       intent(in)    :: trans
      !! Direct or Adjoint?                
      
      ! internal variables
      logical                                :: transpose

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
               if (trans) then
                  call A%rmatvec(vec_in, vec_out)
               else
                  call A%matvec(vec_in, vec_out)
               end if 
            end select
         end select
      end select

   end subroutine exptA

   !--------------------------------------------------------
   !-----     TYPE BOUND PROCEDURES FOR LTI SYSTEMS    -----
   !--------------------------------------------------------

   subroutine initialize_lti_system(self, A, B, CT, D)
      class(lti_system),       intent(inout) :: self
      class(abstract_linop),   intent(in)    :: A
      class(abstract_vector),  intent(in)    :: B(:)
      class(abstract_vector),  intent(in)    :: CT(:)
      real(kind=wp), optional, intent(in)    :: D(:,:)

      ! internal variables
      integer                                :: rk_b, rk_c

      ! Operator
      select type (A)
      type is (exponential_prop)
         allocate(self%A, source=A)
      end select
      ! Input
      select type (B)
      type is (state_vector)
         rk_b = size(B)
         allocate(self%B(1:rk_b), source=B(1:rk_b))
      end select
      ! Output
      select type (CT)
         type is (state_vector)
         rk_c = size(CT)
         allocate(self%CT(1:rk_c), source=CT(1:rk_c))
      end select
      ! Throughput
      allocate(self%D(1:rk_c, 1:rk_b))
      if (present(D)) then
         call assert_shape(D, (/ rk_c, rk_b /), 'initialize_lti_system', 'D')
         self%D = D
      else
         self%D = 0.0_wp
      end if
      return
   end subroutine initialize_lti_system

end module Ginzburg_Landau_RKlib