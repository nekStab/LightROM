module Ginzburg_Landau_Operators
   ! LightKrylov for linear algebra.
   use LightKrylov
   use LightKrylov, only : wp => dp
   use LightKrylov_AbstractVectors ! zero_basis
   use LightKrylov_utils, only : assert_shape, svd
   ! LightROM
   use LightROM_AbstractLTIsystems
   ! Ginzburg Landau
   use Ginzburg_Landau_Base
   ! RKLIB module for time integration.
   use rklib_module
   ! Standard Library.
   use stdlib_optval, only : optval
   use stdlib_linalg, only : diag
   implicit none

   private
   public :: exptA, direct_GL, adjoint_GL, sval, generate_random_initial_condition

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
      real(wp) :: tau ! Integration time.
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
      real(wp), dimension(:), intent(in)  :: vec_in
      !> Time-derivative.
      real(wp), dimension(:), intent(out) :: vec_out

      !> Internal variables.
      integer                 :: i
      real(wp), dimension(nx) :: u, v, du, dv
      real(wp)                :: d2u, d2v, cu, cv

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

      return
   end subroutine direct_GL

   !-----------------------------------------------------------
   !-----     Adjoint linear Ginzburg-Landau equation     -----
   !-----------------------------------------------------------

   subroutine adjoint_GL(vec_in, vec_out)
      !> State vector.
      real(wp), dimension(:), intent(in)  :: vec_in
      !> Time-derivative.
      real(wp), dimension(:), intent(out) :: vec_out

      ! Internal variables.
      integer :: i
      real(wp), dimension(nx) :: u, du
      real(wp), dimension(nx) :: v, dv
      real(wp)                :: d2u, d2v, cu, cv

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

      return
   end subroutine adjoint_GL

   !--------------------------------------
   !-----     WRAPPERS FOR RKLIB     -----
   !--------------------------------------

   subroutine rhs(me, t, x, f)
      ! Time-integrator.
      class(rk_class),               intent(inout) :: me
      ! Current time.
      real(wp),                      intent(in)    :: t
      ! State vector.
      real(wp),        dimension(:), intent(in)    :: x
      ! Time-derivative.
      real(wp),        dimension(:), intent(out)   :: f

      f = 0.0_wp
      call direct_GL(x, f)

      return
   end subroutine rhs

   subroutine adjoint_rhs(me, t, x, f)
      ! Time-integrator.
      class(rk_class),               intent(inout) :: me
      ! Current time.
      real(wp),                      intent(in)    :: t
      ! State vector.
      real(wp),        dimension(:), intent(in)    :: x
      ! Time-derivative.
      real(wp),        dimension(:), intent(out)   :: f

      f = 0.0_wp
      call adjoint_GL(x, f)

      return
   end subroutine adjoint_rhs

   !-------------------------------------------------------------
   !-----     TYPE-BOUND PROCEDURES FOR THE GL OPERATOR     -----
   !-------------------------------------------------------------

   subroutine direct_matvec_GL(self, vec_in, vec_out)
      !> Linear Operator.
      class(GL_operator),          intent(in)  :: self
      !> Input vector.
      class(abstract_vector_rdp), intent(in)  :: vec_in
      !> Output vector.
      class(abstract_vector_rdp), intent(out) :: vec_out
      select type(vec_in)
      type is (state_vector)
         select type(vec_out)
         type is (state_vector)
            call direct_GL(vec_in%state, vec_out%state)
         end select
      end select
      return
   end subroutine direct_matvec_GL

   subroutine adjoint_matvec_GL(self, vec_in, vec_out)
      !> Linear Operator.
      class(GL_operator),         intent(in)  :: self
      !> Input vector.
      class(abstract_vector_rdp), intent(in)  :: vec_in
      !> Output vector.
      class(abstract_vector_rdp), intent(out) :: vec_out
      select type(vec_in)
      type is (state_vector)
         select type(vec_out)
         type is (state_vector)
            call adjoint_GL(vec_in%state, vec_out%state)
         end select
      end select
      return
   end subroutine adjoint_matvec_GL

   !------------------------------------------------------------------------
   !-----     TYPE-BOUND PROCEDURES FOR THE EXPONENTIAL PROPAGATOR     -----
   !------------------------------------------------------------------------

   subroutine direct_solver(self, vec_in, vec_out)
      ! Linear Operator.
      class(exponential_prop),     intent(in)  :: self
      ! Input vector.
      class(abstract_vector_rdp),  intent(in)  :: vec_in
      ! Output vector.
      class(abstract_vector_rdp),  intent(out) :: vec_out

      ! Time-integrator.
      type(rks54_class) :: prop
      real(wp)          :: dt = 1.0_wp

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
      class(exponential_prop),     intent(in)  :: self
      ! Input vector.
      class(abstract_vector_rdp),  intent(in)  :: vec_in
      ! Output vector.
      class(abstract_vector_rdp),  intent(out) :: vec_out

      ! Time-integrator.
      type(rks54_class) :: prop
      real(wp)          :: dt = 1.0_wp

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
      class(abstract_vector_rdp),  intent(out)   :: vec_out
      !! Output vector
      class(abstract_linop_rdp),   intent(inout) :: A
      !! Linear operator
      class(abstract_vector_rdp),  intent(in)    :: vec_in
      !! Input vector.
      real(wp),                    intent(in)    :: tau
      !! Integration horizon
      integer,                     intent(out)   :: info
      !! Information flag
      logical, optional,           intent(in)    :: trans
      logical                                    :: transpose
      !! Direct or Adjoint?

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
               !A%tau = tau
               if (transpose) then
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

   subroutine initialize_lti_system(self, A, prop, B, CT, D)
      class(lti_system),           intent(inout) :: self
      class(abstract_linop_rdp),   intent(in)    :: A
      class(abstract_linop_rdp),   intent(in)    :: prop
      class(abstract_vector_rdp),  intent(in)    :: B(:)
      class(abstract_vector_rdp),  intent(in)    :: CT(:)
      real(wp),          optional, intent(in)    :: D(:,:)

      ! internal variables
      integer                                :: rk_b, rk_c

      ! Operator
      select type (A)
      type is (GL_operator)
         allocate(self%A, source=A)
      end select
      ! Exp prop
      select type (prop)
      type is (exponential_prop)
         allocate(self%prop, source=prop)
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

   subroutine sval(X, svals)
      real(wp), intent(in) :: X(:,:)
      real(wp)             :: svals(min(size(X, 1), size(X, 2)))
      ! internals
      real(wp)             :: U(size(X, 1), size(X, 1))
      real(wp)             :: VT(size(X, 2), size(X, 2))
  
      ! Perform SVD
      call svd(X, U, svals, VT)
    
   end subroutine

   subroutine generate_random_initial_condition(U, S, rk)
      class(state_vector),   intent(out) :: U(:)
      real(wp),              intent(out) :: S(:,:)
      integer,               intent(in)  :: rk
      ! internals
      class(state_vector),   allocatable :: Utmp(:)
      integer,               allocatable :: perm(:)
      ! SVD
      real(wp)                           :: U_svd(rk,rk)
      real(wp)                           :: S_svd(rk)
      real(wp)                           :: V_svd(rk,rk)
      integer                            :: i, info

      if (size(U) < rk) then
         write(*,*) 'Input krylov basis size incompatible with requested rank', rk
         STOP 1
      else
         call zero_basis(U)
         do i = 1,rk
            call U(i)%rand(.false.)
         end do
      end if
      if (size(S,1) < rk) then
         write(*,*) 'Input coefficient matrix size incompatible with requested rank', rk
         STOP 1
      else if (size(S,1) /= size(S,2)) then
         write(*,*) 'Input coefficient matrix must be square.'
         STOP 2
      else
         S = 0.0_wp
      end if
      ! perform QR
      allocate(perm(1:rk)); perm = 0
      allocate(Utmp(1:rk), source=U(1:rk))
      call qr(Utmp, S, perm, info, verbosity=.false.)
      if (info /= 0) write(*,*) '  [generate_random_initial_condition] Info: Colinear vectors detected in QR, column ', info
      ! perform SVD
      call svd(S(:,1:rk), U_svd(:,1:rk), S_svd(1:rk), V_svd(1:rk,1:rk))
      S = diag(S_svd)
      block
         class(abstract_vector_rdp), allocatable :: Xwrk(:)
         call linear_combination(Xwrk, Utmp, U_svd)
         call copy_basis(U, Xwrk)
      end block
      
   end subroutine

end module Ginzburg_Landau_Operators