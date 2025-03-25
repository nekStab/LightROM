module LightROM_TestUtils
    ! Stdlib
    use stdlib_io_npy, only: save_npy
    use stdlib_linalg, only: eye, diag
    use stdlib_math, only : linspace
    use stdlib_stats_distribution_normal, only: normal => rvs_normal
    use stdlib_optval, only: optval
    ! RKLIB module for time integration.
    use rklib_module
    ! LightKrylov
    use LightKrylov
    use LightKrylov_Utils
    use LightKrylov_Logger
    use LightKrylov_Constants
    use LightKrylov_TestUtils
    
    implicit none
    
    private

    character(len=128), parameter, private :: this_module = 'LightROM_TestUtils'

    public :: initialize_GL_parameters
    public :: solve_lyapunov
    public :: get_state

    !---------------------------------------------------
    !-----     GINZBURG-LANDAU TYPE DEFINITION     -----
    !---------------------------------------------------

    ! Mesh related parameters.
    real(dp), parameter :: L  = 50.0_dp ! Domain length
    integer,  parameter :: nx = 128       ! Number of grid points (excluding boundaries).
    integer,  parameter :: N = 2*nx
    real(dp)            :: dx           ! Grid size.

    ! Physical parameters.
    complex(dp), parameter :: nu    = cmplx(2.0_dp, 0.2_dp, kind=dp)
    complex(dp), parameter :: gamma = cmplx(1.0_dp, -1.0_dp, kind=dp)
    real(dp),    parameter :: mu_0  = 0.38_dp
    real(dp),    parameter :: c_mu  = 0.2_dp
    real(dp),    parameter :: mu_2  = -0.01_dp
    real(dp)               :: mu(nx)

    ! Input parameters
    real(dp)               :: weight(N)       ! integration weights
    real(dp), parameter    :: x_b = -11.0_dp     ! location of input Gaussian
    real(dp), parameter    :: s_b = 1.0_dp       ! variance of input Gaussian

    !-------------------------------------------
    !-----     LIGHTKRYLOV VECTOR TYPE     -----
    !-------------------------------------------

    type, extends(abstract_vector_rdp), public :: state_vector
        real(dp) :: state(N) = 0.0_dp
    contains
        private
        procedure, pass(self), public :: zero => zero_state
        procedure, pass(self), public :: dot => dot_state
        procedure, pass(self), public :: scal => scal_state
        procedure, pass(self), public :: axpby => axpby_state
        procedure, pass(self), public :: rand => rand_state
        procedure, pass(self), public :: get_size => get_size_state
    end type state_vector

    !------------------------------------------
    !-----     EXPONENTIAL PROPAGATOR     -----
    !------------------------------------------

    type, extends(abstract_linop_rdp), public :: GL_exponential_prop
        real(dp), public :: tau ! Integration time.
    contains
        private
        procedure, pass(self), public :: matvec => GL_direct_solver
        procedure, pass(self), public :: rmatvec => GL_direct_solver ! dummy
    end type GL_exponential_prop

contains

    !---------------------------------------------------
    !-----     GINZBURG-LANDAU TYPE DEFINITION     -----
    !---------------------------------------------------

    subroutine initialize_GL_parameters(X0, A, BBT)
        type(state_vector), allocatable, intent(out) :: X0(:)
        real(dp), allocatable, intent(out) :: A(:,:)
        real(dp), allocatable, intent(out) :: BBT(:,:)
        ! internal
        real(dp), allocatable :: x(:)
        real(dp)              :: x2(1:N)
        real(dp)              :: X0mat(N,2)
        integer               :: i

        ! Construct mesh.
        x = linspace(-L/2, L/2, nx+2)
        dx = x(2)-x(1)

        ! Construct mu(x)
        mu(:) = (mu_0 - c_mu**2) + (mu_2 / 2.0_dp) * x(2:nx+1)**2

        ! Define integration weights
        weight = dx

        ! Construct Impulse
        ! X0 = [ [ X0r, -X0i ], [ X0i, X0r ] ]
        ! where X0i = 0

        ! The Impulse is a Guassian centered just upstream of branch I
        allocate(X0(2))
        ! column 1
        x2       = 0.0_dp
        x2(1:nx) = x(2:nx+1)
        X0(1)%state = exp(-((x2 - x_b)/s_b)**2)
        ! column 2
        x2            = 0.0_dp
        x2(nx+1:N) = x(2:nx+1)
        X0(2)%state = exp(-((x2 - x_b)/s_b)**2)

        allocate(BBT(N,N))
        call get_state(X0mat, X0)
        BBT = matmul(X0mat, dx*transpose(X0mat))

        ! Build the GL opreator
        allocate(A(N,N))
        do i = 1, N
            x2 = 0.0_dp
            x2(i) = 1.0_dp
            call GL_operator(x2, A(:,i))
        end do

    end subroutine initialize_GL_parameters
    
    !---------------------------------------------------------
    !-----      LINEARIZED GINZBURG-LANDAU EQUATIONS     -----
    !---------------------------------------------------------

    subroutine GL_operator(vec_in, vec_out)

        !> State vector.
        real(dp), dimension(:), intent(in)  :: vec_in
        !> Time-derivative.
        real(dp), dimension(:), intent(out) :: vec_out
  
        !> Internal variables.
        integer                 :: i
        real(dp), dimension(nx) :: u, v, du, dv
        real(dp)                :: d2u, d2v, cu, cv
  
        u = vec_in(1:nx)     
        v = vec_in(nx+1:N)
  
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
        vec_out(nx+1:N) = dv
  
    end subroutine GL_operator
    
    subroutine GL_rhs(me, t, x, f)
        ! Time-integrator.
        class(rk_class), intent(inout)             :: me
        ! Current time.
        real(dp), intent(in)                :: t
        ! State vector.
        real(dp), dimension(:), intent(in)  :: x
        ! Time-derivative.
        real(dp), dimension(:), intent(out) :: f
        
        f = 0.0_dp
        call GL_operator(x, f)

    end subroutine GL_rhs

    !----------------------------------------------------
    !-----     TYPE-BOUND PROCEDURE FOR VECTORS     -----
    !----------------------------------------------------

    subroutine zero_state(self)
        class(state_vector), intent(inout) :: self
        self%state = 0.0_dp
    end subroutine zero_state

    real(dp) function dot_state(self, vec) result(alpha)
        class(state_vector), intent(in) :: self
        class(abstract_vector_rdp), intent(in) :: vec
        select type (vec)
        type is (state_vector)
            alpha = dot_product(self%state, weight*vec%state)
        class default
            call stop_error("The intent [IN] argument 'vec' must be of type 'state_vector'", this_module, 'dot')
        end select
    end function dot_state

    subroutine scal_state(self, alpha)
        class(state_vector), intent(inout) :: self
        real(dp),            intent(in)    :: alpha
        self%state = self%state * alpha
    end subroutine scal_state

    subroutine axpby_state(self, alpha, vec, beta)
        class(state_vector),        intent(inout) :: self
        class(abstract_vector_rdp), intent(in)    :: vec
        real(dp),                   intent(in)    :: alpha, beta
        select type(vec)
        type is(state_vector)
            self%state = alpha*self%state + beta*vec%state
        class default
            call stop_error("The intent [IN] argument 'vec' must be of type 'state_vector'", this_module, 'axpby')
        end select
    end subroutine axpby_state

    integer function get_size_state(self) result(N)
        class(state_vector), intent(in) :: self
        N = N
    end function get_size_state

    subroutine rand_state(self, ifnorm)
        class(state_vector), intent(inout) :: self
        logical, optional, intent(in) :: ifnorm
        real(dp) :: tmp(nx, 2)
        ! internals
        logical :: normalize
        real(dp) :: alpha
        real(dp), dimension(N) :: mean, std
        normalize = optval(ifnorm,.true.)
        mean = 0.0_dp
        std  = 1.0_dp
        self%state = normal(mean,std)
        if (normalize) then
            alpha = self%norm()
            call self%scal(1.0/alpha)
        endif
    end subroutine rand_state

    !------------------------------------------------------------------------
    !-----     TYPE-BOUND PROCEDURES FOR THE EXPONENTIAL PROPAGATOR     -----
    !------------------------------------------------------------------------

    subroutine GL_direct_solver(self, vec_in, vec_out)
        ! Linear Operator.
      class(GL_exponential_prop),  intent(inout)  :: self
      ! Input vector.
      class(abstract_vector_rdp),  intent(in)  :: vec_in
      ! Output vector.
      class(abstract_vector_rdp),  intent(out) :: vec_out

      ! Time-integrator.
      type(rks54_class) :: prop
      real(dp)          :: dt = 1.0_dp

      select type(vec_in)
      type is(state_vector)
         select type(vec_out)
         type is(state_vector)

            ! Initialize propagator.
            call prop%initialize(n=N, f=GL_rhs)
            ! Integrate forward in time.
            call prop%integrate(0.0_dp, vec_in%state, dt, self%tau, vec_out%state)

         class default
            call stop_error('vec_out must be a state_vector', this_module, 'direct_solver')
         end select
      class default
         call stop_error('vec_in must be a state_vector', this_module, 'direct_solver')
      end select
    end subroutine GL_direct_solver

    subroutine reconstruct_TQ(T, Q, A, D, E, tw)
        !! Reconstruct tridiagonal matrix T and orthogonal projector Q from dsytd2 output (A, D, E)
        real(dp), intent(out) :: T(N,N)
        real(dp), intent(out) :: Q(N,N)
        real(dp), intent(in)  :: A(N,N)
        real(dp), intent(in)  :: D(N)
        real(dp), intent(in)  :: E(N-1)
        real(dp), intent(in)  :: tw(N-1)
  
        ! internal variables
        real(dp)  :: Hi(N,N)
        real(dp)  :: vec(N,1)
        integer :: i
  
        ! Build orthogonal Q = H(1) @  H(2) @ ... @ H(n-1)
        Q = eye(N)
        do i = 1, N - 1
           vec          = 0.0_dp
           vec(i+1,1)   = 1.0_dp
           vec(i+2:N,1) = A(i+2:N,i)
           Hi           = eye(N) - tw(i) * matmul( vec, transpose(vec) )
           Q            = matmul( Q, Hi )
        end do
  
        ! Build tridiagonal T
        T = 0.0_dp
        do i = 1, N
           T(i,i) = D(i)
        end do
        do i = 1, N - 1
           T(i,i+1) = E(i)
           T(i+1,i) = E(i)
        end do
  
    end subroutine reconstruct_TQ
  
    subroutine solve_lyapunov(X, A, P)
        !! Solve the Lyapunov equation directly
        real(dp), allocatable, intent(out) :: X(:,:)
        !! Solution
        real(dp), intent(in)  :: A(N,N)
        !! Operator
        real(dp), intent(in)  :: P(N,N)
        !! Inhomogeneity
        ! Internal
        real(dp), dimension(N,N) :: T, Q, Z, V, W, Y
        real(dp), dimension(N)   :: Dm, wrk, wr, wi
        real(dp), dimension(N-1) :: E, tw
        real(dp)                 :: scale
        integer                  :: isgn, info

        allocate(X(N,N))
  
        ! Transform operator to tridiagonal form
        call dsytd2('L', N, A, N, Dm, E, tw, info)
  
        ! Reconstruct T and Q
        call reconstruct_TQ(T, Q, A, Dm, E, tw)
  
        ! compute real Schur form A = Z @ T @ Z.T
        call dhseqr('S', 'I', N, 1, N, T, N, wr, wi, Z, N, wrk, N, info )
  
        ! Change RHS Basis: base --> Q --> Z
        V = matmul(transpose(Q), matmul(-P, Q))
        W = matmul(transpose(Z), matmul( V, Z))
  
        ! Compute solution of Lyapunov equation for Schur decomposition
        isgn = 1; scale = 0.1_dp
        call dtrsyl('N', 'T', isgn, N, N, T, N, T, N, W, N, scale, info)
  
        ! Return to original basis to obtain X_ref: Z --> Q --> base
        Y = matmul(Z, matmul(W, transpose(Z)))
        X = matmul(Q, matmul(Y, transpose(Q)))
  
    end subroutine solve_lyapunov

    !--------------------------------------------------------------------
    !-----     UTILITIES FOR STATE_VECTOR AND STATE MATRIX TYPES    -----
    !--------------------------------------------------------------------

    subroutine get_state(mat_out, state_in)
        !! Utility function to transfer data from a state vector to a real array
        real(dp), intent(out) :: mat_out(:,:)
        type(state_vector), intent(in)  :: state_in(:)
        ! internal variables
        integer :: k, kdim
        mat_out = 0.0_dp
        kdim = size(state_in)
        call assert_shape(mat_out, [ N, kdim ], 'mat_out', this_module, 'get_state')
        do k = 1, kdim
            mat_out(:,k) = state_in(k)%state
        end do
    end subroutine get_state

end module
