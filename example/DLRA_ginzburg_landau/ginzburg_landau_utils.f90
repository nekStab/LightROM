module Ginzburg_Landau_Utils
   ! Standard Library.
   use stdlib_math, only : linspace
   use stdlib_optval, only : optval
   use stdlib_linalg, only : eye, diag
   use stdlib_io_npy, only : save_npy, load_npy
   !use fortime
   ! LightKrylov for linear algebra.
   use LightKrylov
   use LightKrylov, only : wp => dp
   use LightKrylov_AbstractVectors
   use LightKrylov_Utils, only : svd, assert_shape
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
   ! mesh construction
   public  :: initialize_parameters
   ! utilities for state_vectors
   public  :: set_state, get_state, init_rand
   ! initial conditions
   public  :: generate_random_initial_condition
   ! logfiles
   public  :: stamp_logfile_header
   ! misc
   public  :: CALE, CARE, sval

   character*128, parameter :: this_module = 'Ginzburg_Landau_Utils'

contains

   !--------------------------------------------------------------
   !-----     CONSTRUCT THE MESH AND PHYSICAL PARAMETERS     -----
   !--------------------------------------------------------------

   subroutine initialize_parameters()
      implicit none
      ! Mesh array.
      real(wp), allocatable :: x(:)
      real(wp)              :: x2(1:2*nx)
      real(wp)              :: tmpv(N, 2)
      integer               :: i

      ! Construct mesh.
      x = linspace(-L/2, L/2, nx+2)
      dx = x(2)-x(1)

      ! Construct mu(x)
      mu(:) = (mu_0 - c_mu**2) + (mu_2 / 2.0_wp) * x(2:nx+1)**2

      ! Define integration weights
      weight     = dx
      weight_mat = dx

      ! Construct B & C
      ! B = [ [ Br, -Bi ], [ Bi, Br ] ]
      ! B = [ [ Cr, -Ci ], [ Ci, Cr ] ]
      ! where Bi = Ci = 0

      ! actuator is a Guassian centered just upstream of branch I
      ! column 1
      x2       = 0.0_wp
      x2(1:nx) = x(2:nx+1)
      B(1)%state = 0.5*exp(-((x2 - x_b)/s_b)**2)*sqrt(weight)
      ! column 2
      x2            = 0.0_wp
      x2(nx+1:2*nx) = x(2:nx+1)
      B(2)%state = 0.5*exp(-((x2 - x_b)/s_b)**2)*sqrt(weight)

      ! the sensor is a Gaussian centered at branch II
      ! column 1
      x2       = 0.0_wp
      x2(1:nx) = x(2:nx+1)
      CT(1)%state = 0.5*exp(-((x2 - x_c)/s_c)**2)*sqrt(weight)
      ! column 2
      x2            = 0.0_wp
      x2(nx+1:2*nx) = x(2:nx+1)
      CT(2)%state = 0.5*exp(-((x2 - x_c)/s_c)**2)*sqrt(weight)

      ! Note that we have included the integration weights into the actuator/sensor definitions

      ! RK lyap & riccati
      Qc   = eye(rk_c)
      Rinv = eye(rk_b)
      tmpv = 0.0_wp
      call get_state(tmpv(:,1:rk_b), B(1:rk_b))
      BBTW_flat(1:N**2)     = reshape(matmul(tmpv, transpose(tmpv)), shape(BBTW_flat))
      BRinvBTW_mat(1:N,1:N) = matmul(matmul(tmpv, Rinv), transpose(tmpv))
      call get_state(tmpv(:,1:rk_c), CT(1:rk_c))
      CTCW_flat(1:N**2)     = reshape(matmul(tmpv, transpose(tmpv)), shape(CTCW_flat))
      CTQcCW_mat(1:N,1:N)   = matmul(matmul(tmpv, Qc), transpose(tmpv))

      return
   end subroutine initialize_parameters

   !--------------------------------------------------------------------
   !-----     UTILITIES FOR STATE_VECTOR AND STATE MATRIX TYPES    -----
   !--------------------------------------------------------------------

   subroutine get_state(mat_out, state_in)
      !! Utility function to transfer data from a state vector to a real array
      real(wp),                   intent(out) :: mat_out(:,:)
      class(abstract_vector_rdp), intent(in)  :: state_in(:)
      ! internal variables
      integer :: k, kdim
      mat_out = 0.0_wp
      select type (state_in)
      type is (state_vector)
         kdim = size(state_in)
         call assert_shape(mat_out, (/ N, kdim /), 'get_state -> state_vector', 'mat_out')
         do k = 1, kdim
            mat_out(:,k) = state_in(k)%state
         end do
      type is (state_matrix)
         call assert_shape(mat_out, (/ N, N /), 'get_state -> state_matrix', 'mat_out')
         mat_out = reshape(state_in(1)%state, (/ N, N /))
      end select
      return
   end subroutine get_state

   subroutine set_state(state_out, mat_in)
      !! Utility function to transfer data from a real array to a state vector
      class(abstract_vector_rdp), intent(out) :: state_out(:)
      real(wp),                   intent(in)  :: mat_in(:,:)
      ! internal variables
      integer       :: k, kdim
      select type (state_out)
      type is (state_vector)
         kdim = size(state_out)
         call assert_shape(mat_in, (/ N, kdim /), 'set_state -> state_vector', 'mat_in')
         call zero_basis(state_out)
         do k = 1, kdim
            state_out(k)%state = mat_in(:,k)
         end do
      type is (state_matrix)
         call assert_shape(mat_in, (/ N, N /), 'set_state -> state_matrix', 'mat_in')
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

   !------------------------------------
   !-----     INITIAL CONDIIONS    -----
   !------------------------------------

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

   !-----------------------------
   !-----      LOGFILES     -----
   !-----------------------------

   subroutine stamp_logfile_header(iunit, problem, rk, tau, Tend, torder)
      integer,       intent(in) :: iunit
      character(*),  intent(in) :: problem
      integer,       intent(in) :: rk
      real(wp),      intent(in) :: tau
      real(wp),      intent(in) :: Tend
      integer,       intent(in) :: torder

      write(iunit,*) '-----------------------'
      write(iunit,*) '    GINZBURG LANDAU'
      write(iunit,*) '-----------------------'
      write(iunit,*) 'nu    = ', nu
      write(iunit,*) 'gamma = ', gamma
      write(iunit,*) 'mu_0  = ', mu_0
      write(iunit,*) 'c_mu  = ', c_mu
      write(iunit,*) 'mu_2  = ', mu_2
      write(iunit,*) '-----------------------'
      write(iunit,*) problem
      write(iunit,*) '-----------------------'
      write(iunit,*) 'nx    = ', nx
      write(iunit,*) 'rk_b  = ', rk_b
      write(iunit,*) 'x_b   = ', x_b
      write(iunit,*) 's_b   = ', s_b
      write(iunit,*) 'rk_c  = ', rk_c
      write(iunit,*) 'x_c   = ', x_c
      write(iunit,*) 's_c   = ', s_c
      write(iunit,*) '-----------------------'
      write(iunit,*) 'Time Integration: DLRA'
      write(iunit,*) '-----------------------'
      write(iunit,*) 'Tend   =', Tend
      write(iunit,*) 'torder =', torder
      write(iunit,*) 'tau    =', tau
      write(iunit,*) 'rk     =', rk
      write(iunit,*) '---------------------'
      write(iunit,*) '---------------------'
      return
   end subroutine stamp_logfile_header

   !-------------------------
   !-----      MISC     -----
   !-------------------------

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

   subroutine sval(X, svals)
      real(wp), intent(in) :: X(:,:)
      real(wp)             :: svals(min(size(X, 1), size(X, 2)))
      ! internals
      real(wp)             :: U(size(X, 1), size(X, 1))
      real(wp)             :: VT(size(X, 2), size(X, 2))
  
      ! Perform SVD
      call svd(X, U, svals, VT)
    
   end subroutine


end module Ginzburg_Landau_Utils