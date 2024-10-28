module Laplacian2D_LTI_Riccati_Utils
   ! Standard Library.
   use stdlib_math, only : linspace
   use stdlib_optval, only : optval
   use stdlib_linalg, only : eye, diag, svd, inv
   ! RKLIB module for time integration.
   use rklib_module
   ! LightKrylov for linear algebra.
   use LightKrylov
   use LightKrylov, only : wp => dp
   use LightKrylov_AbstractVectors ! zero_basis
   use LightKrylov_Utils, only : assert_shape
   ! Laplacian
   use Laplacian2D_LTI_Riccati_Base
   use laplacian2D_LTI_Riccati_Operators
   implicit none

   private :: this_module
   ! initialisation
   public  :: initialize_problem
   ! utils for state vector/matrix
   public  :: get_state, set_state, init_rand
   ! initial condition
   public  :: generate_random_initial_condition
   ! misc
   public  :: CARE

   character(len=128), parameter :: this_module = 'Laplacian2D_LTI_Riccati_Utils'
   
contains

   !---------------------------------------
   !-----     CONSTRUCT THE MESH      -----
   !---------------------------------------

   subroutine initialize_problem(magQ, magR)
      implicit none
      real(wp), intent(in)  :: magQ, magR
      ! internals
      real(wp), allocatable :: x(:)
      real(wp)              :: Bwrk(N,rk_b)
      integer :: i

      !> Construct mesh.
      x = linspace(-L/2, L/2, nx)

      ! Define C, Qc & compute CTQcC
      Qc = magQ*eye(rk_c)

      call init_rand(CT, ifnorm = .false.)
      call get_state(CTdata, CT, 'Extract CT')
      CTQcCdata = matmul(CTdata, matmul( Qc, transpose(CTdata)))
      CTQcC(:N**2) = reshape(CTQcCdata, shape(CTQcC))

      ! Define B, Rinv & compule BRinvBT
      if (magR .lt. atol_dp) then
         Rinv = 0.0_wp
      else
         Rinv = 1/magR*eye(rk_b)
      endif

      call init_rand(B, ifnorm = .false.)
      Bdata = 0.0_wp
      call get_state(Bdata, B, 'Extract B')
      Bwrk = 0.0_wp
      Bwrk = matmul(Bdata, Rinv)
      BRinvBTdata = matmul( Bwrk, transpose(Bdata) )

      return
   end subroutine initialize_problem

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
      integer :: k, kdim
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

   !--------------------------------------
   !-----     INITIAL CONDITIONS     -----
   !--------------------------------------

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
      character(len=128) :: msg

      if (size(U) < rk) then
         write(msg,'(A,I0)') 'Input krylov basis size incompatible with requested rank ', rk
         call stop_error(msg, module=this_module, procedure='generate_random_initial_condition')
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

   !------------------------
   !-----     MISC     -----
   !------------------------

   function CARE(X,A,Q,BRinvBT) result(res)
      real(wp), dimension(n,n) :: X, A, Q, BRinvBT, res
      res = matmul(A, X) + matmul(X, transpose(A)) + Q - matmul(X, matmul(BRinvBT, X))
   end function CARE

   subroutine solve_riccati(X, A, BRinvBT, CTQC)
      !! Solve the Lyapunov equation directly
      real(wp), intent(out) :: X(N,N)
      !! Solution
      real(wp), intent(in)  :: A(N,N)
      !! Operator
      real(wp), intent(in)  :: BRinvBT(N,N)
      !! Nonlinearity
      real(wp), intent(in)  :: CTQC(N,N)
      !! Inhomogeneity
      ! Internal
      real(wp), dimension(2*N,2*N) :: H, VR
      real(wp), dimension(2*N,N)   :: UR, UI
      real(wp), dimension(N,N)     :: F, Ginv
      real(wp), dimension(2*N)     :: wr, wi
      integer, parameter           :: lwork = 1040
      real(wp)                     :: work(lwork)
      logical                      :: flag
      integer                      :: icnt, i, info

      ! construct Hmat
      H = 0.0_wp
      H(  1:N,    1:N  ) =  A
      H(N+1:2*N,N+1:2*N) = -transpose(A)
      H(  1:N,  N+1:2*N) = -BRinvBT
      H(N+1:2*N,  1:N  ) = -CTQC

      call dgeev('N', 'V', 2*N, H, 2*N, wr, wi, VR, 2*N, VR, 2*N, work, lwork, info)

      ! extract stable eigenspace
      UR = 0.0_wp
      UI = 0.0_wp
      icnt = 0
      flag = .true.
      do i = 1, 2*N
         if ( wr(i) .lt. 0.0 ) then
            icnt = icnt + 1
            if ( wi(i) .eq. 0.0 ) then ! real
               UR(:,icnt) = VR(:,i)
            else                       ! complex
               if (flag) then
                  UR(:,icnt)   =  VR(:,i)
                  UI(:,icnt)   =  VR(:,i+1)
                  UR(:,icnt+1) =  VR(:,i)
                  UI(:,icnt+1) = -VR(:,i+1)
                  flag = .false.
               else
                  flag = .true.
               end if
            end if
         end if
      end do
      ! construct solution
      F    = matmul(UR(n+1:2*N,:), transpose(UR(1:N,:))) + matmul(UI(n+1:2*n,:), transpose(UI(1:N,:)))
      Ginv = matmul(UR(  1:N,  :), transpose(UR(1:N,:))) + matmul(UI(  1:N,  :), transpose(UI(1:N,:)))
      X    = matmul(F, inv(Ginv))

      return
   end subroutine solve_riccati

end module Laplacian2D_LTI_Riccati_Utils