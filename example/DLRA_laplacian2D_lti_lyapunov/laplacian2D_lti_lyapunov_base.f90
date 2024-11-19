module laplacian2D_LTI_Lyapunov_Base
   ! Standard Library.
   use stdlib_stats_distribution_normal, only: normal => rvs_normal
   use stdlib_io_npy, only: save_npy
   use stdlib_strings, only: replace_all
   use stdlib_linalg, only: svdvals
   use stdlib_optval, only : optval
   ! LightKrylov for linear algebra.
   use LightKrylov
   use LightKrylov, only : wp => dp
   use LightKrylov_Logger
   use LightKrylov_Utils, only : assert_shape
   use LightKrylov_BaseKrylov, only : orthogonalize_against_basis, orthonormalize_basis
   use LightKrylov_AbstractVectors ! zero_basis
   ! LightROM
   use LightROM_AbstractLTIsystems ! LR_state
   implicit none

   private :: this_module
   character*128, parameter :: this_module = 'Laplacian2D_LTI_Lyapunov_Base'
   
   ! problem parameters
   public  :: N, nx, dx, dx2, L, rk_b, B, BBT

   !------------------------------
   !-----     PARAMETERS     -----
   !------------------------------

   ! --> Mesh related parameters.
   real(wp),      parameter :: L  = 1.0_wp  !> Domain length
   integer,       parameter :: nx = 4       !> Number of grid points per direction
   integer,       parameter :: N  = nx**2   !> total number of grid points
   real(wp),      parameter :: dx = L/nx    !> Grid size.
   real(wp),      parameter :: dx2= dx**2   !> Grid size.
   integer,       parameter :: rk_b = 1     !> rank of the RHS
   integer,       parameter :: rk_c = 1     !

   !-------------------------------------------------------
   !-----     LIGHTKRYLOV SYM LOW RANK STATE TYPE     -----
   !-------------------------------------------------------

   type, extends(abstract_sym_low_rank_state_rdp), public :: LR_state
   contains
      private
      procedure, pass(self), public :: init => initialize_LR_state
   end type LR_state

   !-------------------------------------------
   !-----     LIGHTKRYLOV VECTOR TYPE     -----
   !-------------------------------------------

   type, extends(abstract_vector_rdp), public :: state_vector
      real(wp) :: state(N) = 0.0_wp
      contains
      private
      procedure, pass(self), public :: zero => vector_zero
      procedure, pass(self), public :: rand => vector_rand
      procedure, pass(self), public :: scal => vector_scal
      procedure, pass(self), public :: axpby => vector_axpby
      procedure, pass(self), public :: dot => vector_dot
      procedure, pass(self), public :: get_size => vector_get_size
   end type state_vector

   !-------------------------------------------
   !-----     LIGHTKRYLOV VECTOR TYPE     -----
   !-------------------------------------------

   type, extends(abstract_vector_rdp), public :: state_matrix
      real(wp)      :: state(N**2) = 0.0_wp
   contains
      private
      procedure, pass(self), public :: zero => matrix_zero
      procedure, pass(self), public :: rand => matrix_rand
      procedure, pass(self), public :: scal => matrix_scal
      procedure, pass(self), public :: axpby => matrix_axpby
      procedure, pass(self), public :: dot => matrix_dot
      procedure, pass(self), public :: get_size => matrix_get_size  
   end type state_matrix

   type(state_vector)       :: B(rk_b)
   real(wp)                 :: BBT(N**2)

contains

   !-----     TYPE-BOUND PROCEDURE FOR VECTORS     -----

   subroutine vector_zero(self)
      class(state_vector), intent(inout) :: self
      self%state = 0.0_wp
      return
   end subroutine vector_zero

   real(wp) function vector_dot(self, vec) result(alpha)
      class(state_vector)       , intent(in) :: self
      class(abstract_vector_rdp), intent(in) :: vec
      select type(vec)
      type is (state_vector)
         alpha = dot_product(self%state, vec%state)
      end select
      return
   end function vector_dot

   integer function vector_get_size(self) result(ntot)
     class(state_vector), intent(in) :: self
     ntot = nx
     return
   end function vector_get_size

   subroutine vector_scal(self, alpha)
      class(state_vector), intent(inout) :: self
      real(wp)           , intent(in)    :: alpha
      self%state = self%state * alpha
      return
   end subroutine vector_scal

   subroutine vector_axpby(self, alpha, vec, beta)
      class(state_vector)       , intent(inout) :: self
      class(abstract_vector_rdp), intent(in)    :: vec
      real(wp)                  , intent(in)    :: alpha, beta
      select type(vec)
      type is (state_vector)
         self%state = alpha*self%state + beta*vec%state
      end select
      return
   end subroutine vector_axpby
   
   subroutine vector_rand(self, ifnorm)
      class(state_vector), intent(inout) :: self
      logical, optional,   intent(in)    :: ifnorm
      ! internals
      logical :: normalize
      real(wp) :: alpha
      real(wp), dimension(N) :: mean, std
      normalize = optval(ifnorm,.true.)
      mean = 0.0_wp
      std  = 1.0_wp
      self%state = normal(mean,std)
      if (normalize) then
         alpha = self%norm()
         call self%scal(1.0/alpha)
      endif
      return
   end subroutine vector_rand

   !-----     TYPE-BOUND PROCEDURE FOR MATRICES     -----

   subroutine matrix_zero(self)
      class(state_matrix), intent(inout) :: self
      self%state = 0.0_wp
      return
   end subroutine matrix_zero

   real(wp) function matrix_dot(self, vec) result(alpha)
      class(state_matrix)       , intent(in) :: self
      class(abstract_vector_rdp), intent(in) :: vec
      select type(vec)
      type is(state_matrix)
          alpha = dot_product(self%state, vec%state)
      end select
      return
   end function matrix_dot

   integer function matrix_get_size(self) result(N)
     class(state_matrix), intent(in) :: self
     N = N
     return
   end function matrix_get_size

   subroutine matrix_scal(self, alpha)
      class(state_matrix), intent(inout) :: self
      real(wp)           , intent(in)    :: alpha
      self%state = self%state * alpha
      return
   end subroutine matrix_scal  

   subroutine matrix_axpby(self, alpha, vec, beta)
      class(state_matrix)       , intent(inout) :: self
      class(abstract_vector_rdp), intent(in)    :: vec
      real(wp)                  , intent(in)    :: alpha, beta
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

   !-----------------------------------------------------------------------
   !-----     TYPE BOUND PROCEDURE FOR SYM LOW RANK REPRESENTATION    -----
   !-----------------------------------------------------------------------

   subroutine initialize_LR_state(self, U, S, rk, rkmax, if_rank_adaptive, casename, outpost)
      class(LR_state),            intent(inout) :: self
      class(abstract_vector_rdp), intent(in)    :: U(:)
      real(wp),                   intent(in)    :: S(:,:)
      integer,                    intent(in)    :: rk
      integer, optional,          intent(in)    :: rkmax
      logical, optional,          intent(in)    :: if_rank_adaptive
      logical                                   :: ifrk
      character(len=128), optional, intent(in)  :: casename
      procedure(abstract_outpost_rdp), optional :: outpost

      ! internals
      class(abstract_vector_rdp), allocatable   :: Utmp(:)
      real(wp), allocatable :: R(:, :)
      integer :: i, m, rka, info
      character(len=128) :: msg

      ifrk = optval(if_rank_adaptive, .false.)

      ! reset time
      self%time = 0.0_wp

      select type (U)
      type is (state_vector)
         ! set time and optional args
         self%time = 0.0_wp
         if (present(outpost)) self%outpost => outpost
         self%casename = optval(casename, '')

         m = size(U)
         call assert_shape(S, [m,m], 'S', this_module, 'initialize_LR_state')
         ! optional size argument
         if (present(rkmax)) then
            if (rkmax < rk) then
               call stop_error('rkmax < rk!', this_module, 'initialize_LR_state')
            end if
            self%rk = rk
            rka = rkmax
            if (ifrk) then
               if (rkmax==rk) then
                  call stop_error('rkmax must be larger than rk for rank-adaptive DLRA!', this_module, 'initialize_LR_state')
               end if
               write(msg,'(A,I0,A)') 'Rank-adaptivity enabled. Computation will begin with X%rk = ', self%rk+1, '.'
               call logger%log_information(msg, module=this_module, procedure='initialize_LR_state')
            end if
         else
            self%rk = rk
            if (ifrk) then
               rka = rk + 1
            else
               rka = rk
            end if
         end if

         ! allocate & initialize
         allocate(self%U(rka), source=U(1)); call zero_basis(self%U)
         allocate(self%S(rka,rka)); self%S = 0.0_wp
         write(msg,'(3(A,I0),A)') 'size(X%U) = [ ', rka,' ], X%rk = ', self%rk, ', size(U0) = [ ', m,' ]'
         call logger%log_information(msg, module=this_module, procedure='initialize_LR_state')
         ! copy inputs
         if (self%rk > m) then   ! copy the full IC into self%U
            call copy(self%U(:m), U)
            self%S(:m,:m) = S
            write(msg,'(4X,A,I0,A)') 'Transfer the first ', m, ' columns of U0 to X%U.'
            call logger%log_information(msg, module=this_module, procedure='initialize_LR_state')
         else  ! fill the first self%rk columns of self%U with the first self%rk columns of the IC
            call copy(self%U(:self%rk), U(:self%rk))
            self%S(:self%rk,:self%rk) = S(:self%rk,:self%rk)
            write(msg,'(4X,A,I0,A)') 'Transfer all ', m, ' columns of U0 to X%U.'
            call logger%log_information(msg, module=this_module, procedure='initialize_LR_state')
         end if
         ! top up basis (to rka for rank-adaptivity) with orthonormal columns if needed
         if (rka > m) then
            write(msg,'(4X,A,I0,A)') 'Fill remaining ', rka-m, ' columns with orthonormal noise orthonormal to U0.'
            call logger%log_information(msg, module=this_module, procedure='initialize_LR_state')
            allocate(Utmp(rka-m), source=U(1))
            do i = 1, rka-m
               call Utmp(i)%rand()
            end do
            allocate(R(rka-m,rka-m)); R = 0.0_wp
            call orthogonalize_against_basis(Utmp, self%U, info)
            call check_info(info, 'orthogonalize_against_basis', module=this_module, procedure='initialize_LR_state')
            call qr(Utmp, R, info)
            call check_info(info, 'qr', module=this_module, procedure='initialize_LR_state')
            call copy(self%U(m+1:), Utmp)
         end if
      end select
      return
   end subroutine initialize_LR_state

end module laplacian2D_LTI_Lyapunov_Base
