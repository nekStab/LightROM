module laplacian2D_LTI_Lyapunov_Base
   ! Standard Library.
   use stdlib_stats_distribution_normal, only: normal => rvs_normal
   use stdlib_io_npy, only: save_npy
   use stdlib_strings, only: replace_all
   use stdlib_linalg, only: svdvals
   use stdlib_optval, only : optval
   ! LightKrylov for linear algebra.
   use LightKrylov
   use LightKrylov_Logger
   use LightKrylov_Utils, only : assert_shape
   use LightKrylov_AbstractVectors ! zero_basis
   ! LightROM
   use LightROM_AbstractLTIsystems ! LR_state
   implicit none

   private :: this_module
   character(len=*), parameter :: this_module = 'Laplacian2D_LTI_Lyapunov_Base'
   
   ! problem parameters
   public  :: N, nx, dx, dx2, L, rk_b, B, BBT

   !------------------------------
   !-----     PARAMETERS     -----
   !------------------------------

   ! --> Mesh related parameters.
   real(dp),      parameter :: L  = 1.0_dp  !> Domain length
   integer,       parameter :: nx = 4       !> Number of grid points per direction
   integer,       parameter :: N  = nx**2   !> total number of grid points
   real(dp),      parameter :: dx = L/nx    !> Grid size.
   real(dp),      parameter :: dx2= dx**2   !> Grid size.
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
      real(dp) :: state(N) = 0.0_dp
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
      real(dp)      :: state(N**2) = 0.0_dp
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
   real(dp)                 :: BBT(N**2)

contains

   !-----     TYPE-BOUND PROCEDURE FOR VECTORS     -----

   subroutine vector_zero(self)
      class(state_vector), intent(inout) :: self
      self%state = 0.0_dp
   end subroutine vector_zero

   real(dp) function vector_dot(self, vec) result(alpha)
      class(state_vector)       , intent(in) :: self
      class(abstract_vector_rdp), intent(in) :: vec
      select type(vec)
      type is (state_vector)
         alpha = dot_product(self%state, vec%state)
      class default
         call stop_error('vec must be a state_vector', this_module, 'dot')
      end select
   end function vector_dot

   integer function vector_get_size(self) result(ntot)
     class(state_vector), intent(in) :: self
     ntot = nx
   end function vector_get_size

   subroutine vector_scal(self, alpha)
      class(state_vector), intent(inout) :: self
      real(dp)           , intent(in)    :: alpha
      self%state = self%state * alpha
   end subroutine vector_scal

   subroutine vector_axpby(alpha, vec, beta, self)
      class(state_vector)       , intent(inout) :: self
      class(abstract_vector_rdp), intent(in)    :: vec
      real(dp)                  , intent(in)    :: alpha, beta
      select type(vec)
      type is (state_vector)
         self%state = beta*self%state + alpha*vec%state
      class default
         call stop_error('vec must be a state_vector', this_module, 'axpby')
      end select
   end subroutine vector_axpby
   
   subroutine vector_rand(self, ifnorm)
      class(state_vector), intent(inout) :: self
      logical, optional,   intent(in)    :: ifnorm
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
   end subroutine vector_rand

   !-----     TYPE-BOUND PROCEDURE FOR MATRICES     -----

   subroutine matrix_zero(self)
      class(state_matrix), intent(inout) :: self
      self%state = 0.0_dp
   end subroutine matrix_zero

   real(dp) function matrix_dot(self, vec) result(alpha)
      class(state_matrix)       , intent(in) :: self
      class(abstract_vector_rdp), intent(in) :: vec
      select type(vec)
      type is(state_matrix)
         alpha = dot_product(self%state, vec%state)
      class default
         call stop_error('vec must be a state_matrix', this_module, 'matrix_dot')
      end select
   end function matrix_dot

   integer function matrix_get_size(self) result(N)
     class(state_matrix), intent(in) :: self
     N = N
   end function matrix_get_size

   subroutine matrix_scal(self, alpha)
      class(state_matrix), intent(inout) :: self
      real(dp)           , intent(in)    :: alpha
      self%state = self%state * alpha
   end subroutine matrix_scal  

   subroutine matrix_axpby(alpha, vec, beta, self)
      class(state_matrix)       , intent(inout) :: self
      class(abstract_vector_rdp), intent(in)    :: vec
      real(dp)                  , intent(in)    :: alpha, beta
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

   !-----------------------------------------------------------------------
   !-----     TYPE BOUND PROCEDURE FOR SYM LOW RANK REPRESENTATION    -----
   !-----------------------------------------------------------------------

   subroutine initialize_LR_state(self, U, S, rk, rkmax, if_rank_adaptive)
      class(LR_state),            intent(inout) :: self
      class(abstract_vector_rdp), intent(in)    :: U(:)
      real(dp),                   intent(in)    :: S(:,:)
      integer,                    intent(in)    :: rk
      integer, optional,          intent(in)    :: rkmax
      logical, optional,          intent(in)    :: if_rank_adaptive
      logical                                   :: ifrk

      ! internals
      character(len=*), parameter :: this_procedure = 'initialize_LR_state'
      class(abstract_vector_rdp), allocatable   :: Utmp(:)
      integer :: i, m, rka, info, n_rem, m_init
      character(len=128) :: msg

      ifrk = optval(if_rank_adaptive, .false.)

      select type (U)
      type is (state_vector)
         ! set time and optional args
         call self%reset(full=.true.)

         m = size(U)
         call assert_shape(S, [m,m], 'S', this_module, this_procedure)
         ! optional size argument
         if (present(rkmax)) then
            if (rkmax < rk) then
               call stop_error('rkmax < rk!', this_module, this_procedure)
            end if
            self%rk = rk
            rka = rkmax
            if (ifrk) then
               if (rkmax==rk) then
                  call stop_error('rkmax must be larger than rk for rank-adaptive DLRA!', this_module, this_procedure)
               end if
               write(msg,'(A,I0,A)') 'Rank-adaptivity enabled. Computation will begin with X%rk = ', self%rk+1, '.'
               call log_information(msg, this_module, 'initialize_LR_state')
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
         allocate(self%S(rka,rka)); self%S = 0.0_dp
         write(msg,'(3(A,I0),A)') 'size(X%U) = [ ', rka,' ], X%rk = ', self%rk, ', size(U0) = [ ', m,' ]'
         call log_information(msg, this_module, this_procedure)
         ! copy inputs
         if (self%rk > m) then   ! copy the full IC into self%U
            call copy(self%U(:m), U)
            self%S(:m,:m) = S
            write(msg,'(4X,A,I0,A)') 'Transfer all ', m, ' columns of U0 to X%U.'
            call log_information(msg, this_module, this_procedure)
         else  ! fill the first self%rk columns of self%U with the first self%rk columns of the IC
            call copy(self%U(:self%rk), U(:self%rk))
            self%S(:self%rk,:self%rk) = S(:self%rk,:self%rk)
            write(msg,'(4X,A,I0,A)') 'Transfer the first ', self%rk, ' columns of U0 to X%U.'
            call log_information(msg, this_module, this_procedure)
         end if        

         ! top up basis (to rka for rank-adaptivity) with orthonormal columns if needed
         m_init = min(self%rk, m)
         n_rem = rka - m_init
         if (m > 0) then
            write(msg,'(4X,A,I0,A)') 'Fill remaining ', n_rem, ' columns with orthonormal noise orthonormal to X%U.'
            call log_information(msg, this_module, this_procedure)
            allocate(Utmp(n_rem), source=U(1))
            call initialize_random_orthonormal_basis(Utmp)
            call orthogonalize_against_basis(Utmp, self%U(:m_init), info)
            call copy(self%U(m_init+1:), Utmp)
         end if
      class default
         call type_error('U', 'state_vector', 'IN', this_module, this_procedure)
      end select
   end subroutine initialize_LR_state

end module laplacian2D_LTI_Lyapunov_Base