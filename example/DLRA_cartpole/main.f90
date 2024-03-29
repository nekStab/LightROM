program demo
   use LightKrylov
   use LightKrylov_expmlib
   use lightKrylov_utils
   use LightROM_RiccatiSolvers
   use LightROM_LyapunovSolvers
   use LightROM_LyapunovUtils
   use LightROM_utils
   use LightKrylov_expmlib
   use cartpole
   use cartpole_riccati

   use stdlib_optval, only : optval 
   use stdlib_linalg, only : eye
   use stdlib_math, only : all_close, logspace
   use stdlib_io_npy, only : save_npy
   implicit none

   ! DLRA
   integer, parameter :: rkmax = n
   integer, parameter :: rk_X0 = 4
   logical, parameter :: verb  = .false.
   logical, parameter :: save  = .false.
   character*128      :: oname

   integer  :: nrk, ndt, rk,  torder
   real(wp) :: dt, Tend
   ! vector of dt values
   real(wp), allocatable :: dtv(:)
   ! vector of rank values
   integer, allocatable :: rkv(:)

   ! LTI system
   type(lti_system)                :: LTI
   type(cartpole_linop)            :: A
   type(identity_linop)            :: Q
   type(state_vector), allocatable :: CT(:)
   real(kind=wp), allocatable      :: D(:,:)
   integer                         :: p

   ! LR representation
   type(state_vector), allocatable :: U(:)
   type(state_vector), allocatable :: Utest(:)
   real(wp) , allocatable          :: S(:,:)

   real(kind=wp)      :: Hmat(2*n,2*n)
   ! LAPACK
   real(kind=wp)      :: wr(2*n), wi(2*n)
   real(kind=wp)      :: VR(2*n,2*n)
   integer, parameter :: lwork = 1040
   real(kind=wp)      :: work(lwork)
   integer            :: info
   real(kind=wp)      :: UR(2*n,n)
   real(kind=wp)      :: UI(2*n,n)
   logical            :: flag
   real(kind=wp)      :: F(n,n)
   real(kind=wp)      :: Ginv(n,n)
   real(kind=wp)      :: Xref(n,n)

   ! initial condition
   real(wp)                :: U0(n, rkmax)
   real(wp)                :: S0(rkmax,rkmax)
   real(wp)                :: X0(n,n)

   ! OUTPUT
   real(wp)                        :: U_out(N,rkmax)
   real(wp)                        :: X_out(N,N)

   ! SVD
   real(wp)  :: U_svd(N,N)
   real(wp)  :: S_svd(rkmax)
   real(wp)  :: V_svd(rkmax,rkmax)

   ! RKlib
   ! Exponential propagator
   type(rklib_riccati_mat), allocatable :: RK_propagator
   ! state matrix
   type(state_matrix)      :: Xmat_RKlib(2)
   real(wp), allocatable   :: X_RKlib(:,:,:)
   real(wp)                :: X_RKlib_ref(n,n)

   ! misc
   integer                 :: i, j, k, icnt, irep, nrep
   ! timer
   integer   :: clock_rate, clock_start, clock_stop

   call system_clock(count_rate=clock_rate)

   call initialize_cartpole()

   p = 1
   LTI = lti_system()
   allocate(LTI%A,      source=A)
   allocate(LTI%B(1:m), source=B(1:m));
   allocate(LTI%CT(1:p),source=B(1)); call mat_zero(LTI%CT)
   allocate(LTI%D(1:p,1:m)); LTI%D = 0.0_wp

   ! construct Hmat
   Hmat = 0.0_wp
   Hmat(  1:n,    1:n  ) =  Amat
   Hmat(n+1:2*n,n+1:2*n) = -transpose(Amat)
   Hmat(  1:n,  n+1:2*n) = -BRinvBTmat
   Hmat(n+1:2*n,  1:n  ) = -Qmat

   call dgeev('N', 'V', 2*n, Hmat, 2*n, wr, wi, VR, 2*n, VR, 2*n, work, lwork, info)

   ! extract stable eigenspace
   UR = 0.0_wp
   UI = 0.0_wp
   icnt = 0
   flag = .true.
   do i = 1, 2*n
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
   F    = matmul(UR(n+1:2*n,:), transpose(UR(1:n,:))) + matmul(UI(n+1:2*n,:), transpose(UI(1:n,:)))
   Ginv = matmul(UR(  1:n,  :), transpose(UR(1:n,:))) + matmul(UI(  1:n,  :), transpose(UI(1:n,:)))
   call inv(Ginv)
   Xref = matmul(F, Ginv)

   call print_mat(n,n,Xref)

   X0 = matmul(Amat, Xref) + matmul(Xref, transpose(Amat)) + Qmat - matmul(Xref, matmul(BRinvBTmat, Xref))

   call print_mat(n,n,X0)

   ! Define initial condition
   U0 = 0.0_wp
   call random_number(U0(:, 1:rk_X0))
   ! Compute SVD to get low-rank representation
   call svd(U0(:,1:rk_X0), U_svd(:,1:N), S_svd(1:rk_X0), V_svd(1:rk_X0,1:rk_X0))
   S0 = 0.0_wp
   do i = 1,rk_X0
      S0(i,i) = S_svd(i)
   end do
   U0(:,1:rk_X0) = U_svd(:,1:rk_X0)

   ! Compute the full initial condition X0 = U_in @ S0 @ U_in.T
   X0 = matmul( U0(:,1:rk_X0), matmul(S0(1:rk_X0,1:rk_X0), transpose(U0(:,1:rk_X0))))

   ! initialize exponential propagator
   nrep = 10
   Tend = 2.0_wp
   RK_propagator = rklib_riccati_mat(Tend)

   allocate(X_RKlib(N, N, nrep))
   call set_state_mat(Xmat_RKlib(1:1), X0)
   write(*,'(A10,A26,A26,A20)') 'RKlib:','Tend','|| dX/dt ||_2/N', 'Elapsed time'
   write(*,*) '         ------------------------------------------------------------------------'
   do irep = 1, nrep
      call system_clock(count=clock_start)     ! Start Timer
      ! integrate
      call RK_propagator%matvec(Xmat_RKlib(1), Xmat_RKlib(2))
      ! recover output
      call get_state_mat(X_RKlib(:,:,irep), Xmat_RKlib(2:2))
      ! replace input
      call set_state_mat(Xmat_RKlib(1:1), X_RKlib(:,:,irep))
      call system_clock(count=clock_stop)      ! Stop Timer
      write(*,'(I10,F26.4,E26.8,F18.4," s")') irep, irep*Tend, &
                     & norm2(CARE(X_RKlib(:,:,irep), Amat, Qmat, BRinvBTmat))/N, &
                     & real(clock_stop-clock_start)/real(clock_rate)
   end do

   Tend = 10.0_wp 

   ! Choose input ranks and integration steps
   nrk = 4; ndt = 2
   allocate(rkv(1:nrk)); allocate(dtv(1:ndt)); 
   rkv = (/ 1,2,3,4 /)
   dtv = logspace(-3.0_wp, -2.0_wp, ndt)

   do torder = 1, 1 !2
      do i = 1, nrk
         rk = rkv(i)

         allocate(U(1:rk)); call mat_zero(U)
         allocate(S(1:rk,1:rk)); S = 0.0_wp
         write(*,'(A10,I1)') ' torder = ', torder

         do j = ndt, 1, -1
            dt = dtv(j)
            if (verb) write(*,*) '    dt = ', dt, 'Tend = ', Tend

            ! Reset input
            call set_state(U(1:rk), U0(:,1:rk))
            S(1:rk,1:rk) = S0(1:rk,1:rk)

            ! run step
            call system_clock(count=clock_start)     ! Start Timer
            call numerical_low_rank_splitting_integrator_riccati( &
                     & U(1:rk), S(1:rk,1:rk), &
                     & LTI, Q, Rinv, &
                     & Tend, dt, torder, info)
            call system_clock(count=clock_stop)      ! Stop Timer

            ! Reconstruct solution
            call get_state(U_out(:,1:rk), U)
            X_out = matmul(U_out(:,1:rk), matmul(S(1:rk,1:rk), transpose(U_out(:,1:rk))))
      
            write(*,'(A10,I4," TO",I1,F10.6,F8.4,E26.8,F18.4," s")') 'OUTPUT', &
                              & rk, torder, dt, Tend, &
                              & norm2(CARE(X_out, Amat, Qmat, BRinvBTmat))/N, &
                              & real(clock_stop-clock_start)/real(clock_rate)
         end do

         deallocate(U);
         deallocate(S);

      end do
   end do

contains

   !--------------------------------------------------------------------
   !-----     UTILITIES FOR STATE_VECTOR AND STATE MATRIX TYPES    -----
   !--------------------------------------------------------------------

   subroutine get_state_mat(mat_out, state_in)
      !! Utility function to transfer data from a state vector to a real array
      real(kind=wp),          intent(out) :: mat_out(:,:)
      class(abstract_vector), intent(in)  :: state_in(:)
      ! internal variables
      integer :: k, kdim
      mat_out = 0.0_wp
      select type (state_in)
      type is (state_matrix)
         call assert_shape(mat_out, (/ N, N /), 'get_state -> state_matrix', 'mat_out')
         mat_out = reshape(state_in(1)%state, (/ N, N /))
      end select
      return
   end subroutine get_state_mat

   subroutine set_state_mat(state_out, mat_in)
      !! Utility function to transfer data from a real array to a state vector
      class(abstract_vector), intent(out) :: state_out(:)
      real(kind=wp),          intent(in)  :: mat_in(:,:)
      ! internal variables
      integer       :: k, kdim
      select type (state_out)
      type is (state_matrix)
         call assert_shape(mat_in, (/ N, N /), 'set_state -> state_matrix', 'mat_in')
         call mat_zero(state_out)
         state_out(1)%state = reshape(mat_in, shape(state_out(1)%state))
      end select
      return
   end subroutine set_state_mat

end program demo