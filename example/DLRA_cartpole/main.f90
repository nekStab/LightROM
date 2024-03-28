program demo
   use LightKrylov
   use LightKrylov_expmlib
   use lightKrylov_utils
   use LightROM_RiccatiSolvers
   use LightROM_LyapunovUtils
   use LightROM_utils
   use cartpole
   use cartpole_riccati

   use stdlib_optval, only : optval 
   use stdlib_linalg, only : eye
   use stdlib_math, only : all_close, logspace
   use stdlib_io_npy, only : save_npy
   implicit none

   ! DLRA
   integer, parameter :: rkmax = n
   integer, parameter :: rk_X0 = 1
   logical, parameter :: verb  = .false.
   logical, parameter :: save  = .false.
   character*128      :: oname

   integer  :: nrk, ndt, rk,  torder
   real(wp) :: dt, Tend
   ! vector of dt values
   real(wp), allocatable :: dtv(:)
   ! vector of rank values
   real(wp), allocatable :: rkv(:)

   real(kind=wp)      :: H(2*n,2*n)
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
   ! RKlib
   ! Exponential propagator
   type(rklib_riccati_mat), allocatable :: RK_propagator
   ! state matrix
   type(state_matrix)      :: Xmat_RKlib(2)
   real(wp), allocatable   :: X_RKlib(:,:,:)
   real(wp)                :: X_RKlib_ref(n,n)
   ! misc
   integer                 :: i, j, k, icnt, irep, nrep

   call initialize_cartpole()

   ! construct Hamiltonian
   H = 0.0_wp
   H(  1:n,    1:n  ) = A
   H(n+1:2*n,n+1:2*n) = -transpose(A)
   H(  1:n,  n+1:2*n) = -BRinvBT
   H(n+1:2*n,  1:n  ) = -CTQC

   call dgeev('N', 'V', 2*n, H, 2*n, wr, wi, VR, 2*n, VR, 2*n, work, lwork, info)

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

   !call init_rand(X0)

   ! initialize exponential propagator
   nrep = 10
   Tend = 0.1_wp
   RK_propagator = rklib_riccati_mat(Tend)

   !allocate(X_RKlib(N, N, nrep))
   !call set_state(X_mat_RKlib(1:1), X0)
   !write(*,'(A10,A26,A26,A20)') 'RKlib:','Tend','|| X_RK - X_ref ||_2/N', 'Elapsed time'
   !write(*,*) '         ------------------------------------------------------------------------'
   !do irep = 1, nrep
   !   call system_clock(count=clock_start)     ! Start Timer
   !   ! integrate
   !   call RK_propagator%matvec(X_mat_RKlib(1), X_mat_RKlib(2))
   !   ! recover output
   !   call get_state(X_RKlib(:,:,irep), X_mat_RKlib(2:2))
   !   ! replace input
   !   call set_state(X_mat_RKlib(1:1), X_RKlib(:,:,irep))
   !   call system_clock(count=clock_stop)      ! Stop Timer
   !   write(*,'(I10,F26.4,E26.8,F18.4," s")') irep, irep*Tend, &
   !                  & norm2(X_RKlib(:,:,irep) - Xref)/N, &
   !                  & real(clock_stop-clock_start)/real(clock_rate)
   !end do

end program demo