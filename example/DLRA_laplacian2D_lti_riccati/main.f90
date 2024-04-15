program demo
   use LightKrylov
   use LightKrylov_expmlib
   use LightKrylov_utils

   use LightROM_AbstractLTIsystems
   use LightROM_utils
   
   use LightROM_LyapunovSolvers
   use LightROM_LyapunovUtils
   use LightROM_RiccatiSolvers

   use Laplacian2D_LTI_Riccati_Base
   use Laplacian2D_LTI_Riccati_Operators
   use Laplacian2D_LTI_Riccati_RKlib

   use stdlib_optval, only : optval 
   use stdlib_linalg, only : eye
   use stdlib_math, only : all_close, logspace
   use stdlib_io_npy, only : save_npy
   implicit none

   ! DLRA
   integer, parameter :: rkmax = 14
   integer, parameter :: rk_X0 = 10
   logical, parameter :: verb  = .false.
   logical, parameter :: save  = .false.
   character*128      :: oname
   ! rk_B is set in laplacian2D_lti_base.f90

   integer  :: nrk, ndt, rk,  torder
   real(wp) :: dt, Tend
   ! vector of dt values
   real(wp), allocatable :: dtv(:)
   ! vector of rank values
   integer,  allocatable :: rkv(:)

   ! Exponential propagator (RKlib).
   type(rklib_riccati_mat), allocatable  :: RK_prop_ricc

   ! LTI system
   type(lti_system)                :: LTI

   ! Laplacian
   type(laplace_operator),   allocatable :: A

   ! LR representation
   type(LR_state)                  :: X
   type(state_vector), allocatable :: U(:)
   real(wp) , allocatable          :: S(:,:)
   
   !> STATE MATRIX (RKlib)
   type(state_matrix)              :: X_mat_RKlib(2)
   real(wp), allocatable           :: X_RKlib(:,:,:)
   real(wp)                        :: X_RKlib_ref(N,N)

    ! Initial condition
   real(wp)                        :: U0(N, rkmax)
   real(wp)                        :: S0(rkmax,rkmax)
   real(wp)                        :: X0(N,N)

   ! OUTPUT
   real(wp)                        :: U_out(N,rkmax)
   real(wp)                        :: X_out(N,N)

   !> Information flag.
   integer                         :: info
   integer                         :: i, j, irep, nrep

   ! PROBLEM DEFINITION
   real(wp)  :: Adata(N,N)

   ! SVD
   real(wp)  :: U_svd(N,N)
   real(wp)  :: S_svd(rkmax)
   real(wp)  :: V_svd(rkmax,rkmax)

   ! LAPACK SOLUTION RICATTI
   real(kind=wp)      :: Hdata(2*N,2*N)
   real(kind=wp)      :: wr(2*N), wi(2*N)
   real(kind=wp)      :: VR(2*N,2*N)
   integer, parameter :: lwork = 1040
   real(kind=wp)      :: work(lwork)
   real(kind=wp)      :: UR(2*N,N)
   real(kind=wp)      :: UI(2*N,N)
   logical            :: flag
   real(kind=wp)      :: F(N,N)
   real(kind=wp)      :: Ginv(N,N)
   real(kind=wp)      :: Xref(N,N)
   integer :: icnt

   ! timer
   integer   :: clock_rate, clock_start, clock_stop
   
   call system_clock(count_rate=clock_rate)

   ! Set up problem
   call initialize_problem(1.0_wp, 0.0_wp)

   ! Define LTI system
   LTI = lti_system()
   allocate(LTI%A,          source=A)
   allocate(LTI%B(1:rk_b),  source=B(1:rk_b))
   allocate(LTI%CT(1:rk_c), source=CT(1:rk_c))
   allocate(LTI%D(1:rk_c,1:rk_b)); LTI%D = 0.0_wp

   write(*,*)
   write(*,*) 'RICCATI EQUATION FOR THE 2D LAPLACE OPERATOR:'
   write(*,*)
   write(*,*) '    Algebraic Lyapunov equation:'
   write(*,*) '                0 = A.T @ X + X @ A + C.T @ Qc @ C'
   write(*,*)
   write(*,*) '    Differential Lyapunov equation:'
   write(*,*) '          \dot{X} = A.T @ X + X @ A + C.T @ Qc @ C'
   write(*,*)
   write(*,'(A16,I4,"x",I4)') '  Problem size: ', N, N
   write(*,*)
   write(*,*) '  Initial condition: rank(X0)            =', rk_X0
   write(*,*) '  Inhomogeneity:     rank(C.T @ Qc @ C)  =', rk_C
   write(*,*)
   write(*,*) '---------------------------------------------'
   write(*,*)

   !------------------
   ! COMPUTE EXACT SOLUTION OF THE RICCATI EQUATION WITH LAPACK
   !------------------

   write(*,*) 'I.   Exact solution of the algebraic Riccati equation (LAPACK):'
   call system_clock(count=clock_start)     ! Start Timer

   ! Explicit 2D laplacian
   call build_operator(Adata)

   ! construct Hmat
   Hdata = 0.0_wp
   Hdata(  1:N,    1:N  ) =  Adata
   Hdata(N+1:2*N,N+1:2*N) = -transpose(Adata)
   Hdata(  1:N,  N+1:2*N) = -BRinvBTdata
   Hdata(N+1:2*N,  1:N  ) = -CTQcCdata

   call dgeev('N', 'V', 2*N, Hdata, 2*N, wr, wi, VR, 2*N, VR, 2*N, work, lwork, info)

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
   call inv(Ginv)
   Xref = matmul(F, Ginv)

   ! sanity check
   !call print_mat(N,N,Xref)
   !X0 = CARE(Xref, Adata, CTQcCdata, BRinvBTdata)

   call system_clock(count=clock_stop)      ! Stop Timer
   write(*,'(A40,F10.4," s")') '--> X_ref.    Elapsed time:', real(clock_stop-clock_start)/real(clock_rate)
   write(*,*)

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

   write(*,*)
   write(*,*) 'II.  Compute approximate solution of the differential Riccati equation using RKlib:'
   write(*,*)

   ! initialize exponential propagator
   nrep = 10
   Tend = 0.1_wp
   RK_prop_ricc = rklib_riccati_mat(Tend)

   allocate(X_RKlib(N, N, nrep))
   call set_state(X_mat_RKlib(1:1), X0)
   write(*,'(A10,A26,A26,A20)') 'RKlib:','Tend','|| dX/dt ||_2/N', 'Elapsed time'
   write(*,*) '         ------------------------------------------------------------------------'
   do irep = 1, nrep
      call system_clock(count=clock_start)     ! Start Timer
      ! integrate
      call RK_prop_ricc%matvec(X_mat_RKlib(1), X_mat_RKlib(2))
      ! recover output
      call get_state(X_RKlib(:,:,irep), X_mat_RKlib(2:2))
      ! replace input
      call set_state(X_mat_RKlib(1:1), X_RKlib(:,:,irep))
      call system_clock(count=clock_stop)      ! Stop Timer
      write(*,'(I10,F26.4,E26.8,F18.4," s")') irep, irep*Tend, &
                     & norm2(CARE(X_RKlib(:,:,irep), Adata, CTQcCdata, BRinvBTdata))/N, &
                     & real(clock_stop-clock_start)/real(clock_rate)
   end do

   !------------------
   ! COMPUTE DLRA FOR SHORTEST INTEGRATION TIMES FOR DIFFERENT DT AND COMPARE WITH RK SOLUTION
   !------------------

   write(*,*)
   write(*,*) 'III. Compute approximate solution of the differential Riccati equation using DLRA:'
   write(*,*)
   write(*,'(A10,A4,A4,A10,A8,A26,A20)') 'DLRA:','  rk',' TO','dt','Tend','|| X_DLRA - X_RK ||_2/N', 'Elapsed time'
   write(*,*) '         ------------------------------------------------------------------------'

   ! Choose relevant reference case from RKlib
   X_RKlib_ref = X_RKlib(:,:,1)

   ! Choose input ranks and integration steps
   nrk = 4; ndt = 5
   allocate(rkv(1:nrk)); allocate(dtv(1:ndt)); 
   rkv = (/ 2, 6, 10, 14 /)
   dtv = logspace(-5.0_wp, -1.0_wp, ndt)

   X = LR_state()

   do torder = 1, 2
      do i = 1, nrk
         rk = rkv(i)

         allocate(U(1:rk)); call mat_zero(U)
         allocate(S(1:rk,1:rk)); S = 0.0_wp
         allocate(X%U(1:rk), source=U(1:rk))
         allocate(X%S(1:rk,1:rk)); 
         write(*,'(A10,I1)') ' torder = ', torder

         do j = ndt, 1, -1
            dt = dtv(j)
            if (verb) write(*,*) '    dt = ', dt, 'Tend = ', Tend

            ! Reset input
            call X%set_LR_state(U0(:,1:rk), S0(1:rk,1:rk))

            ! run step
            call system_clock(count=clock_start)     ! Start Timer
            call numerical_low_rank_splitting_integrator_riccati(X, LTI, Qc, Rinv, Tend, dt, torder, info)
            call system_clock(count=clock_stop)      ! Stop Timer

            ! Reconstruct solution
            call get_state(U_out(:,1:rk), X%U)
            X_out = matmul(U_out(:,1:rk), matmul(X%S, transpose(U_out(:,1:rk))))
      
            write(*,'(A10,I4," TO",I1,F10.6,F8.4,E26.8,F18.4," s")') 'OUTPUT', &
                              & rk, torder, dt, Tend, &
                              & norm2(X_RKlib_ref - X_out)/N, &
                              & real(clock_stop-clock_start)/real(clock_rate)
         end do

         if (save) then
            write(oname,'("example/DLRA_laplacian2D_riccati/data_X_DRLA_TO",I1,"_rk",I2.2,".npy")') torder, rk
            call save_npy(oname, X_out)
         end if

         deallocate(X%U);
         deallocate(X%S);
         deallocate(U);
         deallocate(S);

      end do
   end do
   deallocate(rkv); deallocate(dtv);

   ! Set up problem with non-zero R
   call initialize_problem(1.0_wp, 1.0_wp)

   ! Reset LTI system
   call mat_copy(LTI%B, B)
   call mat_copy(LTI%CT, CT)

   write(*,*)
   write(*,*) 'RICCATI EQUATION FOR THE 2D LAPLACE OPERATOR:'
   write(*,*)
   write(*,*) '    Algebraic Lyapunov equation:'
   write(*,*) '                0 = A.T @ X + X @ A + C.T @ Qc @ C - X @ B @ R^{-1} @ B.T @ X'
   write(*,*)
   write(*,*) '    Differential Lyapunov equation:'
   write(*,*) '          \dot{X} = A.T @ X + X @ A + C.T @ Qc @ C - X @ B @ R^{-1} @ B.T @ X'
   write(*,*)
   write(*,'(A16,I4,"x",I4)') '  Problem size: ', N, N
   write(*,*)
   write(*,*) '  Initial condition: rank(X0)            =', rk_X0
   write(*,*) '  Inhomogeneity:     rank(C.T @ Qc @ C)  =', rk_C
   write(*,*) '  Inhomogeneity:     rank(B)             =', rk_B
   write(*,*)
   write(*,*) '---------------------------------------------'
   write(*,*)

   !------------------
   ! COMPUTE EXACT SOLUTION OF THE RICCATI EQUATION WITH LAPACK
   !------------------

   write(*,*) 'I.   Exact solution of the algebraic Riccati equation (LAPACK):'
   call system_clock(count=clock_start)     ! Start Timer

   ! construct Hmat
   Hdata = 0.0_wp
   Hdata(  1:N,    1:N  ) =  Adata
   Hdata(N+1:2*N,N+1:2*N) = -transpose(Adata)
   Hdata(  1:N,  N+1:2*N) = -BRinvBTdata
   Hdata(N+1:2*N,  1:N  ) = -CTQcCdata

   call dgeev('N', 'V', 2*N, Hdata, 2*N, wr, wi, VR, 2*N, VR, 2*N, work, lwork, info)

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
   call inv(Ginv)
   Xref = matmul(F, Ginv)

   ! sanity check
   !call print_mat(N,N,Xref)
   !X0 = CARE(Xref, Adata, CTQcCdata, BRinvBTdata)

   call system_clock(count=clock_stop)      ! Stop Timer
   write(*,'(A40,F10.4," s")') '--> X_ref.    Elapsed time:', real(clock_stop-clock_start)/real(clock_rate)
   write(*,*)

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

   write(*,*)
   write(*,*) 'II.  Compute approximate solution of the differential Riccati equation using RKlib:'
   write(*,*)

   ! initialize exponential propagator
   nrep = 10
   Tend = 0.1_wp
   RK_prop_ricc = rklib_riccati_mat(Tend)

   X_RKlib = 0.0_wp
   call set_state(X_mat_RKlib(1:1), X0)
   write(*,'(A10,A26,A26,A20)') 'RKlib:','Tend','|| dX/dt ||_2/N', 'Elapsed time'
   write(*,*) '         ------------------------------------------------------------------------'
   do irep = 1, nrep
      call system_clock(count=clock_start)     ! Start Timer
      ! integrate
      call RK_prop_ricc%matvec(X_mat_RKlib(1), X_mat_RKlib(2))
      ! recover output
      call get_state(X_RKlib(:,:,irep), X_mat_RKlib(2:2))
      ! replace input
      call set_state(X_mat_RKlib(1:1), X_RKlib(:,:,irep))
      call system_clock(count=clock_stop)      ! Stop Timer
      write(*,'(I10,F26.4,E26.8,F18.4," s")') irep, irep*Tend, &
                     & norm2(CARE(X_RKlib(:,:,irep), Adata, CTQcCdata, BRinvBTdata))/N, &
                     & real(clock_stop-clock_start)/real(clock_rate)
   end do

   write(*,*)
   write(*,*) 'III. Compute approximate solution of the differential Riccati equation using DLRA:'
   write(*,*)
   write(*,'(A10,A4,A4,A10,A8,A26,A20)') 'DLRA:','  rk',' TO','dt','Tend','|| X_DLRA - X_RK ||_2/N', 'Elapsed time'
   write(*,*) '         ------------------------------------------------------------------------'

   ! Choose relevant reference case from RKlib
   X_RKlib_ref = X_RKlib(:,:,1)

   ! Choose input ranks and integration steps
   nrk = 4; ndt = 5
   allocate(rkv(1:nrk)); allocate(dtv(1:ndt));
   rkv = (/ 2, 6, 10, 14 /)
   dtv = logspace(-5.0_wp, -1.0_wp, ndt)

   X = LR_state()

   do torder = 1, 2
      do i = 1, nrk
         rk = rkv(i)

         allocate(U(1:rk)); call mat_zero(U)
         allocate(X%U(1:rk), source=U(1:rk))
         allocate(X%S(1:rk,1:rk)); 
         write(*,'(A10,I1)') ' torder = ', torder

         do j = ndt, 1, -1
            dt = dtv(j)
            if (verb) write(*,*) '    dt = ', dt, 'Tend = ', Tend

            ! Reset input
            call X%set_LR_state(U0(:,1:rk), S0(1:rk,1:rk))

            ! run step
            call system_clock(count=clock_start)     ! Start Timer
            call numerical_low_rank_splitting_integrator_riccati(X, LTI, Qc, Rinv, Tend, dt, torder, info)
            call system_clock(count=clock_stop)      ! Stop Timer

            ! Reconstruct solution
            call get_state(U_out(:,1:rk), X%U)
            X_out = matmul(U_out(:,1:rk), matmul(X%S, transpose(U_out(:,1:rk))))
      
            write(*,'(A10,I4," TO",I1,F10.6,F8.4,E26.8,F18.4," s")') 'OUTPUT', &
                              & rk, torder, dt, Tend, &
                              & norm2(X_RKlib_ref - X_out)/N, &
                              & real(clock_stop-clock_start)/real(clock_rate)
         end do

         if (save) then
            write(oname,'("example/DLRA_laplacian2D_riccati/data_X_DRLA_TO",I1,"_rk",I2.2,".npy")') torder, rk
            call save_npy(oname, X_out)
         end if

         deallocate(X%U);
         deallocate(X%S);
         deallocate(U);

      end do
   end do

end program demo