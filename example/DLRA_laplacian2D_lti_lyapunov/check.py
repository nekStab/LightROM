import sys, os
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm, qr, svdvals

from numpy.linalg import norm as nrm

sys.path.append(os.path.dirname(os.path.abspath('/home/skern/projects/pyGL')))
from solvers.lyapunov_ProjectorSplittingIntegrator import LR_OSI_test, LR_OSI_base
from core.utils import p

# naive Runge-Kutta time integration
def Xdot(t,Xv,A,Q):
    n = A.shape[0]
    X = Xv.reshape((n,n))
    dXdt = A @ X + X @ A.T + Q
    return dXdt.flatten()

def res(X, A, Q):
    return Xdot(0, X.flatten(), A, Q)

run_full     = False
print_ref    = False
print_SVD    = True
print_subSVD = True
print_hlp    = True
print_deb    = True
print_orth   = True
ifp = False
print_from = 0
print_to   = 4
    
Xref = np.load('Xref.npy')
A = np.load('A.npy')
B = np.load('B.npy')
BBT = B @ B.T
N = Xref.shape[0]

if (print_ref):
    print('\nA:')
    for i in range(N):
        print(' '.join(f'{x:7.4f}' for x in A[i,:]))
    print('\nBBT:')
    for i in range(N):
        print(' '.join(f'{x:7.4f}' for x in (BBT)[i,:]))
    print('\nXref:')
    for i in range(N):
        print(' '.join(f'{x:7.4f}' for x in Xref[i,:]))

print('\nDirect solution Xref:\n')
print(f'\t |res| = {nrm(res(Xref, A, BBT))/N:16.12f}')
print(f'\t|Xref| = {nrm(Xref)/N:16.12f}')

U0 = np.load('U0_full.npy')
S0 = np.load('S0_full.npy')
X0 = U0 @ S0 @ U0.T

# RK45
tol  = 1e-12
sol = solve_ivp( Xdot, (0,10), X0.flatten(), args=(A,BBT), dense_output=True, atol=tol, rtol=tol) 
X = sol.sol

if run_full:
    
    Tend = 0.01
    dtv = np.logspace(-2.0, -5.0, 4)
    TO  = [ 1, 2 ]
    rkv = [ 4, 10 ]
    
    print('\nRunge-Kutta solution Xrk:\n')
    tv = np.sort(np.append(np.linspace(0, 1, 11, endpoint=True), Tend))
    for i, t in enumerate(tv):
        print(f'\tt = {t:6.4f}: |res| = {nrm(res(X(t).reshape((N,N)), A, BBT))/N:16.12f}', end='')
        if t==Tend:
            print(' < reference')
        else:
            print('')
    
    print(f'\nLow-rank projector-splitting integrator:\tTend = {Tend:6.4f}\n')
    # LR_OSI
    Xrk = X(Tend).reshape((N,N))
    ifexpm = False
    for torder in TO:
        for rk in rkv:
            print(f'\t rk = {rk}, TO = {torder}')
            for dt in dtv:
                U, S, svals, res_rk, res_rf = LR_OSI_test(A, BBT, X0, Xrk, Xref, Tend, dt, rk, torder=torder, verb=0, ifexpm=ifexpm)
                X = U @ S @ U.T
                print(f' dt = {dt:.0e}: |X| = {nrm(X)/N:.8e} |X-Xrk| = {nrm(X-Xrk)/N:.8e} |X-Xref| = {nrm(X-Xref)/N:.8e} |res| = {np.linalg.norm(res(X, A, BBT))/N:.8e}')
        print('')
    sys.exit()
            
pivot = True

p(U0, 'U0', f='F')
p(S0, 'S0', f='F')

ir = 4
rk = 8

print("\n#   INIT: SVD S0");
ss = svdvals(S0)[:rk].reshape((ir, int(rk/ir)), order='F')
for i in range(ir):
    print('#   ',' '.join(f'{x:19.12f}' for x in ss[i,:]))

Tend = 0.01

#
# RK
#

Xrk = X(Tend).reshape((N,N))
print(f"\n#    Xrk, T = {Tend:6.4f}:");
ss = svdvals(Xrk).reshape((4, int(N/4)), order='F')
for i in range(4):
    print('#   ',' '.join(f'{x:19.12f}' for x in ss[i,:]))

tau = 0.002

U0m = U0.copy()
S0m = S0.copy()

nsteps = round(Tend/tau)

#
# OSI
#

U, S = LR_OSI_base(U0m, S0m, A, BBT, nsteps*tau, tau, torder=1, verb=0)
print("\n#    OSI:");
ss = svdvals(S).reshape((4, int(rk/4)), order='F')
for i in range(4):
    print('#   ',' '.join(f'{x:19.12f}' for x in ss[i,:]))
    
print('\nSolver parameters:')
print(f'\tdt = {tau:8.6f}')
print(f'\tTO = {1}')
print(f'\trk = {rk}\n')
    
t = 0.0

nsteps = 2

for j in range(nsteps):
    
    if j >= print_from and j <= print_to:
        ifp = True
    
    if ifp:
        print(f'\n\n\nSTEP {j+1}:\n\n\n')
   
    if ifp and print_deb: 
        p(U0m, 'U0', f='F')
    
    U1 = expm(tau*A) @ U0m
    
    if ifp and print_deb: 
        p(U1, 'U1', f='F')
    
    if (pivot):
        print(' '.join(f'{U1[:,i].T @ U1[:,i]:16.12f}' for i in range(rk)))
        UA, R, P = qr(U1, mode='economic', pivoting = True)
        print('\n  Pivoting order M QR:', P,'\n')
        print('Rii pre reord:\n',' '.join(f'{x:16.12f}' for x in np.diag(R)))
        R = R[:,np.argsort(P)]
    else:
        UA, R = qr(U1, mode='economic')
    
    if j==1:
        sys.exit()
        
    if ifp and print_deb:
        print('M pre QR: U1.T @ U1\tRii:\n',' '.join(f'{x:16.12f}' for x in np.diag(U1.T @ U1)))
        print('K post QR: Rii:\n',' '.join(f'{x:16.12f}' for x in np.diag(R)))
        
    if ifp and print_orth:
        print(f'Orthogonality: Q1.T @ Q - I = {(UA.T @ UA - np.eye(rk)).max()}')
        p(UA, 'Q', f='F')
        p(R, 'R', f='F')
        print('Rii:\n',' '.join(f'{x:16.12f}' for x in np.diag(R)))
    
    SA    = R @ S0m @ R.T
    
    
    '''
    fname = os.path.join(ffldr,'U_after_M.npy')
    print(f"\nRead U after M from file:\n  {fname}" )
    Uam = np.load(fname)
    print('max(U (after M, fortran) - UA (after M python)) = ')
    for i in range(rk):
        print(f'{(Uam[:,i] - UA[:,i]).max():6.3f}', end=' ')
    fname = os.path.join(ffldr,'S_after_M.npy')
    print(f"\nRead A after M from file:\n  {fname}" )
    Sam = np.load(fname)
    print('max(S (after M, fortran) - SA (after M python)) = ', (Sam - SA).max())
    '''
    if print_SVD:
        print("\nSVD  M step:"); 
        ss = svdvals(SA).reshape((4, int(rk/4)), order='F')
        for i in range(4):
            print(' '.join(f'{x:19.12f}' for x in ss[i,:]))
        print('')
    #
    #
    # GGGGG
    #
    #
      
    '''
    fname = os.path.join(ffldr,'GstepKstep_BBTU.npy')
    print(f"\nRead GstepKstep from file:\n  {fname}" )
    BBTUf = np.load(fname)
    
    print((BBTUf - Qc @ UA).max())
 
    
    sys.exit() 
    '''
    # solve Kdot = Q @ UA with K0 = UA @ SA for one step tau
    K1 = UA @ SA + tau*(BBT @ UA)
    
    if ifp and print_deb: 
        p(UA @ SA, 'K0', f='F')
        p(BBT @ UA, 'Kdot', f='F')
        p(BBT @ UA * tau, 'Kdot*dt', f='F')
        
    if ifp and print_deb: 
        p(K1, 'K1', f='F')
        
    if (pivot):
        U1, Sh, P = qr(K1, pivoting=True, mode='economic')
        print('\n  Pivoting order K QR:', P,'\n')
        print('Rii pre reord:\n',' '.join(f'{x:16.12f}' for x in np.diag(Sh)))
        Sh = Sh[:, np.argsort(P)]
    else:
        U1, Sh = qr(K1, mode='economic')
        
    if ifp and print_orth:
        print(f'Orthogonality: Q1.T @ Q - I = {(UA.T @ UA - np.eye(rk)).max()}')
        p(U1, 'Q', f='F')
        p(Sh, 'Sh = R', f='F')
        print('Rii:\n',' '.join(f'{x:16.12f}' for x in np.diag(Sh)))
        
    if ifp and print_subSVD:
        print("\n\tSVD  G - K step:"); 
        ss = svdvals(Sh).reshape((4, int(rk/4)), order='F')
        for i in range(4):
            print('\t',' '.join(f'{x:19.12f}' for x in ss[i,:]))
        print('')
        
    '''
    # orthonormalise K1
    print('K pre QR: Rii:')
    ip = K1.T @ K1
    print(' '.join(f'{ip[i,i]:16.12f}' for i in range(8)))
    ip = K1.T @ U1
    print('K1.T @ U1:')
    print(' '.join(f'{ip[i,i]:16.12f}' for i in range(8)))
    print('K post QR: Rii:')
    print(' '.join(f'{Sh[i,i]:16.12f}' for i in range(8)))
    '''
    
    # solve Sdot = - U1.T @ Q @ UA with S0 = Sh for one step tau
    if ifp and print_deb: 
        p(U1.T @ BBT @ UA, 'Sdot', f='F')
        p(U1.T @ BBT @ UA * tau, 'Sdot*dt', f='F')
        p(Sh - tau*( U1.T @ BBT @ UA ), 'St', f='F')
        
    St = Sh - tau*( U1.T @ BBT @ UA )
    
    if ifp and print_subSVD:
        print("\n\tSVD  G - S step:"); 
        ss = svdvals(St).reshape((4, int(rk/4)), order='F')
        for i in range(4):
            print('\t',' '.join(f'{x:19.12f}' for x in ss[i,:]))
        print('')
    
    # solve Ldot = U1.T @ Q with L0 = St @ UA.T for one step tau
    if ifp and print_deb: 
        p(St.T, 'St.T', f='F')
        p(UA, 'UA', f='F')
        p((St @ UA.T).T, 'L0.T', f='F')
        p((U1.T @ BBT).T, 'Ldot.T', f='F')
        p((St @ UA.T + tau*( U1.T @ BBT )).T, 'L1.T', f='F')
    
    L1  = St @ UA.T + tau*( U1.T @ BBT )
     
    # update S
    S1  = L1 @ U1
    
    if print_SVD:
        print("\nSVD G step:");
        ss = svdvals(S1).reshape((4, int(rk/4)), order='F')
        for i in range(4):
            print(' '.join(f'{x:19.12f}' for x in ss[i,:]))
    
    U0m = U1.copy()
    S0m = S1.copy()
    
    ifp = False
    t += tau
    
#
#
# DONE
#

print(f"\n#    EXIT nsteps = {nsteps}, T = {t:6.4f}:");
ss = svdvals(S1).reshape((4, int(rk/4)), order='F')
for i in range(4):
    print('#   ',' '.join(f'{x:19.12f}' for x in ss[i,:]))
