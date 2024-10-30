import sys, os, glob, math, re
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm, qr, svdvals

from numpy.linalg import norm as nrm

sys.path.append(os.path.dirname(os.path.abspath('/home/skern/projects/pyGL')))
from solvers.lyapunov_ProjectorSplittingIntegrator import LR_OSI_test, LR_OSI_base
from core.utils import p

def fsc(number):
    if number == 0:
        return "0.00E+00"
    
    exponent = int(math.floor(math.log10(abs(number))))
    mantissa = number / (10 ** exponent)
    
    # Shift decimal point one place left and adjust exponent
    mantissa /= 10
    exponent += 1
    
    # Format the result
    return f"{mantissa:.2f}E{exponent:+03d}"

fldr = 'local'
base = 'DATA_IVa'

N = 16

rk = 7
TO = 2
dt = 0.0001
tolv = [ 1e-4, 1e-6, 1e-8 ]

for tol in tolv:
    pattern = f'{base}_rk{rk:02d}_TO{TO}_dt{dt:8.6f}_inctol{fsc(tol)}*.npy'
    files = glob.glob(os.path.join(fldr, pattern))
    nfiles = len(files)
    ncase  = int(nfiles/2.0)
    print(pattern)
    Ud = np.empty((ncase, N,rk))
    Sd = np.empty((ncase, rk,rk))
    tv    = np.empty((ncase,))
    for i in range(ncase):
        print(f'\t{i+1}')
        pattern2 = f'{base}_rk{rk:02d}_TO{TO}_dt{dt:8.6f}_inctol{fsc(tol)}_{i+1:03d}_s*.npy'
        fs = glob.glob(os.path.join(fldr, pattern2))
        if len(fs) == 2:
            m = re.search('_s[0-9]+', fs[0])
            step = int(re.search(r'\d+', m.group())[0])
            tv = np.append(tv, step*dt)
            for ifile in fs:
                if ifile.endswith('_S.npy'):
                    Sd[i,:,:] = np.load(ifile)
                else:
                    Ud[i,:,:] = np.load(ifile)
    print('\tCross-check:')
    for i, (U,S) in enumerate(zip(Ud,Sd)):
        p(U.T @ U, 'inner_prod')
        p(S, 'S')
        print(np.allclose(S, S.T))
        print(svdvals(S))
            
            
    

