#!/ccs/proj/geo112/hzfmer/summit/opt/anaconda3/bin/python
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 17 2018

@author: Zhifeng Hu <zhh076@ucsd.edu>
"""

import numpy as np
from numpy import sin, cos, pi, sqrt
from scipy.signal import resample
from scipy.interpolate import interp1d
from mpi4py import MPI
import os
import sys
import glob
import time


def upsample(s1, dt1, dt2, nt2):
    '''
    For upsample only, since no anti-aliasing is implemented.
    '''
    s2 = np.zeros((nt2,))
    for n in range(nt2):
        t_n = dt2 * n
        n_l = int(np.floor(t_n / dt1))
        n_r = int(np.ceil(t_n / dt1))
        if n_l==n_r:
            s2[n] = s1[n_l] 
        elif n_r > len(s1) - 1:
            s2[n] = s1[n_l]
        else:
            t_l = n_l * dt1
            t_r = n_r * dt1
            s2[n] = s1[n_l] + (s1[n_r] - s1[n_l]) * (t_n - t_l) / (t_r - t_l)
    return s2

nt_ref = 2000
nt_des = 10 * nt_ref
theta_rot = 35

f = open(glob.glob('./*.srf')[0],'r')
f.readline()
f.readline()
token = f.readline()
nx = int(token.split()[2])
nz = int(token.split()[3])
f.close()

comm = MPI.COMM_WORLD
ncpus = comm.Get_size()
rank = comm.Get_rank()

if not os.path.isfile('sliprate2.bin'):
    final_sr = np.fromfile('sliprate.bin', dtype='float32').reshape(nt_ref, 2 * nz, 2 * nx)
    final_sr2 = np.zeros((2 * nx, 2 * nz * nt_des), dtype='float32') # (nx, nz, nt_afterinterpolation) --> (nt, nz, nx) to save
    for i in range(rank, 2 * nx, ncpus):
        for j in range(2 * nz):
            final_sr2[i, j * nt_des : (j +1) * nt_des] = upsample(final_sr[:, j, i], 0.05, 0.005, nt_des)
            # final_sr2[i, j * nt_des : (j +1) * nt_des] = np.zeros((nt_des, ))
        if rank == 0:
            print(f'\r Upsampling column {i} from {final_sr2.shape}', end=' ', flush=True)
            for p in range(1, ncpus):
                comm.Recv(final_sr2[i + p, :], source=p, tag=100)
        else:
            comm.Send(final_sr2[i, :], dest=0, tag=100)
            # print(f'\rSent data from {rank}', end=' ', flush=True)
        comm.Barrier()
    if rank == 0:
        print("\nTo reshape")
        new_sr = np.reshape(np.ravel(final_sr2), (2 * nx, 2 * nz, nt_des))
        print("\nReshaped")
        new_sr.T.tofile('sliprate2.bin', format='float32')
    comm.Barrier()
MPI.Finalize()

#final_sr = np.fromfile('sliprate.bin', dtype='float32').reshape(nt_ref, 2 * nz, 2 * nx)
final_sr = np.fromfile('sliprate2.bin', dtype='float32').reshape(nt_des, 2 * nz, 2 * nx)
test_sr = final_sr[:, 16, 540]
test_sr.tofile('test_sr.bin')
del final_sr
wn = np.random.normal(0, 1, nt_ref)
res = np.convolve(test_sr, wn, 'same')

