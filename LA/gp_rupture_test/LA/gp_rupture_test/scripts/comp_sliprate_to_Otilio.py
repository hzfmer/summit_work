#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 17 2018

@author: Zhifeng Hu <zhh076@ucsd.edu>

Input:  *.srf, srf source file from CyberShake
        awp.USC.media, mesh file including vp, vs, rho
        fault_idx.bin, indices of subfaults in the mesh

Output: subfaults.idx, fault information
        sliprate.bin, sliprate time history on each subfault
"""

import numpy as np
from numpy import sin, cos, pi, sqrt
from scipy.signal import resample
from scipy.interpolate import interp1d
from struct import pack
import os
import sys
import glob
import time


f = open(glob.glob('./*.srf')[0],'r')
f_sliprate = open("sliprate.bin", 'wb')
f_fault = open("subfaults.idx", 'w')
f.readline()
f.readline()
token = f.readline()
nx = int(token.split()[2])
nz = int(token.split()[3])
npt = nx * nz
final_sr = np.zeros((2 * nx, 2 * nz, nt_ref), dtype='float32')   # (nx, nz, nt) --> (nt, nz, nx) to save
fault_idx = np.fromfile('fault_idx.bin', dtype='int32').reshape((2 * nz, 2 * nx, 2))

# vs, rho, mu
para = np.fromfile('../../data_for_Zhifeng/mesh/awp.USC.media', dtype='f').reshape(400, 4200, 6320, 3)
vs = para[:,:,:,1]
rho = para[:,:,:,2]
del para
mu = np.ones((2 * nz, 2 * nx))
for j in range(2 * nz):
        for i in range(2 * nx):
            mu[j, i] = rho[j, fault_idx[j, i, 1], fault_idx[j, i, 0]] * \
                       vs[j, fault_idx[j, i, 1], fault_idx[j, i, 0]] ** 2

del vs
del rho


f.readline()
f.readline()

for j in range(nz):
    sys.stdout.write("\rreading subfault %d of %d" % (j+1, nz))
    sys.stdout.flush()
    last_pos = f.tell()
    for i in range(nx):
        sr = np.zeros((nt_ref,))
        nl1 = f.readline().split()
        nl2 = f.readline().split()
        stk = float(nl1[3]) - theta_rot
        dip = float(nl1[4])
        rake = float(nl2[0])
        area = float(nl1[5]) / 1.e4
        tinit= float(nl1[6])
        dt = float(nl1[7])
        nt1 = int(nl2[2])
        nskip1 = int(np.ceil(nt1 / 6.))
        sliprate1 = np.zeros((nt1), dtype=float)
        p1 = 0
        for l in range(nskip1):
            tmp = f.readline()
            nt = len(tmp.split())
            for n in range(nt):
                sliprate1[p1] = float(tmp.split()[n]) / 1.e2
                p1 += 1
        t1 = int(np.floor(tinit/dt))
        sr[t1 : t1+len(sliprate1)] = sliprate1.copy()
        f_fault.write("{0} {1} {2} {3} {4} {5} {6} {7}\n".format(fault_idx[j, i, 0], fault_idx[j, i, 1], j, area, mu[j, i] , stk, dip, rake))
        final_sr[i, j, :] = sr
final_sr = np.float32(final_sr.T)
final_sr.tofile(f_sliprate, format='float32')
f_fault.close()
f_sliprate.close()

