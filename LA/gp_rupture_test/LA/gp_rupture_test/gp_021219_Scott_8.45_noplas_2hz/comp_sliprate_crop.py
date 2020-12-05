#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 17 2018

@author: Zhifeng Hu <zhh076@ucsd.edu>
Input:  *.srf, srf source file from CyberShake
        awp.USC.media, mesh file, including vp, vs, rho
        fault_idx.bin indices of subfaults in the mesh

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

def upsample(s1, dt1, dt2, nt2):
    '''
    For upsample only, since no anti-aliasing is implemented.
    '''
    s2 = np.zeros((nt2,))
    for n in range(nt2):
        t_n = dt2 * n
        n_l = np.floor(t_n / dt1)
        n_r = np.ceil(t_n / dt1)
        if n_l==n_r:
            s2[n] = s1[n_l] 
        elif n_r > len(s1) - 1:
            s2[n] = s1[n_l]
        else:
            t_l = n_l * dt1
            t_r = n_r * dt1
            s2[n] = s1[n_l] + (s1[n_r] - s1[n_l]) * (t_n - t_l) / (t_r - t_l)
    return s2

nt_ref = 4000
dx = 0.2
theta_rot = 35
NX, NY, NZ = 6048, 4200, 400

f = open(glob.glob('./*.srf')[0],'r')
f.readline()
f.readline()
token = f.readline()
nx = int(token.split()[2])
nz = int(token.split()[3])
npt = nx * nz
final_sr = np.zeros((2 * nx, 2 * nz, nt_ref), dtype='float32') # (nx, nz, nt) --> (nt, nz, nx) to save
fault_idx = np.fromfile('fault_idx.bin', dtype='int32').reshape((2 * nz, 2 * nx, 2))
fault_end = np.array(np.loadtxt("fault_loc.idx"))
x1 = int(fault_end[0,0])
x2 = int(fault_end[1,0])

dtop = int(float(f.readline().split()[2]) / dx)
print(f.readline())

f_sliprate = open("sliprate.bin", 'wb')
f_fault = open("subfaults.idx", 'w')

para = np.fromfile('../../data_for_Zhifeng/mesh/awp.USC.crop.media', dtype='f').reshape(NZ, NY, NX, 3)
vs = para[:,:,:,1]
rho = para[:,:,:,2]
del para
mu = np.ones((2 * nz, 2 * nx))
for j in range(2 * nz):
        for i in range(2 * nx):
            mu[j, i] = rho[j, fault_idx[j, i, 1], fault_idx[j, i, 0]] * \
                       vs[j, fault_idx[j, i, 1], fault_idx[j, i, 0]] ** 2

print(f"At [20,20], x = {fault_idx[20, 20, 0]}, y = {fault_idx[20, 20, 1]}\n", 
      f"vs = {vs[20, fault_idx[20, 20, 1], fault_idx[20, 20, 0]]}",
      f"rho = {rho[20, fault_idx[20, 20, 1], fault_idx[20, 20, 0]]}")
del vs
del rho


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
        # trise = len(sliprate1) * dt
        t1 = int(np.floor(tinit/dt))
        sr[t1 : t1+len(sliprate1)] = sliprate1.copy()
        f_fault.write("{0} {1} {2} {3} {4} {5} {6} {7}\n".format(fault_idx[2 * j, 2 * i, 0], fault_idx[2 * j, 2 * i, 1], 2 * j + 1 + dtop, area / 4, mu[2 * j, 2 * i] , stk, dip, rake))
        f_fault.write("{0} {1} {2} {3} {4} {5} {6} {7}\n".format(fault_idx[2 * j, 2 * i + 1, 0], fault_idx[2 * j, 2 * i + 1, 1], 2 * j + 1 + dtop, area / 4, mu[2 * j, 2 * i + 1] , stk, dip, rake))
        final_sr[2 * i, 2 * j, :] = sr
        final_sr[2 * i + 1, 2 * j, :] = sr
    f.seek(last_pos)
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
        # trise = len(sliprate1) * dt
        t1 = int(np.floor(tinit/dt))
        sr[t1 : t1+len(sliprate1)] = sliprate1.copy()
        f_fault.write("{0} {1} {2} {3} {4} {5} {6} {7}\n".format(fault_idx[2 * j + 1, 2 * i, 0], fault_idx[2 * j + 1, 2 * i, 1], 2 * j + 2 + dtop, area / 4, mu[2 * j + 1, 2 * i] , stk, dip, rake))
        f_fault.write("{0} {1} {2} {3} {4} {5} {6} {7}\n".format(fault_idx[2 * j + 1, 2 * i + 1, 0], fault_idx[2 * j + 1, 2 * i + 1, 1], 2 * j + 2 + dtop, area / 4, mu[2 * j + 1, 2 * i + 1] , stk, dip, rake))
        final_sr[2 * i, 2 * j + 1, :] = sr
        final_sr[2 * i + 1, 2 * j + 1, :] = sr

# reverse the fault horizontally when it starts from the right portion.
if x1 > x2:
    final_sr = final_sr[::-1, :, :]
final_sr = np.float32(final_sr.T)
final_sr.tofile(f_sliprate, format='float32')
f_fault.close()
f_sliprate.close()

if x1 > x2:
    print ("Reversed!")
    f_in = open('subfaults.idx', 'r')
    data = f_in.readlines()
    f_in.close()
    time.sleep(10)
    f_out = open('subfaults.idx', 'w')
    for j in range(2 * nz):
        f_out.writelines(reversed(data[2 * nx * j : 2 * nx * (j + 1)]))
    f_out.close()
