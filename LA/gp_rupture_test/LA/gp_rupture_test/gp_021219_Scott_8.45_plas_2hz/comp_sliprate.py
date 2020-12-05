#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 17 2018

@author: Zhifeng Hu <zhh076@ucsd.edu>
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
dx = 0.1
#ypos = 1173
#xoff = 20
theta_rot = 35


f = open(glob.glob('./*.srf')[0],'r')
f_sliprate = open("sliprate.bin", 'wb')
f_fault = open("subfaults.idx", 'w')
f.readline()
f.readline()
token = f.readline()
nx = int(token.split()[2])
nz = int(token.split()[3])
npt = nx * nz
final_sr = np.zeros((2 * nx, 2 * nz, nt_ref), dtype='float32') # (nx, nz, nt) --> (nt, nz, nx) to save
fault_loc = np.array(np.loadtxt("fault_loc.idx"))
x1 = int(fault_loc[0,0])
x2 = int(fault_loc[1,0])
y1 = int(fault_loc[0,1])
y2 = int(fault_loc[1,1])
x_tmp = np.linspace(x1, x2, np.abs(x2 - x1) + 1) 
y_tmp = [np.float((y2-y1)/(x2-x1))*(x-x1) + y1 for x in x_tmp]
f_interp=interp1d(x_tmp, y_tmp, fill_value='extrapolate')
if x1 < x2:
    new_x = np.arange(x1, x1 + nx * 2 )
    new_y = [np.int(i) for i in f_interp(new_x)]
else:
    new_x = np.arange(x1 + 1 - nx * 2, x1 + 1)
    new_y = [np.int(i) for i in f_interp(new_x)]
    new_x = new_x[::-1]
    new_y = new_y[::-1]

para = np.fromfile('/u/sciteam/hu2/scratch/LA/data_for_Zhifeng/mesh/awp.USC.media', dtype='f').reshape(400, 4200, 6320, 3)
vs = para[:,:,:,1]
rho = para[:,:,:,2]
del para
mu = np.ones((2 * nz, 2 * nx))
for j in range(2 * nz):
        for i in range(2 * nx):
            mu[j, i] = rho[j, new_y[i], new_x[i]] * vs[j, new_y[i], new_x[i]] ** 2
print(f"At [20,20], x = {new_x[20]}, y = {new_y[20]}\n", 
      f"vs = {vs[20, new_y[20], new_x[20]]}",
      f"rho = {rho[20, new_y[20], new_x[20]]}")
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
        # trise = len(sliprate1) * dt
        t1 = int(np.floor(tinit/dt))
        sr[t1 : t1+len(sliprate1)] = sliprate1.copy()
        f_fault.write("{0} {1} {2} {3} {4} {5} {6} {7}\n".format(new_x[2 * i], new_y[2 * i], 2 * j + 1, area / 4, mu[2 * j, 2 * i] , stk, dip, rake))
        f_fault.write("{0} {1} {2} {3} {4} {5} {6} {7}\n".format(new_x[2 * i + 1], new_y[2 * i + 1], 2 * j + 1, area / 4, mu[2 * j, 2 * i + 1], stk, dip, rake))
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
        f_fault.write("{0} {1} {2} {3} {4} {5} {6} {7}\n".format(new_x[2 * i], new_y[2 * i], 2 * j + 2, area / 4, mu[2 * j + 1, 2 * i], stk, dip, rake))
        f_fault.write("{0} {1} {2} {3} {4} {5} {6} {7}\n".format(new_x[2 * i + 1], new_y[2 * i + 1], 2 * j + 2, area / 4, mu[2 * j + 1, 2 * i + 1], stk, dip, rake))
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
