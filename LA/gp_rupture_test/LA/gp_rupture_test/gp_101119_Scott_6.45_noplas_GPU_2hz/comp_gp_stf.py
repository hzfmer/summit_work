#!/ccs/proj/geo112/hzfmer/summit/opt/anaconda3/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 5 2019

@author: Zhifeng Hu <zhh076@ucsd.edu>
Input:  *.srf, srf source file from CyberShake
        fault_idx.bin indices of subfaults in the mesh

Output: gp_src.bin, gp-type source for IFAULT=5
"""

import numpy as np
from numpy import sin, cos, pi, sqrt
from scipy.signal import resample
import os
import sys
import glob
import struct


nt_ref = 4000
dx = 0.2
theta_rot = 35
NX, NY, NZ = 6320, 4200, 400

f = open(glob.glob('./*.srf')[0],'r')
f.readline()
f.readline()
token = f.readline()
f.readline()
print(f.readline())
nx = int(token.split()[2])
nz = int(token.split()[3])
npt = nx * nz
fault_end = np.array(np.loadtxt("fault_loc.idx"))
x1 = int(fault_end[0,0])
x2 = int(fault_end[1,0])


res = np.zeros((2 * nz, 2 * nx, 12), dtype='float32')
idx = np.genfromtxt('subfaults.idx')
idx = idx[:, : 3].astype('int').reshape(4 * nz *  nx, 3)
stf_type = 3.0

for j in range(nz):
    sys.stdout.write("\rreading subfault %d of %d" % (j+1, nz))
    sys.stdout.flush()
    last_pos = f.tell()
    for i in range(nx):
        nl1 = f.readline().split()
        nl2 = f.readline().split()
        stk = float(nl1[3]) - theta_rot
        dip = float(nl1[4])
        area = float(nl1[5]) / 1.e4 / 4  # reduce the dx from 200 to 100
        tinit= float(nl1[6])
        dt = float(nl1[7])
        rake = float(nl2[0])
        slip = float(nl2[1]) / 1.e2
        nt1 = int(nl2[2]) 
        nskip1 = int(np.ceil(nt1 / 6.))
        for l in range(nskip1):
            tmp = f.readline()
        trise = nt1 * dt
        res[2 * j, 2 * i, :] = [stf_type, slip, tinit, trise, area, 0.0, stk, dip, rake, 0.0, 0.0, 0.0]
        res[2 * j + 1, 2 * i, :] = [stf_type, slip, tinit, trise, area, 0.0, stk, dip, rake, 0.0, 0.0, 0.0]
        res[2 * j, 2 * i + 1, :] = [stf_type, slip, tinit, trise, area, 0.0, stk, dip, rake, 0.0, 0.0, 0.0]
        res[2 * j + 1, 2 * i + 1, :] = [stf_type, slip, tinit, trise, area, 0.0, stk, dip, rake, 0.0, 0.0, 0.0]

# reverse the fault horizontally when it starts from the right portion.
if x1 > x2:
    res = res[:, ::-1, :]
res = res.reshape((-1, 12))

# res.tofile('tp_src_noindex.bin') 
f_out = open('gp_src.bin', 'wb')
fmt = "<3i12f"
for i in range(len(res)):
    f_out.write(struct.pack(fmt, *idx[i, :], *res[i, :])) 
f_out.close()
try:
    os.symlink('gp_src.bin', 'gp_src.bin_0')
except FileExistsError:
    os.remove('gp_src.bin_0')
    os.symlink('gp_src.bin', 'gp_src.bin_0')
    pass
