#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 17 2018

@author: Zhifeng Hu <zhh076@ucsd.edu>
"""

import numpy as np
from numpy import sin, cos, pi, sqrt
from scipy.interpolate import interp1d
from struct import pack
import os
import sys
import glob

f = open(glob.glob('./*.srf')[0],'r')
#f_out = open("fault_loc.txt", 'w')
f.readline()
f.readline()
token = f.readline()
nx = int(token.split()[2])
nz = int(token.split()[3])
npt = nx * nz

f.readline()
f.readline()

lon = np.zeros((nz, nx), dtype='float32')
lat = np.zeros((nz, nx), dtype='float32')
tinit = np.zeros((2 * nz, 2 * nx), dtype='float32')
ll = np.zeros((2 * nz, 2 * nx, 2), dtype='float32')
for i in range(nz):
    for j in range(nx):
        nl1 = f.readline().split()
        nl2 = f.readline().split()
        tinit[2 * i, 2 * j] = np.float32(nl1[6])
        tinit[2 * i, 2 * j + 1] = np.float32(nl1[6])
        tinit[2 * i, 2 * j] = np.float32(nl1[6])
        tinit[2 * i + 1, 2 * j + 1] = np.float32(nl1[6])
        lon[i, j] = np.float32(nl1[0])
        lat[i, j] = np.float32(nl1[1])
        nt1 = int(nl2[2])
        nskip1 = int(np.ceil(nt1 / 6.))
        for l in range(nskip1):
            f.readline()
tinit.tofile('tinit.bin', format='float32')
if lat.min() == 0:
    print("Error, original zeros in lat\n")

x = np.arange(1, 2 * nx, 2)
x_new = np.arange(0, 2 * nx)
for i in range(nz):
    f_interp = interp1d(x, lon[i, :], fill_value='extrapolate')
    ll[2 * i, :, 0] = f_interp(x_new)
    ll[2 * i, 0, 0] = 2 * ll[2 * i, 1, 0] - ll[2 * i, 2, 0]  # avoid zeros at the boundary
    # ll[2 * i, :, 0] = np.interp(x_new, x, lon[i, :])
    ll[2 * i + 1, :, 0] = ll[2 * i, :, 0]
    
    f_interp = interp1d(x, lat[i, :], fill_value='extrapolate')
    ll[2 * i, :, 1] = f_interp(x_new)
    ll[2 * i, 0, 1] = 2 * ll[2 * i, 1, 1] - ll[2 * i, 2, 1]  # avoid zeros at the boundary
    # ll[2 * i, :, 1] = np.interp(x_new, x, lat[i, :])

    ll[2 * i + 1, :, 1] = ll[2 * i, :, 1]
    if ll[2 * i + 1, :, 1].min() < 1:
        print(ll[2 * i + 1, :, 1])
        sys.exit(-1)


ll.tofile('fault_loc.bin', format='float32')
np.savetxt('fault_surf_loc.txt', ll[0, :, :], fmt='%f')
f.close()
