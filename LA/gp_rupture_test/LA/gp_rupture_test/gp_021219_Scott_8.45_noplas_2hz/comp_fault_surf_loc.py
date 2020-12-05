#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 17 2018

@author: Zhifeng Hu <zhh076@ucsd.edu>
"""

import numpy as np
from numpy import sin, cos, pi, sqrt
from scipy.signal import resample
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

lon = np.zeros((nx, 1))
lat = np.zeros((nx, 1))
for k in range(nx):
    nl1 = f.readline().split()
    nl2 = f.readline().split()
    lon[k] = np.float32(nl1[0])
    lat[k] = np.float32(nl1[1])
    nt1 = int(nl2[2])
    nskip1 = int(np.ceil(nt1 / 6.))
    for l in range(nskip1):
        f.readline()
ll = np.hstack((lon, lat))
np.savetxt('fault_surf_loc.txt', ll, fmt='%f', delimiter=' ', newline='\n')
f.close()
