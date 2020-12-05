#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 17 2018

@author: Zhifeng Hu <zhh076@ucsd.edu>
"""

import numpy as np
from numpy import sin, cos, pi, sqrt
import os
import sys
import glob
import time


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


if not os.path.isfile('fault_full_loc.txt'):
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

    mx = 6320
    my = 4200
    ll = np.fromfile('../scripts/surf.grid', dtype='float64', count=2 * my * mx).reshape(my, mx, 2)
    ll_fault = [np.float32((ll[new_y[i], new_x[i], 0], ll[new_y[i], new_x[i], 1])) for i in range(len(new_x))]
    np.savetxt('fault_full_loc.txt', ll, fmt='%f')
    # np.array(ll_fault).tofile('latlon_fault.bin')

