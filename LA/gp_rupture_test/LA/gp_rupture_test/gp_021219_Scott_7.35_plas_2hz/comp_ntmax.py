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

ntmax = 1
dt = 0.05
f = open(glob.glob('./*.srf')[0],'r')
f.readline()
f.readline()
token = f.readline()
nx = int(token.split()[2])
nz = int(token.split()[3])
npt = nx * nz

f.readline()
f.readline()

for k in range(npt):
    nl1 = f.readline().split()
    nl2 = f.readline().split()
    tinit= float(nl1[6])
    nt1 = int(nl2[2])
    nskip1 = int(np.ceil(nt1 / 6.))
    for l in range(nskip1):
        f.readline()
    ntmax = np.max((ntmax, tinit//dt+nt1))
      
print(ntmax)
f.close()
