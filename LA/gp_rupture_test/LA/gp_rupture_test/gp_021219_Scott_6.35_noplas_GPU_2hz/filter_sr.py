#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import scipy.fftpack
from scipy.signal import butter, sosfilt, sosfiltfilt, resample
from numpy import sin, pi
from filter_BU import filt_B

dt = 0.05
nt = 4000
fs = 1 / dt
nx = 354
nz = 212

sr = np.fromfile('sliprate.bin', dtype='float32').reshape(nt, nz, nx).T
print(sr.shape)

for j in range(nz):
    for i in range(nx):
        sr[i, j, :] = filt_B(sr[i, j, :], fs, highcut=2, causal=False)

sr.T.tofile('sliprate_2hz.bin', format='%f')
            
