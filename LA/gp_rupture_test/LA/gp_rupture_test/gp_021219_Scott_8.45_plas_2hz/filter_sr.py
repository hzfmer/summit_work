#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import scipy.fftpack
from scipy.signal import butter, sosfiltfilt, resample
from numpy import sin, pi

def filt_B(data_in, fs, lowcut=0, highcut=11.2, order=4, causal=True):
    if highcut>=fs/2:
        return data_in
    sz = data_in.shape
    if len(sz)>1 and sz[0]<sz[1]:
        data_in = data_in.T
    data_out = np.zeros_like(data_in)
    nyq = fs / 2  
    low = lowcut / nyq
    high = highcut / nyq
    if low==0:  # lowpass
        sos = butter(order, high, analog=False, btype='low', output='sos')
    else:   # bandpass
        sos = butter(order, [low, high], analog=False, btype='band', output='sos')
    if len(sz)>1:
        for i in range(sz[-1]):
            if causal == True:
                data_out[:,i] = sosfilt(sos, data_in[:,i])
                data_out[i, :] = sosfilt(sos, data_out[i, :])
            else:
                data_out[i, :] = sosfiltfilt(sos, data_in[i, :])
    else:
        data_out = sosfiltfilt(sos, data_in)
    return data_out



dt = 0.05
nt = 4000
fs = 1 / dt
nx = 5116
nz = 220

sr = np.fromfile('sliprate.bin', dtype='float32').reshape(nt, nz, nx).T
print(sr.shape)

for j in range(nz):
    for i in range(nx):
        sr[i, j, :] = filt_B(filt_B(sr[i, j, :], fs), fs)

sr.T.tofile('sliprate_1hz.bin', format='%f')
            
