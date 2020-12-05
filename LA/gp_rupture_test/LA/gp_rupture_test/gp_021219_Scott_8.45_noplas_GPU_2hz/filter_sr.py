#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from filter_BU import filt_B
from read_params import read_params


args = read_params()
nt = args['NST'] / args['NTISKP']
dt = args['DT'] * args['NTISKP']
fs = 1 / dt
ft_dim = {'6.35': (354, 212), '7.35': (980, 240), '8.45': (5116, 220)}
M = re.findall(r"\d+\.\d+", os.getcwd().split('/')[-1])[0]
nx, nz = ft_dim[M]

sr = np.fromfile('sliprate.bin', dtype='float32').reshape(nt, nz, nx).T
print(sr.shape)

for j in range(nz):
    for i in range(nx):
        sr[i, j, :] = filt_B(sr[i, j, :], fs, highcut=2, causal=False)

sr.T.tofile('sliprate_2hz.bin', format='%f')
            
