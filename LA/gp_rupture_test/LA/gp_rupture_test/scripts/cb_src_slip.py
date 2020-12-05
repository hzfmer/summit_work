#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import os
import re
from read_params import read_params 
import matplotlib.pyplot as plt

args = read_params()

nsrc = args['NSRC']
ntskp = args['NTISKP']
tmax = args['TMAX']
dt = args['DT']
nt_mom = args['NST']
ft_dim = {'6.35': (354, 212), '6.45': (144, 230), '7.35': (980, 240), '8.15': (5116, 220), '8.45': (5116, 220)}
M = re.findall(r"\d+\.\d+", os.getcwd().split('/')[-1])[0]
nx, nz = ft_dim[M]
print(nx, nz)
slipr = np.fromfile('sliprate.bin', 'f').reshape((-1, nz, nx))
nt_slip = slipr.shape[0]
fig, ax = plt.subplots(4, 1, figsize=(6, 6))
for i in range(4):
    np.random.seed()
    ix = np.random.randint(nx)
    iz = np.random.randint(nz)
    disp = iz * nx + ix
    print(ix, iz)
    buf_slip = slipr[:, iz, ix]
    max_slip = max(abs(buf_slip))
    tot_slip = np.sum(buf_slip)
    ruptt = np.argwhere(buf_slip > 1e-9)[0]
    riset = np.sum(buf_slip > 1e-9)
    maxt = np.argmax(buf_slip)
    print(f'rupture time = {ruptt}, rise time = {riset}, max time = {maxt}',
          f'ratio = {(maxt - ruptt) / riset}')
    idx = np.arange(np.max((0, ruptt - 20)), np.min((nt_slip, ruptt + riset + 20)))
    ax[i].plot(idx * dt * ntskp, buf_slip[idx], linewidth=2)
fig.savefig(f'sample_sliprate_{M}.png', dpi=300, bbox='tight', pad_inches=0.1)
