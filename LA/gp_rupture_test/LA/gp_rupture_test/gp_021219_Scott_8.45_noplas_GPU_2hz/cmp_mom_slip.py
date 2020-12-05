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
ft_dim = {'6.35': (354, 212), '7.35': (980, 240), '8.45': (5116, 220)}
M = re.findall(r"\d+\.\d+", os.getcwd().split('/')[-1])[0]
nx, nz = ft_dim[M]
print(nx, nz)
slipr = np.fromfile('sliprate.bin', 'f').reshape((-1, nz, nx))
nt_slip = slipr.shape[0]
with open("momrate.dat", 'rb') as f_mom:
    fig, ax = plt.subplots(4, 1, figsize=(6, 6))
    for i in range(4):
        np.random.seed(i)
        ix = np.random.randint(nx)
        iz = np.random.randint(nz)
        disp = iz * nx + ix
        print(disp)
        f_mom.seek((12 + 4 * 6 * nt_mom) * disp + 12, 0)
        # momrate has 6 components
        buf_mom = np.frombuffer(f_mom.read(4 * 6 * nt_mom), dtype='f')[::6]
        max_mom = max(abs(buf_mom))
        sum_mom = buf_mom.sum() * dt
        buf_mom = buf_mom / max_mom
        buf_slip = slipr[:, iz, ix]
        max_slip = max(abs(buf_slip))
        buf_slip = buf_slip / max_slip
        print("mom, max_momrate, max_sliprate: ", sum_mom, max_mom, max_slip)
        print("index: ", abs(buf_mom.argmax()), abs(buf_slip).argmax())

        ax[i].plot(np.arange(0, nt_mom * dt, dt), buf_mom, 'k', linewidth=1.5)
        ax[i].plot(np.arange(0, nt_slip * dt * ntskp, dt * ntskp), buf_slip, 'r--', linewidth=2)
        ax[i].set_ylabel('Amplitude')
    ax[0].legend(['momrate', 'sliprate'], loc=1)
    ax[3].set_xlabel('Time (s)')
    fig.savefig(f'comp_momrate_sliprate_{M}.png', dpi=300, bbox='tight', pad_inches=0.1)
