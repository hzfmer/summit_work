#!/ccs/home/hzfmer/file_back/programs/anaconda3/bin/python

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, sosfiltfilt
from scipy.fftpack import fft
import os
import sys

plt.tight_layout()
rcparams = {'font.size': 16,
            'xtick.labelsize': 10,
            'ytick.labelsize': 10,
            'legend.fontsize': 14,
            'axes.titlesize': 16,
            'axes.labelsize': 14,
            'lines.linewidth': 2}
plt.rcParams.update(rcparams)

def comp_z1(vs_profile):
    if vs_profile.min() > 1000:
        return 0
    else:
        return np.argwhere(vs_profile <= 1000)[-1][0] * dh


nx = 6320
ny = 4200
nz = 400
dh = 0.1
vs = np.fromfile('../../data_for_Zhifeng/mesh/vs.media',
        dtype='float32').reshape((nz, ny, nx))


for mag in [6.35, 6.45, 7.35, 8.15, 8.45]:
    M = str(mag).replace('.', '_')
    ll = np.genfromtxt(f'stat_M{M}.txt', dtype='int', skip_header=1)
    ll = ll - 1  # start from zero
    nsite = len(ll)
    dist = np.fromfile(f'site_rjb_M{M}.bin', dtype='float32')
    z1 = np.zeros((nsite, 1), dtype='float32')
    vs30 = np.zeros((nsite, 1), dtype='float32')
    vs0 = np.zeros((nsite, 1), dtype='float32')
    for i in range(nsite):
        vs_local = vs[:, ll[i][1], ll[i][0]]
        vs0[i] = vs_local[0]
        z1[i] = comp_z1(vs_local)
        if vs_local[0] < 501:
            vs30[i] = 250.0
        else:
            vs30[i] = 760.0
    data = np.column_stack((dist, z1, vs30, vs0))
    data.tofile(f'site_info_for_GMPE_M{M}.bin', format='float32')


