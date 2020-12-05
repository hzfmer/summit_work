#!/ccs/proj/geo112/hzfmer/summit/opt/anaconda3/bin/python

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, sosfilt, sosfiltfilt
from scipy.fftpack import fft
import os
from filter_BU import filt_B
from read_params import read_params


rcparams = {'font.size': 16,
            'xtick.labelsize': 10,
            'ytick.labelsize': 10,
            'legend.fontsize': 10,
            'axes.titlesize': 14,
            'axes.labelsize': 14,
            'lines.linewidth': 1.5}
plt.rcParams.update(rcparams)

args = read_params()

ntskp, wstep = args['NTISKP'], args['WRITE_STEP']
skip = wstep * ntskp
dt = args['DT']
nt = np.int(args['TMAX'] / dt)
t = np.linspace(0, nt * dt, nt // ntskp)
nx, ny, nz = args['X'], args['Y'], args['Z']
base_x = args['SXRGO'].strip('\'')
base_y = args['SYRGO'].strip('\'')
base_z = args['SZRGO'].strip('\'')
print(base_x, type(base_x))


N = 4
np.random.seed(128)
ix = np.random.randint(nx, size=N)
iy = np.random.randint(ny, size=N)
vx = np.zeros((N, nt // ntskp))
for i in range(1, nt // skip):
    print("i=", i)
    idx = np.arange((i - 1) * wstep, i * wstep)
    v = np.fromfile(base_x + f'{i * skip:07d}', dtype='f').reshape(-1, ny, nx)
    for j in range(N):
        vx[j, idx] = v[:, iy[j], ix[j]]

vx = filt_B(vx, 1 / dt, 0, 10)
fig, ax = plt.subplots(N, 1, figsize=(6, 6))
plt.tight_layout()
fig.suptitle(f'VX at {N} sites')
for j in range(N):
    ax[j].plot(t, vx[j, :], label=f'X={ix[j]}, Y={iy[j]}')
    ax[j].set_ylabel('VX (m/s)')
    ax[j].legend(loc=1)
ax[N - 1].set_xlabel('Time (s)')
fig.savefig('ssmgrm_rand_sites.png', dpi=400, bbox='tight', pad_inches=0.1)

