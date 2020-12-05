#!/usr/bin/env python

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import struct
from scipy.signal import butter, sosfilt, sosfiltfilt
import sys

def filt_B(data_in, fs, lowcut=0, highcut=1, order=4, causal=True):
    if highcut >= fs/2:
        return data_in
    sz = data_in.shape
    if len(sz) > 1 and sz[0] > sz[1]:
        data_in = data_in.T
    data_out = np.zeros_like(data_in)
    nyq = fs / 2  
    low = lowcut / nyq
    high = highcut / nyq
    if low == 0:  # lowpass
        sos = butter(order, high, analog=False, btype='low', output='sos')
    else:   # bandpass
        sos = butter(order, [low, high], analog=False, btype='band', output='sos')
    if len(sz)>1:
        for i in range(sz[-1]):
            if causal == True:
                data_out[i, :] = sosfilt(sos, data_in[i, :])
                data_out[i, :] = sosfilt(sos, data_out[i, :])
            else:
                data_out[i, :] = sosfiltfilt(sos, data_in[i, :])
        if sz[0] > sz[1]:
            data_out = data_out.T
    else:
        data_out = sosfiltfilt(sos, data_in)
    return data_out


nx = 6320
ny = 4200
nz = 400
nt = 10000
dh = 0.1
dt = 0.05

# sx = np.fromfile(f'output_sfc/SX{nt}', dtype='float32').reshape((nt, ny, nx))
# sy = np.fromfile(f'output_sfc/SY{nt}', dtype='float32').reshape((nt, ny, nx))
# sz = np.fromfile(f'output_sfc/SZ{nt}', dtype='float32').reshape((nt, ny, nx))

fname = f'output_sfc/SX{nt:07d}'
fid = open(fname, 'rb')
for i in range(0, 400, 10):
    break
    bskip = 4 * i * ny * nx
    fid.seek(bskip, 0)
    sx = struct.unpack('%df' % (nx * ny), fid.read(4 * ny * nx))
    sx = np.array(sx).reshape((ny, nx))
    fig, ax = plt.subplots(figsize=(8,8))
    plt.rc('image', origin='lower')
    c = ax.imshow(sx, extent=[0, nx * dh, 0, ny * dh], cmap='RdBu_r', 
            norm=mpl.colors.SymLogNorm(linthresh=0.03, linscale=0.03))
    ax.set_xlim((300, 450))
    ax.set_ylim((130, 260))
    ax.set_xlabel("X (km)")
    ax.set_ylabel("Y (km)")
    plt.colorbar(c, ax=ax, orientation='horizontal')
    fig.savefig(f'Snapshot_vx_timestep{i}.png', dpi=400)
    plt.close('all')


# site BVH
ix = 3455
iy = 1563
nt = nt // 10
vx = np.fromfile(fname, dtype='float32').reshape((nt, ny, nx))
vx = vx[:, iy, ix]
vx_extrts = np.loadtxt('output_sfc/SX3455_1563_0001.dat')
np.column_stack((vx, vx_extrts)).tofile('BVH2.txt', format='%.4f %.4f\n')
# vx.tofile('BVH.txt', sep='\n', format='%.4f')
# sys.exit(-1)
fig, ax = plt.subplots(3, 1, figsize=(8,8))
plt.subplots_adjust(hspace=0.3)
ax[0].plot(np.arange(nt) * dt, vx, label='Python')
ax[0].plot(np.arange(nt) * dt, vx_extrts, label='extrts')
ax[1].plot(np.arange(nt) * dt, filt_B(vx, 20, 0.01, 1, True))
ax[1].plot(np.arange(nt) * dt, filt_B(vx_extrts, 20, 0.01, 1, True))
ax[2].plot(np.arange(nt) * dt, filt_B(vx, 20, 0.01, 1, False))
for i in range(3):
    ax[i].set_ylabel('V (m/s)')
    ax[i].set_xlabel('Time (s)')
ax[0].legend()
fig.savefig(f"Time_History_BVH.png", dpi=600, bbox_inches='tight', pad_inches=0.05)
plt.close('all')
fid.close()
