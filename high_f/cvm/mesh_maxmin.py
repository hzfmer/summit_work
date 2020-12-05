#!/usr/bin/env python
'''
Max and min of the mesh.
'''
import matplotlib
import matplotlib.pyplot as plt
from pylab import *
import numpy as np
import sys
from sys import stdout
import os

X = 180000
Y = 135000
Z = 61440
dh20 = 20
dh = 8

# large-20m
nx = 9000
ny = 6750
# ext-large-20m
nz = 3072
nx = 9504
ny = 7020

# large-8m
nx, ny, nz = 19440, 14904, 160


fmax = 5
ppw = 7
start=0
nz = 50
skip=1
force = 0
f_mesh = "la_habra_large_cvmsi_20m_0.media"

if len(sys.argv) > 1:
    nx, ny, nz = map(int, (sys.argv[1], sys.argv[2], sys.argv[3]))
    start, skip = map(int, (sys.argv[4], sys.argv[5]))
    f_mesh = sys.argv[6]
    force = int(sys.argv[7])

print(nx, ny, nz, start, skip)

if force or not os.path.isfile('v_depth_ext_large.dat'):
    fidin = open(f_mesh,'rb')
    length = len(np.arange(start, nz, skip))
    vpmax = np.zeros(length)
    vsmax = np.zeros(length)
    vsmin = np.zeros(length)
    rhomin = np.zeros(length)
    depth = np.zeros(length)
    vpvsmax = np.zeros(length)
    vpvsmin = np.zeros(length)
    for ii, idx in enumerate(np.arange(start, nz, skip)):
        fidin.seek(idx * nx * ny * 12, 0) 
        print(idx * nx * ny * 12)
        dat = np.fromfile(fidin, count=3*nx*ny, dtype='f').reshape(-1, 3)
        depth[ii] = idx * dh
        vpmax[ii] = dat[:, 0].max()
        #vsmax[ii] = np.percentile(dat[:, 1], 0.01)
        vsmax[ii] = dat[:, 1].max()
        vpvsmax[ii] = np.max(dat[:, 0] / dat[:, 1])
        vpvsmin[ii] = np.min(dat[:, 0] / dat[:, 1])
        vsmin[ii] = dat[:, 1].min()
        rhomin[ii] = dat[:, 2].min()
        print(f"\rReading layer {ii}/{nz}: vsmin = {vsmin[ii]}")
    fidin.close()
    v=np.vstack((depth, vpmax, vsmax, vsmin, rhomin, vpvsmax, vpvsmin)).T
    np.savetxt('v_depth_ext_large.dat',v, fmt='%f',
            header="Depth (m) Vpmax (m/s) Vsmax (m/s) Vsmin (m/s) Rhomin\n")
else:
    v=np.genfromtxt('v_depth_ext_large.dat')
    vpmax = v[:, 0]
    vsmax = v[:, 0]
    vsmin = v[:, 1]

sys.exit(-1)
z1 = np.argwhere(vsmin > dh * fmax * 3 * ppw)[0] * dh20 /1000;
z2 = np.argwhere(vsmin > dh * fmax * 9 * ppw)[0] * dh20 /1000 - z1;
#vs.tofile('vs_tofile.dat', sep=" ", format='%s')
fig, ax = plt.subplots()
ax.invert_yaxis()
ax.plot(vsmin/fmax/dh, np.arange(nz) / 1000 * dh, linewidth=2)
ax.plot(np.full((nz,), ppw), np.arange(nz) * dh / 1000, '--', linewidth=3, label=f'ppw={ppw}')
ax.plot(np.full((nz,), 3*ppw), np.arange(nz) * dh / 1000, '--', linewidth=3, label=f'ppw={3 * ppw}, depth = {z1}km')
ax.plot(np.full((nz,), 9*ppw), np.arange(nz) * dh / 1000, '--', linewidth=3, label=f'ppw={9 * ppw}, depth = {(z1 + z2)}km')
ax.text(0.5, 0.5, f'Total number of grids: {int(X / dh * Y / dh * z1 * 1000/ dh)}, '
                  f'{int(X / dh * Y / dh * z2 * 1000 / dh / 27)}, '
                  f'{int(X / dh * Y / dh * (nz * dh20 / 1000 - z2) * 1000 / dh / 729)}, ')
ax.legend()
ax.set_xlabel('PPW')
ax.set_ylabel('Depth (km)')
fig.savefig('ppw_depth.png', dpi=600, bbox_inches='tight', pad_inches=0.1)

