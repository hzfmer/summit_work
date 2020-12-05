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
nz = 160
skip=1

if not os.path.isfile('v_depth_ext_large.dat'):
    fidin = open('mesh_large_8m_orig.bin_0','rb')
    vpmax = np.zeros(nz)
    vsmax = np.zeros(nz)
    vsmin = np.zeros(nz)
    rhomin = np.zeros(nz)
    depth = np.zeros(nz)
    for ii in np.arange(nz):
        idx = ii * skip + start
        fidin.seek(idx * nx * ny * 12, 0) 
        buf = fidin.read(12 * nx * ny)
        dat = np.frombuffer(buf, dtype='f', count=3*nx*ny).reshape(-1,3)
        depth[ii] = idx * dh
        vpmax[ii] = dat[:, 0].max()
        vsmax[ii] = np.percentile(dat[:, 1], 0.01)
        vsmin[ii] = dat[:, 1].min()
        rhomin[ii] = dat[:, 2].min()
        print(f"\rReading layer {ii}/{nz}: vsmin = {vsmin[ii]}")
    fidin.close()
    v=np.vstack((depth, vpmax, vsmax, vsmin, rhomin)).T
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

