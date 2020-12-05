#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import sys
from numpy import pi, exp


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
        new_cmap = colors.LinearSegmentedColormap.from_list(
                        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
                                cmap(np.linspace(minval, maxval, n)))
        return new_cmap


plt.style.use('seaborn-talk')
plt.tight_layout()

nx = 400
ny = 400
nz = 400
dh = 0.0025
nvar = 5
sz = 0.1
site_x = 200
site_y = 200

mesh = np.fromfile('mesh_3d', dtype='f', count=nvar*nx*ny*nz).reshape(nz,
        ny, nx, nvar)

sigma = 20
gauss = np.zeros((nx, ny))
for j in range(nx):
    for i in range(ny):
        gauss[i, j] = 1 / (2 * pi * sigma ** 2) * exp(-((j - site_x) ** 2 + (i - site_y) ** 2) / (2 * sigma ** 2))
   
print(mesh[::10, 200, 200, 1])
vp = np.zeros((nz, 1))
vs = np.zeros((nz, 1))
ro = np.zeros((nz, 1))
qp = np.zeros((nz, 1))
qs = np.zeros((nz, 1))
for k in range(nz):
    vp[k] = np.tensordot(mesh[k, :, :, 0], gauss)
    vs[k] = np.tensordot(mesh[k, :, :, 1], gauss)
#vs[k] = mesh[k, site_y, site_x, 1] * np.mean(mesh[k, :, :, 1] / np.tensordot(mesh[k, :, :, 1], gauss) 
    ro[k] = np.tensordot(mesh[k, :, :, 2], gauss)
    qp[k] = np.tensordot(mesh[k, :, :, 3], gauss)
    qs[k] = np.tensordot(mesh[k, :, :, 4], gauss)

print(vs[::10])
vs_surf = 140
scale = vs_surf / vs[0]
[vp, vs, ro, qp, qs] = [vp, vs, ro, qp, qs] * scale
np.savetxt('mesh_TKCH05_gauss40_vsurf140_1d.txt', np.hstack([vp, vs, ro, qp, qs]), fmt='%f')
nxx = 64
nyy = 64
site_xx = 33
site_yy = 32
[vp, vs, ro, qp, qs] = [np.tile(x.T, [nxx, nyy, 1]).swapaxes(0, 2) for x in [vp, vs, ro, qp, qs]]

dat = np.vstack(list(map(np.ravel, (vp, vs, ro, qp, qs)))).T.astype('f')
dat.tofile('mesh_TKCH05_gauss40_vsurf140_1d.bin')


fig, ax = plt.subplots(figsize=(3,6))
plt.subplots_adjust(left=0, right=0.7)
cax1 = ax.imshow(np.squeeze(vs[:,site_yy,:]), cmap='coolwarm_r',
                            interpolation='none', extent=[0, nxx * dh, nz * dh, 0], aspect='auto')
cb = plt.colorbar(cax1, ax=ax, orientation='vertical', shrink=0.7)
cb.set_label('Vs (m/s)')
ax.scatter(site_xx*dh, sz, s=100, c='g', marker=(5,2,0))
ax.set_xlabel('West-East (km)')
ax.set_ylabel('Depth (km)')
fig.savefig('vs_horizon_TKCH05_gauss40_vsurf140_1d.svg', dpi=600, bbox_inches='tight', pad_inches=0.1)
