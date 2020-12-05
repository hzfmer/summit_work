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

nx = 1000
ny = 1000
nz = 100
dh = 0.004
nvar = 5
sz = 0.064
site_x = 501
site_y = 505

mesh = np.fromfile('mesh_0', dtype='f', count=nvar*nx*ny*nz).reshape(nz,
        ny, nx, nvar)

sigma = 30
gauss = np.zeros((nx, ny))
for j in range(nx):
    for i in range(ny):
        gauss[i, j] = 1 / (2 * pi * sigma ** 2) * exp(-((j - site_x) ** 2 + (i - site_y) ** 2) / (2 * sigma ** 2))
    
vp = np.zeros((nz, 1))
vs = np.zeros((nz, 1))
ro = np.zeros((nz, 1))
qp = np.zeros((nz, 1))
qs = np.zeros((nz, 1))
print(mesh.shape)
print(mesh[0,:,:,0].shape)
print(gauss.shape)
print(np.tensordot(mesh[1, :, :, 0], gauss))
for k in range(nz):
    vp[k] = np.tensordot(mesh[k, :, :, 0], gauss)
    vs[k] = np.tensordot(mesh[k, :, :, 1], gauss)
#vs[k] = mesh[k, site_y, site_x, 1] * np.mean(mesh[k, :, :, 1] / np.tensordot(mesh[k, :, :, 1], gauss) 
    ro[k] = np.tensordot(mesh[k, :, :, 2], gauss)
    qp[k] = np.tensordot(mesh[k, :, :, 3], gauss)
    qs[k] = np.tensordot(mesh[k, :, :, 4], gauss)

vs_surf = 200
scale = vs_surf / vs[0]
[vp, vs, ro, qp, qs] = [vp, vs, ro, qp, qs] * scale
np.savetxt('mesh_gvda_4m_gauss100_vsurf200_1d.txt', np.hstack([vp, vs, ro, qp, qs]), fmt='%f')
nxx = 80
nyy = 80
site_xx = 40
site_yy = 41
[vp, vs, ro, qp, qs] = [np.tile(x.T, [nxx, nyy, 1]).swapaxes(0, 2) for x in [vp, vs, ro, qp, qs]]

dat = np.vstack(list(map(np.ravel, (vp, vs, ro, qp, qs)))).T.astype('f')
dat.tofile('mesh_gvda_4m_gauss100_vsurf200_1d.bin')


fig, ax = plt.subplots(figsize=(4,5))
plt.subplots_adjust(left=0.1, right=0.9, hspace=0.5)
cax1 = ax.imshow(np.squeeze(vs[:, ::-1, site_xx]), cmap='coolwarm_r',
                            interpolation='none', extent=[0, nyy * dh, nz * dh, 0], aspect='auto')
cmap = plt.get_cmap('coolwarm_r')
new_cmap = truncate_colormap(cmap, 0, 0.5)
cb = plt.colorbar(cax1, ax=ax, orientation='horizontal', shrink=0.75, pad=0.2,
                ticks=np.linspace(vs.min(), vs.max(),5))
cb.set_label('Vs (m/s)')
ax.scatter(site_yy*dh, sz, s=100, c='g', marker=(5,2,0))
ax.set_xlabel('North-South (km)')
ax.set_ylabel('Depth (km)')
ax.set_aspect(2)
fig.savefig('vs_perpend_GVDA_gauss100_vsurf200_1d.svg', dpi=600, bbox_inches='tight',
                        pad_inches=0.1)


fig, ax = plt.subplots(figsize=(4,5))
plt.subplots_adjust(left=0.15, right=0.9, hspace=0.5)
cax1 = ax.imshow(np.squeeze(vs[:,site_yy,:]), cmap='coolwarm_r',
                            interpolation='none', extent=[0, nxx * dh, nz * dh, 0], aspect='auto')
cb = plt.colorbar(cax1, ax=ax, orientation='horizontal', shrink=0.75, pad=0.2,
                ticks=np.linspace(vs.min(), vs.max(),5))
cb.set_label('Vs (m/s)')
ax.scatter(site_xx*dh, sz, s=100, c='g', marker=(5,2,0))
ax.set_xlabel('West-East (km)')
ax.set_ylabel('Depth (km)')
ax.set_aspect(2)
fig.savefig('vs_horizon_GVDA_gauss100_vsurf200_1d.svg', dpi=600, bbox_inches='tight', pad_inches=0.1)
