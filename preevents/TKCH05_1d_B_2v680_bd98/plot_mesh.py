#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import os


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
        new_cmap = colors.LinearSegmentedColormap.from_list(
                        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
                                cmap(np.linspace(minval, maxval, n)))
        return new_cmap


plt.style.use('seaborn-talk')
plt.tight_layout()

nx = 64  
ny = 64 
nz = 400 
dh = 0.0025
nvar = 5
site_x = 33
site_y = 32
sz = 0.1

model = os.getcwd().split('/')[-1]
mesh = np.fromfile('mesh_0', dtype='f', count=nvar*nx*ny*nz).reshape(nz,
        ny, nx, nvar)
vs = mesh[:,:,:,1].swapaxes(0,2)


fig, ax = plt.subplots(figsize=(3,6))
plt.subplots_adjust(left=0, right=0.7)
cax1 = ax.imshow(np.squeeze(vs[site_x, ::-1,:]).T,cmap='coolwarm_r',
                interpolation='none',extent=[0, ny * dh, nz * dh, 0])
cb = plt.colorbar(cax1, ax=ax, orientation='vertical', shrink=0.7)
ax.scatter(site_y*dh, sz,s=100, c='g', marker=(5,2,0))
cb.set_label('Vs (m/s)')
ax.set_xlabel('North-South (km)')
ax.set_ylabel('Depth (km)')
fig.savefig(f'vs_perpend_{model}.png', dpi=600, bbox_inches='tight',
            pad_inches=0.1)

fig, ax = plt.subplots(figsize=(3,6))
plt.subplots_adjust(left=0, right=0.7)
cax1 = ax.imshow(np.squeeze(vs[:,site_y,:]).T,cmap='coolwarm_r',
                interpolation='none',extent=[0, nx * dh, nz * dh, 0])
cb = plt.colorbar(cax1, ax=ax, orientation='vertical', shrink=0.7)
ax.scatter(site_x*dh, sz, s=100, c='g', marker=(5,2,0))
cb.set_label('Vs (m/s)')
ax.set_xlabel('West-East (km)')
ax.set_ylabel('Depth (km)')
fig.savefig(f'vs_horizon_{model}.png', dpi=600, bbox_inches='tight',
            pad_inches=0.1)
