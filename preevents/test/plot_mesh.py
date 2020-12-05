#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
        new_cmap = colors.LinearSegmentedColormap.from_list(
                        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
                                cmap(np.linspace(minval, maxval, n)))
        return new_cmap


plt.style.use('seaborn-talk')
plt.tight_layout()

nx = 80  
ny = 80  
nz = 75
dh = 0.004
nvar = 5
sz = 0.064
site_x = 40 
site_y = 41 

mesh = np.fromfile('mesh_0', dtype='f', count=nvar*nx*ny*nz).reshape(nz,
        ny, nx, nvar)
vs = mesh[:,:,:,1].swapaxes(0,2)

fig, ax = plt.subplots(figsize=(4,5)) 
plt.subplots_adjust(left=0.1, right=0.9, hspace=0.5)
cax1 = ax.imshow(np.squeeze(vs[site_x,::-1,:]).T, cmap='coolwarm_r',
                    interpolation='none', extent=[0, ny * dh, nz * dh, 0], aspect='auto')
cmap = plt.get_cmap('coolwarm_r')
new_cmap = truncate_colormap(cmap, 0, 0.5)
cb = plt.colorbar(cax1, ax=ax, orientation='horizontal', shrink=0.75, pad=0.2, 
        ticks=np.linspace(vs.min(), vs.max(),5))
cb.set_label('Vs (m/s)')
ax.scatter(site_y*dh, sz, s=100, c='g', marker=(5,2,0))
ax.set_xlabel('North-South (km)')
ax.set_ylabel('Depth (km)')
ax.set_aspect(2)
#top_vs = np.squeeze(vs[site_x, ::-1, 0:(nz//6)]).T
#cax0 = ax[0].imshow(top_vs, interpolation='none', extent=[0, ny * dh, nz // 6 * dh, 0]
#                    , aspect='auto', cmap=new_cmap)
#ax[0].set_xlabel('North-South (km)')
#ax[0].set_ylabel('Depth (km)')
#cb = plt.colorbar(cax0, ax=ax[0], ticks=np.linspace(top_vs.min(),
#            top_vs.max(), 3), orientation='vertical', shrink=1.5,
#            aspect=9)
#cb.set_label('Vs (m/s)')
fig.savefig('vs_perpend_GVDA_1d_h005l100.svg', dpi=600, bbox_inches='tight',
            pad_inches=0.1)

fig, ax = plt.subplots(figsize=(4,5))
plt.subplots_adjust(left=0.15, right=0.9, hspace=0.5)
cax1 = ax.imshow(np.squeeze(vs[:,site_y,:]).T, cmap='coolwarm_r',
                    interpolation='none', extent=[0, nx * dh, nz * dh, 0], aspect='auto')
cb = plt.colorbar(cax1, ax=ax, orientation='horizontal', shrink=0.75, pad=0.2, 
        ticks=np.linspace(vs.min(), vs.max(),5))
cb.set_label('Vs (m/s)')
ax.scatter(site_x*dh, sz, s=100, c='g', marker=(5,2,0))
ax.set_xlabel('West-East (km)')
ax.set_ylabel('Depth (km)')
ax.set_aspect(2)
#top_vs = np.squeeze(vs[:, site_y, 0:(nz//6)]).T
#cax0 = ax[0].imshow(top_vs, interpolation='none', extent=[0, nx * dh, nz // 6 * dh, 0]
#                    , aspect='auto', cmap=new_cmap)
#cb = plt.colorbar(cax0, ax=ax[0], ticks=np.linspace(top_vs.min(),
#                  top_vs.max(), 3), orientation='vertical', shrink=1.5,
#                  aspect=9)
#ax[0].set_xlabel('West-East (km)')
#ax[0].set_ylabel('Depth (km)')
#cb.set_label('Vs (m/s)')
fig.savefig('vs_horizon_GVDA_1d_h005l100.svg', dpi=600, bbox_inches='tight', pad_inches=0.1)
#fig.savefig('vs_het_horizon_TKCH05_h005l100_std10_sm.png', dpi=600, bbox_inches='tight', pad_inches=0.1)
