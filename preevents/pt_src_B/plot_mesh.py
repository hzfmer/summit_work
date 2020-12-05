#!/usr/bin/env python
#!/ccs/home/hzfmer/file_back/programs/anaconda2/bin python

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

nx = 1000
ny = 480
nz = 800
dh = 0.005
nvar = 5
site_x = 600
site_y = 237
sz = 600

mesh = np.fromfile('mesh_0', dtype='f', count=nvar*nx*ny*nz).reshape(nz,
        ny, nx, nvar)

vs = mesh[:,:,:,1].swapaxes(0,2)

fig, ax = plt.subplots(2, 1, figsize=(5,6), 
                       gridspec_kw=dict(width_ratios=[1],
                       height_ratios=[1,4]))
plt.subplots_adjust(left=0.15, right=0.85)
cax1 = ax[1].imshow(np.squeeze(vs[site_x,::-1,:]).T, cmap='coolwarm_r',
                    interpolation='none', extent=[0, 2.4, 4, 0], aspect='auto')
cmap = plt.get_cmap('coolwarm_r')
new_cmap = truncate_colormap(cmap, 0, 0.5)
cb = plt.colorbar(cax1, ax=ax[1], orientation='vertical', shrink=0.7)
cb.set_label('Vs (m/s)')
#ax[1].scatter(site_y*dh, sz*dh, s=300, c='g', marker=(5,2,0))
ax[1].set_xlabel('South-North (km)')
ax[1].set_ylabel('Depth (km)')

top_vs = np.squeeze(vs[site_x, ::-1, 0:(nz//8)]).T
cax0 = ax[0].imshow(top_vs, interpolation='none', extent=[0, 2.4, 0.5, 0]
                    , aspect='auto', cmap=new_cmap)
ax[0].set_xlabel('South-North (km)')
ax[0].set_ylabel('Depth (km)')
cb = plt.colorbar(cax0, ax=ax[0], ticks=np.linspace(top_vs.min(),
            top_vs.max(), 3), orientation='vertical', shrink=1.5,
            aspect=9)
cb.set_label('Vs (m/s)')
fig.savefig('vs_het_perpend_TKCH05.png', dpi=600, bbox_inches='tight',
            pad_inches=0.1)

fig, ax = plt.subplots(2, 1, figsize=(6,6),
                       gridspec_kw=dict(width_ratios=[1],
                       height_ratios=[1,3]))
plt.subplots_adjust(left=0.15, right=0.9, hspace=0.5)
cax1 = ax[1].imshow(np.squeeze(vs[:,site_y,:]).T, cmap='coolwarm_r',
                    interpolation='none', extent=[0, 5, 4, 0], aspect='auto')
cb = plt.colorbar(cax1, ax=ax[1], orientation='vertical', shrink=0.7)
cb.set_label('Vs (m/s)')
#ax[1].scatter(site_x*dh, sz*dh, s=300, c='g', marker=(5,2,0))
ax[1].set_xlabel('West-East (km)')
ax[1].set_ylabel('Depth (km)')

top_vs = np.squeeze(vs[:, site_y, 0:(nz//8)]).T
cax0 = ax[0].imshow(top_vs, interpolation='none', extent=[0, 5, 0.5, 0]
                    , aspect='auto', cmap=new_cmap)
cb = plt.colorbar(cax0, ax=ax[0], ticks=np.linspace(top_vs.min(),
                  top_vs.max(), 3), orientation='vertical', shrink=1.5,
                  aspect=9)
ax[0].set_xlabel('West-East (km)')
ax[0].set_ylabel('Depth (km)')
cb.set_label('Vs (m/s)')
fig.savefig('vs_het_horizon_TKCH05.png', dpi=600, bbox_inches='tight',
            pad_inches=0.1)
