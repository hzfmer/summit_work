#!/ccs/home/hzfmer/file_back/programs/anaconda3/bin/python

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

nx = 100  
ny = 100  
nz = 1000
dh = 0.002
nvar = 5
sz = 30 * dh
site_x = 50 
site_y = 50 
zz = 200

mesh = np.fromfile('mesh_0', dtype='f', count=nvar*nx*ny*nz).reshape(nz,
        ny, nx, nvar)
vs = mesh[:,:,:,1].swapaxes(0,2)

fig, ax = plt.subplots(figsize=(4,5)) 
plt.subplots_adjust(left=0.1, right=0.9, hspace=0.5)
cax1 = ax.imshow(np.squeeze(vs[site_x,::-1,:zz]).T, cmap='coolwarm_r',
                    interpolation='none', extent=[0, ny * dh, zz * dh, 0], aspect='auto')
cmap = plt.get_cmap('coolwarm_r')
new_cmap = truncate_colormap(cmap, 0, 0.5)
cb = plt.colorbar(cax1, ax=ax, orientation='horizontal', shrink=0.75, pad=0.2, 
        ticks=np.linspace(vs.min(), vs.max(),5))
cb.set_label('Vs (m/s)')
ax.scatter(site_y*dh, sz, s=100, c='g', marker=(5,2,0))
ax.set_xlabel('North-South (km)')
ax.set_ylabel('Depth (km)')
#ax.set_aspect(2)
fig.savefig('mesh_perpend.png', dpi=600, bbox_inches='tight', pad_inches=0.1)

fig, ax = plt.subplots(figsize=(4,5))
plt.subplots_adjust(left=0.15, right=0.9, hspace=0.5)
cax1 = ax.imshow(np.squeeze(vs[:,site_y,:zz]).T, cmap='coolwarm_r',
                    interpolation='none', extent=[0, nx * dh, zz * dh, 0], aspect='auto')
cb = plt.colorbar(cax1, ax=ax, orientation='horizontal', shrink=0.75, pad=0.2, 
        ticks=np.linspace(vs.min(), vs.max(),5))
cb.set_label('Vs (m/s)')
ax.scatter(site_x*dh, sz, s=100, c='g', marker=(5,2,0))
ax.set_xlabel('West-East (km)')
ax.set_ylabel('Depth (km)')
#ax.set_aspect(2)
fig.savefig('mesh_horizontal.png', dpi=600, bbox_inches='tight', pad_inches=0.1)
