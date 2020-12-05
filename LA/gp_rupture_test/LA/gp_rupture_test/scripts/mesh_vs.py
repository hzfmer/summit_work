#!/ccs/proj/geo112/hzfmer/summit/opt/anaconda3/bin/python

import matplotlib.pyplot as plt
import numpy as np
import sys

nx = 6320
ny = 4200
nz = 400
dh = 0.1

vs = np.fromfile('mesh_info_for_GMPE_M6.35.bin',
        dtype='float32').reshape((ny, nx, 3))
vs = np.array(vs[:, :, 2])
print(np.sum(vs<500))

fig, ax = plt.subplots(figsize=(8,8))
plt.rc('image', origin='lower')
c = ax.imshow(vs, extent=[0, nx * dh, 0, ny * dh])
plt.colorbar(c, ax=ax, label="Vs (m/s)")
plt.tight_layout()
ax.set_xlabel('X (km)')
ax.set_ylabel('Y (km)')
fig.savefig(f"vs_gmpe_surface.png", dpi=600, bbox_inches='tight', pad_inches=0.05)

mesh = np.fromfile('../../data_for_Zhifeng/mesh/vs.media',
        dtype='float32').reshape((nz, ny, nx))
# ll = np.fromfile('/mnt/c/scratch/sciteam/hu2/LA/gp_rupture_test/gp_021219_Scott_6.35/surf.grid',
#         dypte='float64', count=2 * ny * nx).reshape(2, ny, nx)
vs = mesh[0, :, :]
del mesh
fig, ax = plt.subplots(figsize=(8,8))
plt.rc('image', origin='lower')
c = ax.imshow(vs, extent=[0, nx * dh, 0, ny * dh], vmax=2600)
plt.colorbar(c, ax=ax, label="Vs (m/s)")
plt.tight_layout()
ax.set_xlabel('X (km)')
ax.set_ylabel('Y (km)')
fig.savefig(f"vs_surface.png", dpi=600, bbox_inches='tight', pad_inches=0.05)
