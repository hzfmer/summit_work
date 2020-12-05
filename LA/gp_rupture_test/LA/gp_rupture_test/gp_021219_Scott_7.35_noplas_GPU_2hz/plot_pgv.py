#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import os
import re

def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:

    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)

M = re.findall("\d+\.\d+", os.getcwd().split('/')[-1])[0]

params = {'6.35': (354, 212), '7.35': (980, 240), '8.45': (5116, 220)}
nx, ny = 6320, 4200
nxf, nzf = params[M]
dx = 0.1

with open('cities_name.txt', 'r') as f_name:
    cities_name = f_name.readlines()
cities_idx = dx * np.loadtxt('cities.idx', dtype='int')
vx = np.fromfile('./peak_velocity_H_01.0Hz.bin', dtype='float32').reshape((ny, nx))
vy = np.fromfile('./peak_velocity_Z_01.0Hz.bin', dtype='float32').reshape((ny, nx))
vx = np.flip(vx, 0)
vy = np.flip(vy, 0)
trace = np.fromfile('./fault_idx.bin', dtype='int32').reshape((nzf, nxf, 2))

fig, ax = plt.subplots(figsize=(6,6))
c = ax.imshow(vx, extent=[0, nx * dx, 0, ny * dx], cmap=discrete_cmap(20, 'RdBu_r'),
        norm=LogNorm(vmin=0.01, vmax=3))
# c2 = ax[1].imshow(vy, extent=[0, nx * dx, 0, ny * dx], cmap='RdBu_r', vmax=10)
ax.scatter(trace[0, :, 0] * dx, trace[0, :, 1] * dx, 50, 'g', marker='*')
cb = plt.colorbar(c, ax=ax, format='%.2f', label='PGV (m/s)', orientation='horizontal')
ax.scatter(cities_idx[:, 0], cities_idx[:, 1], 20, 'w', marker='^')
for i in range(len(cities_name)):
    ax.text(cities_idx[i, 0] - 45, cities_idx[i, 1] - 20, cities_name[i].strip('\n'), color='b')
plt.tight_layout()
ax.set_xlabel('X (km)')
ax.set_ylabel('Y (km)')
fig.savefig(f"PGV_with_trace_{M}.png", dpi=600, bbox_inches='tight', pad_inches=0.05)
