#!/bin/usr/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import os

def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:

    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)

M = os.getcwd().split('/')[-1]

nx, ny = 6320, 4200
nxf, nzf = 5116, 220
dx = 0.1

# trace = np.fromfile('./fault_full_loc.bin', dtype='int32').reshape((nxf, 2))
trace = np.loadtxt('fault_loc.idx', dtype='int32')
trace_x = np.linspace(trace[0, 0], trace[1, 0], nxf)
trace_y = np.linspace(trace[0, 1], trace[1, 1], nxf)


with open('cities_name.txt', 'r') as f_name:
    cities_name = f_name.readlines()
cities_idx = dx * np.loadtxt('cities.idx', dtype='int')

T_range = [2, 3]
for T in T_range:
    vx = np.fromfile(f'gmrotD50_{1 / T:05.2f}Hz.bin', dtype='float32').reshape((ny, nx))
    print(np.max(vx))
    vx = np.flip(vx, 0) / 10

    fig, ax = plt.subplots(figsize=(6,6))
    c = ax.imshow(vx, extent=[0, nx * dx, 0, ny * dx], cmap=discrete_cmap(20, 'RdBu_r'),
            norm=LogNorm(vmin=0.01, vmax=np.max(vx) * 0.75))
    # c2 = ax[1].imshow(vy, extent=[0, nx * dx, 0, ny * dx], cmap='RdBu_r', vmax=10)
    # ax.scatter(trace_x * dx, trace_y * dx, 20, 'g', marker='*')
    cb = plt.colorbar(c, ax=ax, format='%.2f', label=f'SA-{T}s (g)', orientation='horizontal')
    ax.scatter(cities_idx[:, 0], cities_idx[:, 1], 50, 'w', marker='^')
    for i in range(len(cities_name)):
        ax.text(cities_idx[i, 0] - 45, cities_idx[i, 1] - 20, cities_name[i].strip('\n'), color='m')
    plt.tight_layout()
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    fig.savefig(f"SA_{T}s_with_trace_{M}.png", dpi=600, bbox_inches='tight', pad_inches=0.05)
