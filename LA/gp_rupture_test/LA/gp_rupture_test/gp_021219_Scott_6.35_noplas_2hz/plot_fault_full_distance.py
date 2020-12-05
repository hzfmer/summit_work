#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import os

nx, ny = 6320, 4200
dh = 0.1
M = os.getcwd().split('_')[3]

dist = np.fromfile('fault_distance.bin', dtype='float32').reshape((ny, nx))
fig, ax = plt.subplots(figsize=(8, 8))
im = ax.imshow(dist, extent=[0, nx * dh, 0, ny * dh], cmap='RdBu_r', origin='lower')
ax.set_xlabel('X (km)')
ax.set_ylabel('Y (km)')
cb = plt.colorbar(im, ax=ax, orientation='horizontal', shrink=0.6)
cb.set_label('Rjb (km)')
fig.savefig(f"fault_nearest_dist_M{M}.png", dpi=600, bbox_inches='tight', pad_inches=0.05)

