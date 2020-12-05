#!/usr/bin/env python

import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import resample
import matplotlib.pyplot as plt
import sys
from mpl_toolkits import mplot3d
import os

M = os.getcwd().split('/')[-1]
nxf, nzf = 5116, 220
dh = 0.1
idx = np.fromfile('fault_idx.bin', dtype='int32').reshape((nzf, nxf, 2))
z = np.tile(np.arange(nzf), (nxf, 1)).T
fig, ax = plt.subplots(figsize=(6,6))
ax = plt.axes(projection='3d')
ax.plot_surface(dh * idx[:,:,0], dh * idx[:,:,1], dh * z)
ax.scatter(3704 * dh, 2131 * dh, s=40, marker='s', label='(3704, 2131)')
ax.scatter(2425 * dh, 1855 * dh, s=40, marker='s', label='(2425, 1855)')
ax.legend()
ax.invert_zaxis()
fig.savefig(f"fault_idx_{M}.png", dpi=600, bbox_inches='tight', pad_inches=0.05)

loc = np.loadtxt('fault_surf_loc.txt', dtype='float32')
loc2 = np.fromfile('./fault_full_loc.bin', dtype='float32')
loc2 = loc2.reshape((-1, 2))
fig, ax = plt.subplots(figsize=(6,6))
ax.scatter(loc[:,0], loc[:,1]) 
ax.scatter(loc2[:,0], loc2[:,1])
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
fig.savefig(f"fault_surf_loc_{M}.png", dpi=600, bbox_inches='tight', pad_inches=0.05)


