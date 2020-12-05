#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on Tues Oct 23 2018

@author: Zhifeng Hu<hzfmer94@gmail.com>
'''

import numpy as np
from numpy import sin, cos, pi, sqrt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.signal import resample
from struct import pack
import os
import sys
plt.tight_layout()

module = os.getcwd().split('/')[-1]

dt_ref = 0.05 
nt_ref = 4000
dx = 0.1
#ypos = 1173
#xoff = 20

nx = 5116
nz = 220
npt = nx * nz

sr = np.fromfile('../gp_021219_Scott_8.45_noplas_GPU_2hz/sliprate.bin', dtype='float32').reshape(nt_ref, nz, nx).T
tinit = np.fromfile('tinit.bin', dtype='float32').reshape((nz, nx))
fault_end = np.array(np.loadtxt("fault_loc.idx"))
x1 = int(fault_end[0,0])
x2 = int(fault_end[1,0])
if x1 > x2:
    sr = np.flip(sr, 0)
    tinit = np.flip(tinit, 1)
print(sr.shape)
#psr = np.zeros((nz, nx))
#slip = np.zeros((nz, nx))
slip = sr.sum(axis=2).T * dt_ref
psr = sr.max(axis=2).T
#for j in range(nx):
#    for i in range(nz):
#        slip[i,j] = np.sum(sr[j,i,:]) * dt_ref
#        psr[i,j] = sr[j,i,:].max()
#

fig, ax = plt.subplots(figsize=(9,3))
ax.plot(np.arange(nx) * dx, psr[0,:])
ax.set_xlabel('Strike (km)')
ax.set_ylabel('PSR (m/s)')
plt.savefig("surf_psr_{}.png".format(module), dpi=600, bbox_inches='tight', pad_inches=0.05)


fig = plt.figure(figsize=(9, 5))
ax0 = fig.add_axes([0.1, 0.58, 0.65, 0.3])
cax = ax0.imshow(slip, extent=[0, nx * dx, nz * dx, 0], 
                cmap='twilight', vmax=40, aspect=3)
xticks = ax0.get_xticks()
ax0.set_xticks(xticks[:-1])
ax0.set_ylabel("Depth (km)")
cb = plt.colorbar(cax, ax=ax0, orientation='horizontal', shrink = 0.6,
                  aspect = 18, pad=0.15)
cb.set_label("Slip (m)")
ax1 = fig.add_axes([0.78, 0.2, 0.13, 0.6])
ax1.plot(np.mean(slip, axis=1), np.arange(nz) * dx)
ax1.invert_yaxis()
ax1.set_xlabel("Average slip (m)")
#ax1.set_labelpad = 8
ax1.set_ylabel("Depth (km)")
#ax1.set_yticks([])
ax1.yaxis.set_label_position("right")
ax1.yaxis.set_ticks_position('right')
#ax1.tick_params(axis='y', which='both', labelleft='off', labelright='on')
ax2 = fig.add_axes([0.1, 0.15, 0.65, 0.3])
cnt = np.linspace(0, tinit.max(), 25) 
# ct = ax0.imshow(tinit, extent=[0, nx * dx, nz * dx, 0], cmap=cm.hot_r, aspect=4)
ct = ax2.contour(np.arange(0, nx) * dx, np.arange(nz) * dx, tinit, cnt, cmap=cm.hot_r) 
ax2.invert_yaxis()
# ax2.set_xticks(xticks[:-1])
ax2.set_ylabel("Depth (km)")
ax2.set_xlabel('Along strike (km)')
cb = plt.colorbar(ct, ax=ax2, orientation='horizontal', shrink = 0.6,
                  aspect = 18, pad=0.15)
cb.set_label("Rupture Time (s)")

plt.savefig("plot_slip_{}.png".format(module), dpi=600, bbox_inches='tight', pad_inches=0.05)
