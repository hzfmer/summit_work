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
nt_ref = 500
dx = 0.1
ypos = 1173
xoff = 20

nx = 354
nz = 212
npt = nx * nz

sr = np.fromfile('sliprate.bin', dtype='float32').reshape(nt_ref, nz, nx).T
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
#fig = plt.figure(figsize=(9, 5))
#ax0 = fig.add_axes([0.1, 0.1, 0.7, 0.8])
#cax = ax0.imshow(psr, extent=[0, nx * dx, nz * dx, 0], aspect = 5,
#                cmap=cm.hot_r, vmax=12)
##cnt = np.arange(0., 100., 1.)
##ct = ax0.contour(grid[1,:,:]*dx, grid[0,:,:]*dx, slip)#, cnt)
##cb = plt.colorbar(ct, shrink=0.6, extend='both', orientation='horizontal')
##cb.set_label("Rupture start time (s)")
#ax0.set_xlabel("Strike (km)")
#xticks = ax0.get_xticks()
#ax0.set_xticks(xticks[:-1])
#ax0.set_ylabel("Depth (km)")
##ax.set_title("fault peak slip rate (PSR)")
# #   PCM=ax.get_children()[2] #get the mappable, the 1st and the 2nd are the x and y axes
#cb = plt.colorbar(cax, ax=ax0, orientation='horizontal', shrink = 0.6,
#                  aspect = 18, pad=0.15)
#cb.set_label("PSR (m/s)")
#ax1 = fig.add_axes([0.8, 0.34, 0.12, 0.38])
#ax1.plot(np.mean(psr, axis=1), np.arange(nz) * dx)
##print(psr.max(axis=1))
#ax1.invert_yaxis()
#ax1.set_xlabel("Average PSR (m/s)")
##ax1.set_labelpad = 8
#ax1.set_ylabel("Depth (km)")
##ax1.set_yticks([])
#ax1.yaxis.set_label_position("right")
#ax1.yaxis.set_ticks_position('right')
##ax1.tick_params(axis='y', which='both', labelleft='off', labelright='on')
#
#plt.savefig(f"plot_psr_{module}.svg", dpi=600, bbox_inches='tight', pad_inches=0.05)
#
fig = plt.figure(figsize=(9, 5))
ax0 = fig.add_axes([0.1, 0.1, 0.7, 0.8])
cax = ax0.imshow(slip, extent=[0, nx * dx, nz * dx, 0], 
                cmap=cm.hot_r, vmax=0.8)
#cnt = np.arange(0., 100., 1.)
#ct = ax0.contour(grid[1,:,:]*dx, grid[0,:,:]*dx, slip)#, cnt)
#cb = plt.colorbar(ct, shrink=0.6, extend='both', orientation='horizontal')
#cb.set_label("Rupture start time (s)")
ax0.set_xlabel("Strike (km)")
xticks = ax0.get_xticks()
ax0.set_xticks(xticks[:-1])
ax0.set_ylabel("Depth (km)")
#ax.set_title("fault peak slip rate (PSR)")
 #   PCM=ax.get_children()[2] #get the mappable, the 1st and the 2nd are the x and y axes
cb = plt.colorbar(cax, ax=ax0, orientation='horizontal', shrink = 0.6,
                  aspect = 18, pad=0.15)
cb.set_label("Slip (m)")
ax1 = fig.add_axes([0.75, 0.34, 0.15, 0.48])
ax1.plot(np.mean(slip, axis=1), np.arange(nz) * dx)
ax1.invert_yaxis()
ax1.set_xlabel("Average slip (m)")
#ax1.set_labelpad = 8
ax1.set_ylabel("Depth (km)")
#ax1.set_yticks([])
ax1.yaxis.set_label_position("right")
ax1.yaxis.set_ticks_position('right')
#ax1.tick_params(axis='y', which='both', labelleft='off', labelright='on')

plt.savefig("plot_slip_{}.png".format(module), dpi=600, bbox_inches='tight', pad_inches=0.05)
