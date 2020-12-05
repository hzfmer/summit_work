#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  5 00:10:54 2018

@author: zhh076
"""

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
import glob
import re

plt.tight_layout()

def compmt(stk, dip, rake, area, mu):
    # stk, dip, rake: degree
    # area: cm^2
    # rho: g/cm^3
    # vs: cm/s
    stk = (stk -theta) * pi / 180
    rake = rake * pi / 180
    dip = dip * pi / 180
    areamu = mu * area 
    yy = -(sin(dip) * cos(rake) * sin(2 * stk) + sin(2 * dip) \
         * sin(rake) * sin(stk) ** 2)
    xy = sin(dip) * cos(rake) * cos(2 * stk) + 0.5 * sin(2 * dip) \
         * sin(rake) * sin(2 * stk)
    yz = cos(dip) * cos(rake) * cos(stk) + cos(2 * dip) \
         * sin(rake) * sin(stk)
    xx = sin(dip) * cos(rake) * sin(2 * stk) - sin(2 * dip) \
         * sin (rake) * cos(stk) ** 2
    xz = cos(dip) * cos(rake) * sin(stk) - cos(2 * dip) \
         * sin(rake) * cos(stk)
    zz = sin(2 * dip) * sin(rake)
    xx = xx * areamu
    yy = yy * areamu
    zz = zz * areamu
    xy = xy * areamu
    xz = xz * areamu
    yz = yz * areamu
    return xx, yy, zz, xz, yz, xy


M = re.findall("\d+\.\d+", os.getcwd().split('/')[-1])[0]

theta = 35
dx = 0.2
dt = 0.05
f = open(glob.glob('./*.srf')[0],'r')
f.readline()
f.readline()
token = f.readline()
nx = int(token.split()[2])
nz = int(token.split()[3])
print("nx, nz", nx, nz)
npt = nx * nz
f.readline()
f.readline()

fault_end = np.array(np.loadtxt("fault_loc.idx"))
x1 = int(fault_end[0,0])
x2 = int(fault_end[1,0])
param = np.loadtxt('subfaults.idx').reshape((-1, 8))
mu = param[:, 4].reshape((2 * nz, 2 * nx))
slipr = np.fromfile('sliprate.bin', dtype='f').reshape((4000, 2 * nz, 2 * nx))
slip = np.sum(np.abs(slipr), 0) * dt
moment = np.sum(np.multiply(slip, mu)) * (dx * 1000 / 2) ** 2
print(f'Total moment from sliprate.bin is: {np.sum(moment)}')
mag = 2. / 3. * (np.log10(moment) - 9.1)
print(f'Magnitude is {mag}')

if x1 > x2:
    mu = np.flip(mu, 0)
mu = np.ravel(mu[::2, ::2])

mom = np.zeros((npt,), dtype=float)
rt = np.zeros((npt,), dtype=float)
trise = np.zeros((npt,), dtype=float)
area = np.zeros((npt,), dtype=float)
tslip = np.zeros((npt,), dtype=float)
strike = np.zeros((nz, nx), dtype=float)
dip = np.zeros((nz, nx), dtype=float)
rake = np.zeros((nz, nx), dtype='float')

for k in range(npt):
    nl1 = f.readline().split()
    nl2 = f.readline().split()
    strike[k//nx, k%nx] = float(nl1[3]) - theta
    if strike[k//nx, k%nx] < 80:
        print(nl1)
        sys.exit(-1)
    dip[k//nx, k%nx] = float(nl1[4])
    rake[k//nx, k%nx] = float(nl2[0])
     
    area = float(nl1[5]) / 1.e4
    tinit= float(nl1[6])
    dt = float(nl1[7])
    
    slip1 = float(nl2[1]) / 1.e2
    nt1 = int(nl2[2])
    nskip1 = int(np.ceil(nt1 / 6.))

    for l in range(nskip1):
        f.readline()   
    
    mom[k] = slip1 * mu[k] * area
    rt[k] = tinit
    
    
f.close()
print(f'Total moment from srf is: {np.sum(mom)}')
mag = 2. / 3. * (np.log10(sum(mom)) - 9.1)
print(f'Magnitude is {mag}')


fig = plt.figure(figsize=(9, 5))
ax0 = fig.add_axes([0.1, 0.1, 0.6, 0.8])
cax = ax0.imshow(strike, extent=[0, nx * dx, nz * dx, 0], aspect = 'auto',
                cmap=cm.coolwarm_r)
ax0.set_xlabel("Along Strike (km)")
xticks = ax0.get_xticks()
ax0.set_xticks(xticks[:-1])
ax0.set_ylabel("Depth (km)")
cb = plt.colorbar(cax, ax=ax0, orientation='horizontal', shrink = 0.6,
                  aspect = 18, pad=0.15)
cb.set_label("Strike")
ax1 = fig.add_axes([0.78, 0.34, 0.12, 0.48])
ax1.plot(np.mean(strike, axis=1), np.arange(nz) * dx)
ax1.invert_yaxis()
ax1.set_xlabel("Avg Strike Angle")
ax1.set_ylabel("Depth (km)")
ax1.yaxis.set_label_position("right")
ax1.yaxis.set_ticks_position('right')

plt.savefig(f"strike_fault_{M}.png", dpi=600, bbox_inches='tight', pad_inches=0.05)

fig = plt.figure(figsize=(9, 5))
ax0 = fig.add_axes([0.1, 0.1, 0.6, 0.8])
cax = ax0.imshow(dip, extent=[0, nx * dx, nz * dx, 0], aspect = 'auto',
                cmap=cm.hot_r)
ax0.set_xlabel("Along Strike (km)")
xticks = ax0.get_xticks()
ax0.set_xticks(xticks[:-1])
ax0.set_ylabel("Depth (km)")
cb = plt.colorbar(cax, ax=ax0, orientation='horizontal', shrink = 0.6,
                  aspect = 18, pad=0.15)
cb.set_label("Dip")
ax1 = fig.add_axes([0.78, 0.34, 0.12, 0.48])
ax1.plot(np.mean(dip, axis=1), np.arange(nz) * dx)
ax1.invert_yaxis()
ax1.set_xlabel("Avg Dip Angle")
ax1.set_ylabel("Depth (km)")
ax1.yaxis.set_label_position("right")
ax1.yaxis.set_ticks_position('right')
ax1.tick_params(axis='y', which='both', labelleft=False, labelright=True)
fig.savefig(f"dip_fault_{M}.png", dpi=600, bbox_inches='tight', pad_inches=0.05)

fig = plt.figure(figsize=(9, 5))
ax0 = fig.add_axes([0.1, 0.1, 0.6, 0.8])
cax = ax0.imshow(rake, extent=[0, nx * dx, nz * dx, 0], aspect = 'auto',
                cmap=cm.hot_r)
ax0.set_xlabel("Along Strike (km)")
xticks = ax0.get_xticks()
ax0.set_xticks(xticks[:-1])
ax0.set_ylabel("Depth (km)")
cb = plt.colorbar(cax, ax=ax0, orientation='horizontal', shrink = 0.6,
                  aspect = 18, pad=0.15)
cb.set_label("Rake")
ax1 = fig.add_axes([0.78, 0.34, 0.12, 0.48])
ax1.plot(np.mean(rake, axis=1), np.arange(nz) * dx)
ax1.invert_yaxis()
ax1.set_xlabel("Avg Rake Angle")
ax1.set_ylabel("Depth (km)")
ax1.yaxis.set_label_position("right")
ax1.yaxis.set_ticks_position('right')
fig.savefig(f"rake_fault_{M}.png", dpi=600, bbox_inches='tight', pad_inches=0.05)


## To plot psr distribution from the srf file,
## should be consistent with that from sliprate.bin
#fig = plt.figure(figsize=(9, 5))
#ax0 = fig.add_axes([0.1, 0.1, 0.7, 0.8])
#cax = ax0.imshow(psr, extent=[0, nx * dx, nz * dx, 0], aspect = 5,
#                cmap=cm.hot_r, vmax=12)
##cnt = np.arange(0., 100., 1.)
#ct = ax0.contour(grid[1,:,:]*dx, grid[0,:,:]*dx, rt.reshape(nz, nx), 10)#, cnt)
#norm = matplotlib.colors.Normalize(vmin=ct.cvalues.min(), vmax=ct.cvalues.max())
#sm = plt.cm.ScalarMappable(norm=norm, cmap=ct.cmap)
#sm.set_array([])
#cb = plt.colorbar(sm, ticks=ct.levels, shrink=0.6, extend='both', orientation='horizontal')
#cb.set_label("Rupture start time (s)")
#ax0.set_xlabel("Strike (km)")
#xticks = ax0.get_xticks()
#ax0.set_xticks(xticks[:-1])
#ax0.set_ylabel("Depth (km)")
##ax.set_title("fault peak slip rate (PSR)")
# #   PCM=ax.get_children()[2] #get the mappable, the 1st and the 2nd are the x and y axes
#cb = plt.colorbar(cax, ax=ax0, orientation='horizontal', shrink = 0.6, 
#                  aspect = 18, pad=0.15)
#cb.set_label("PSR (m/s)")
#ax1 = fig.add_axes([0.8, 0.507, 0.12, 0.38])
#ax1.plot(np.mean(psr, axis=1), np.arange(nz) * dx)
#ax1.invert_yaxis()
#ax1.set_xlabel("Average PSR (m/s)")
##ax1.set_labelpad = 8
#ax1.set_ylabel("Depth (km)")
##ax1.set_yticks([])
#ax1.yaxis.set_label_position("right")
#ax1.yaxis.set_ticks_position('right')
##ax1.tick_params(axis='y', which='both', labelleft='off', labelright='on')
#
#plt.savefig(f"plot_psr_gp_103018_{module}.svg", dpi=600, bbox_inches='tight', pad_inches=0.05)      
