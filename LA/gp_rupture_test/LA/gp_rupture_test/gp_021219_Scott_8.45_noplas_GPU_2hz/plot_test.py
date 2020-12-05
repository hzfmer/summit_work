#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 17 2018

@author: Zhifeng Hu <zhh076@ucsd.edu>
"""

import numpy as np
from numpy import sin, cos, pi, sqrt
from scipy.signal import resample
from struct import pack
import os
import sys
import numpy as np
from numpy import sin, cos, pi, sqrt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
plt.tight_layout()

def upsample(s1, dt1, dt2, nt2):
    '''
    For upsample only, since no anti-aliasing is implemented.
    '''
    s2 = np.zeros((nt2,))
    for n in range(nt2):
        t_n = dt2 * n
        n_l = np.floor(t_n / dt1)
        n_r = np.ceil(t_n / dt1)
        if n_l==n_r:
            s2[n] = s1[n_l] 
        elif n_r > len(s1) - 1:
            s2[n] = s1[n_l]
        else:
            t_l = n_l * dt1
            t_r = n_r * dt1
            s2[n] = s1[n_l] + (s1[n_r] - s1[n_l]) * (t_n - t_l) / (t_r - t_l)
    return s2

module = os.getcwd()[-1]
nt_ref = 1200
dx = 0.1
ypos = 1173
xoff = 20

f = open("saf.srf",'r')
#f_sliprate = open("sliprate.bin", 'wb')
#f_fault = open("subfaults.idx", 'w')
f.readline()
f.readline()
token = f.readline()
nx = int(token.split()[2])
nz = int(token.split()[3])
npt = nx * nz
final_sr = np.zeros((nx, nz, nt_ref), dtype='float32')
grid = np.mgrid[1 :  nz + 1, xoff + 1 : nx + xoff + 1]
gridz = np.ravel(grid[0,:,:])
gridx = np.ravel(grid[1,:,:])


f.readline()
f.readline()

strike = np.zeros((nz, nx))
dip = np.zeros((nz, nx))
rake = np.zeros((nz, nx))
for k in range(npt):
    sys.stdout.write("\rreading subfault %d of %d" % (k+1, npt))
    sys.stdout.flush()
    sr = np.zeros((nt_ref,))
    nl1 = f.readline().split()
    nl2 = f.readline().split()
    strike[k//nx, k%nx] = float(nl1[3]) - 35 # so that the fault is along y direction, at 90'
    dip[k//nx, k%nx] = float(nl1[4])
    rake[k//nx, k%nx] = float(nl2[0])
    area = float(nl1[5]) / 1.e4
    tinit= float(nl1[6])
    dt = float(nl1[7])
    vs = float(nl1[8]) / 1.e2
    rho = float(nl1[9]) * 1.e3
    nt1 = int(nl2[2])
    nskip1 = int(np.ceil(nt1 / 6.))
#    sliprate1 = np.zeros((nt1), dtype=float)
    p1 = 0
    for l in range(nskip1):
        tmp = f.readline()
#        nt = len(tmp.split())
#        for n in range(nt):
#            sliprate1[p1] = float(tmp.split()[n]) / 1.e2
#            p1 += 1
#    mu = pow(vs, 2.) * rho
#    trise = len(sliprate1) * dt
#    t1 = int(np.floor(tinit/dt))
#    sr[t1 : t1+len(sliprate1)] = sliprate1
#      
#n_sr = int(np.ceil(trise / dt_ref))
#if n_sr>0:
#sr[t1 : t1+n_sr] = resample(sliprate1, n_sr)
#sr[t1 : t1+n_sr] = upsample(sliprate1, dt, dt_ref, n_sr)
fig = plt.figure(figsize=(9, 5))
ax0 = fig.add_axes([0.1, 0.1, 0.7, 0.8])
cax = ax0.imshow(strike, extent=[0, nx * dx, nz * dx, 0], aspect = 5,
                cmap=cm.hot_r)
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
cb.set_label("Strike")
ax1 = fig.add_axes([0.8, 0.34, 0.12, 0.38])
ax1.plot(np.mean(strike, axis=1), np.arange(nz) * dx)
ax1.invert_yaxis()
ax1.set_xlabel("Average Strike angle")
#ax1.set_labelpad = 8
ax1.set_ylabel("Depth (km)")
#ax1.set_yticks([])
ax1.yaxis.set_label_position("right")
ax1.yaxis.set_ticks_position('right')
#ax1.tick_params(axis='y', which='both', labelleft='off', labelright='on')

plt.savefig(f"plot_strike_gp_101618_{module}.svg", dpi=600, bbox_inches='tight', pad_inches=0.05)

fig = plt.figure(figsize=(9, 5))
ax0 = fig.add_axes([0.1, 0.1, 0.7, 0.8])
cax = ax0.imshow(dip, extent=[0, nx * dx, nz * dx, 0], aspect = 5,
                cmap=cm.hot_r)
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
cb.set_label("Dip")
ax1 = fig.add_axes([0.8, 0.34, 0.12, 0.38])
ax1.plot(np.mean(dip, axis=1), np.arange(nz) * dx)
ax1.invert_yaxis()
ax1.set_xlabel("Average dip angle")
#ax1.set_labelpad = 8
ax1.set_ylabel("Depth (km)")
#ax1.set_yticks([])
ax1.yaxis.set_label_position("right")
ax1.yaxis.set_ticks_position('right')
#ax1.tick_params(axis='y', which='both', labelleft='off', labelright='on')

plt.savefig(f"plot_dip_gp_101618_{module}.svg", dpi=600, bbox_inches='tight', pad_inches=0.05)

fig = plt.figure(figsize=(9, 5))
ax0 = fig.add_axes([0.1, 0.1, 0.7, 0.8])
cax = ax0.imshow(rake, extent=[0, nx * dx, nz * dx, 0], aspect = 5,
                cmap=cm.hot_r)
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
cb.set_label("Rake")
ax1 = fig.add_axes([0.8, 0.34, 0.12, 0.38])
ax1.plot(np.mean(rake, axis=1), np.arange(nz) * dx)
ax1.invert_yaxis()
ax1.set_xlabel("Average rake angle")
#ax1.set_labelpad = 8
ax1.set_ylabel("Depth (km)")
#ax1.set_yticks([])
ax1.yaxis.set_label_position("right")
ax1.yaxis.set_ticks_position('right')
#ax1.tick_params(axis='y', which='both', labelleft='off', labelright='on')

plt.savefig(f"plot_rake_gp_101618_{module}.svg", dpi=600, bbox_inches='tight', pad_inches=0.05)
