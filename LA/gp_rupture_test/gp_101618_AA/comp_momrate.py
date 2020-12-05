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
plt.tight_layout()

def compmt(stk, dip, rake, area, mu):
    # stk, dip, rake: degree
    # area: cm^2
    # rho: g/cm^3
    # vs: cm/s
    stk = (stk -35) * pi / 180
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

dt_ref = 0.006
nt_ref = 20000
dx = 0.1
ypos = 1173
xoff = 20

module = os.getcwd()[-1]
os.chdir(os.path.split(os.path.abspath(os.path.realpath(sys.argv[0])))[0])


#f = open(f"saf_{module}.srf",'r')
f = open("saf.srf",'r')
#f_out = open("momrate.dat", 'wb')
f.readline()
f.readline()
token = f.readline()
nx = int(token.split()[2])
nz = int(token.split()[3])
npt = nx * nz


grid = np.mgrid[1 :  nz + 1, xoff + 1 : nx + xoff + 1]
gridz = np.ravel(grid[0,:,:])
gridx = np.ravel(grid[1,:,:])


f.readline()
f.readline()

mom = np.zeros((npt,), dtype=float)
rt = np.zeros((npt,), dtype=float)
trise = np.zeros((npt,), dtype=float)
mu = np.zeros((npt,), dtype=float)
area = np.zeros((npt,), dtype=float)
tslip = np.zeros((npt,), dtype=float)
psr1 = np.zeros((npt,), dtype=float)
psr2 = np.zeros((npt,), dtype=float)
psr3 = np.zeros((npt,), dtype=float)


for k in range(npt):
    sys.stdout.write("\rreading subfault %d of %d" % (k+1, npt))
    sys.stdout.flush()
    sr = np.zeros((nt_ref,), dtype=float)
    
    nl1 = f.readline().split()
    nl2 = f.readline().split()
    stk = float(nl1[3])
    dip = float(nl1[4])
    rake = float(nl2[0])
    
    area = float(nl1[5]) / 1.e4
    tinit= float(nl1[6])
    dt = float(nl1[7])
    vs = float(nl1[8]) / 1.e2
    rho = float(nl1[9]) * 1.e3
    
    slip1 = float(nl2[1]) / 1.e2
    slip2 = float(nl2[3]) / 1.e2
    slip3 = float(nl2[5]) / 1.e2
    nt1 = int(nl2[2])
    nt2 = int(nl2[4])
    nt3 = int(nl2[6])
    nskip1 = int(np.ceil(nt1 / 6.))
    nskip2 = int(np.ceil(nt2 / 6.))
    nskip3 = int(np.ceil(nt3 / 6.))

    sliprate1 = np.zeros((nt1), dtype=float)
    sliprate2 = np.zeros((nt2), dtype=float)
    sliprate3 = np.zeros((nt3), dtype=float)
    p1 = 0
    p2 = 0
    p3 = 0
    for l in range(nskip1):
        tmp = f.readline()   
        nt = len(tmp.split())
        for n in range(nt):
            sliprate1[p1] = float(tmp.split()[n]) / 1.e2
            p1 += 1
    slip1 = sum(sliprate1) * dt 
    if len(sliprate1) > 0:
        psr1[k] = sliprate1.max()
    
#   for l in range(nskip2):
#        tmp = f.readline()  
#        nt = len(tmp.split())
#        for n in range(nt):
#            sliprate2[p2] = float(tmp.split()[n]) / 1.e2
#            p2 += 1
#    slip2 = sum(sliprate2) * dt
#    if len(sliprate2) > 0:
#        psr2[k] = sliprate2.max()
#     
#    for l in range(nskip3):
#        tmp = f.readline()
#        nt = len(tmp.split())
#        for n in range(nt):
#            sliprate3[p3] = float(tmp.split()[n]) / 1.e2
#            p3 += 1
#    slip3 = sum(sliprate3) * dt 
#    if len(sliprate3) > 0:
#        psr3[k] = sliprate3.max() 
#    
    mu[k] = pow(vs, 2.) * rho
    tslip[k] = sqrt(pow(slip1, 2.) + pow(slip2, 2.) + pow(slip3, 2.))
    mom[k] = tslip[k] * mu[k] * area
    rt[k] = tinit
    
    # how about sliprate2, 3?
#   xx, yy, zz, xz, yz, xy = compmt(stk, dip, rake, area, mu[k])
#    trise[k] = len(sliprate1) * dt
#    t1 = int(np.ceil(tinit/dt_ref))
#    n_sr = int(np.ceil(trise[k] / dt_ref))
#    if n_sr>0:
#        sr[t1 : t1+n_sr] = resample(sliprate1, n_sr)
# Do we need multiply dt_ref?
#    momrate = np.ravel(np.array([sr * x for x in (xx, yy, zz, xz, yz, xy)]).T)
#    f_out.write(pack('3i', gridx[k], ypos, gridz[k]))
#    f_out.write(pack(f'{nt_ref * 6}f', *momrate))
    
          
f.close()
#f_out.close()
#mag = 2. / 3. * (np.log10(sum(mom)) - 9.1)
#print(f"Magnitude is {mag}")
#tslip = tslip.reshape((nz, nx))    
psr = psr1.reshape((nz, nx))  


fig = plt.figure(figsize=(9, 5))
ax0 = fig.add_axes([0.1, 0.1, 0.7, 0.8])
cax = ax0.imshow(psr, extent=[0, nx * dx, nz * dx, 0], aspect = 5,
                cmap=cm.hot_r, vmax=12)
#cnt = np.arange(0., 100., 1.)
ct = ax0.contour(grid[1,:,:]*dx, grid[0,:,:]*dx, rt.reshape(nz, nx), 10)#, cnt)
norm = matplotlib.colors.Normalize(vmin=ct.cvalues.min(), vmax=ct.cvalues.max())
sm = plt.cm.ScalarMappable(norm=norm, cmap=ct.cmap)
sm.set_array([])
cb = plt.colorbar(sm, ticks=ct.levels, shrink=0.6, extend='both', orientation='horizontal')
cb.set_label("Rupture start time (s)")
ax0.set_xlabel("Strike (km)")
xticks = ax0.get_xticks()
ax0.set_xticks(xticks[:-1])
ax0.set_ylabel("Depth (km)")
#ax.set_title("fault peak slip rate (PSR)")
 #   PCM=ax.get_children()[2] #get the mappable, the 1st and the 2nd are the x and y axes
cb = plt.colorbar(cax, ax=ax0, orientation='horizontal', shrink = 0.6, 
                  aspect = 18, pad=0.15)
cb.set_label("PSR (m/s)")
ax1 = fig.add_axes([0.8, 0.507, 0.12, 0.38])
ax1.plot(np.mean(psr, axis=1), np.arange(nz) * dx)
ax1.invert_yaxis()
ax1.set_xlabel("Average PSR (m/s)")
#ax1.set_labelpad = 8
ax1.set_ylabel("Depth (km)")
#ax1.set_yticks([])
ax1.yaxis.set_label_position("right")
ax1.yaxis.set_ticks_position('right')
#ax1.tick_params(axis='y', which='both', labelleft='off', labelright='on')

plt.savefig(f"plot_psr_gp_101618_{module}.svg", dpi=600, bbox_inches='tight', pad_inches=0.05)      
