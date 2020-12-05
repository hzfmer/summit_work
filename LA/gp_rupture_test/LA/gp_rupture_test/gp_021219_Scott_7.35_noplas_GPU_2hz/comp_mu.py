#!/usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os

module = os.getcwd().split('/')[-1]

nx = 3000
nz = 180
xoff = 20
dx = 0.1
grid = np.mgrid[1 :  nz + 1, xoff + 1 : nx - xoff + 1]
gridz = np.ravel(grid[0,:,:])
gridx = np.ravel(grid[1,:,:])
ypos = 1173

fl=np.loadtxt("/lustre/atlas/proj-shared/geo112/huzf/LA/source_7.8_so_cvm_lvz/fault_3000.dat")
rho_3d=fl[:,5].reshape(nz,nx)
vs_3d=fl[:,4].reshape(nz,nx)
mu_3d = rho_3d * vs_3d * vs_3d
mu_3d = mu_3d[:, xoff + 1 : nx - xoff + 1]
vs_3d = vs_3d[:, xoff + 1 : nx - xoff + 1]

f = open("saf.srf",'r')
f_sliprate = open("sliprate.bin", 'wb')
f_fault = open("subfaults.idx", 'w')
f.readline()
f.readline()
token = f.readline()
nx = int(token.split()[2])
nz = int(token.split()[3])
npt = nx * nz
mu = np.zeros((nz, nx), dtype='float')
tinit = np.zeros((nz, nx), dtype='float')
vrup = np.zeros((nz, nx), dtype='float')
vs = np.zeros((nz, nx), dtype='float')
f.readline()
f.readline()

for k in range(npt):
    xid = k % nx
    zid = k // nx
    nl1 = f.readline().split()
    nl2 = f.readline().split()
    stk = float(nl1[3]) - 35 # so that the fault is along y direction, at 90'
    dip = float(nl1[4])
    rake = float(nl2[0])
    area = float(nl1[5]) / 1.e4
    tinit[zid, xid] = float(nl1[6])
    dt = float(nl1[7])
    vs[zid, xid] = float(nl1[8]) / 1.e2
    rho = float(nl1[9]) * 1.e3
    nt1 = int(nl2[2])
    nskip1 = int(np.ceil(nt1 / 6.))
    sliprate1 = np.zeros((nt1), dtype=float)
    p1 = 0
    for l in range(nskip1):
        tmp = f.readline()
        nt = len(tmp.split())
        for n in range(nt):
            sliprate1[p1] = float(tmp.split()[n]) / 1.e2
            p1 += 1
    mu[zid, xid] = pow(vs[zid, xid], 2.) * rho
    trise = len(sliprate1) * dt
    t1 = int(np.floor(tinit[zid, xid]/dt))

    f_fault.write(f"{gridx[k]} {ypos} {gridz[k]} {area} "  \
            + f"{mu_3d[zid, xid]} {stk} {dip} {rake}\n")

f_fault.close()

vrup[:, 1:] = np.diff(tinit)

fig, ax = plt.subplots(2, 1, figsize=(8,8), sharex=True)
ax[0].imshow(mu, extent=[0, nx * dx, nz * dx, 0], aspect=5, cmap=cm.hot_r)
ax[0].set_ylabel("depth (km)")

ax[1].imshow(mu_3d, extent=[0, nx * dx, nz * dx, 0], aspect=5, cmap=cm.hot_r)
ax[1].set_xlabel("Strike (km)")
xticks = ax[1].get_xticks()
ax[1].set_xticks(xticks[:-1])
ax[1].set_ylabel("Depth (km)")

plt.savefig(f"plot_mu_{module}.svg", dpi=600, bbox_inches='tight', pad_inches=0.05)

fig, ax = plt.subplots(2, 1, figsize=(8,8), sharex=True)
ax[0].imshow(vrup / vs, extent=[0, nx * dx, nz * dx, 0], aspect=5, cmap=cm.hot_r)
ax[0].set_ylabel("depth (km)")

ax[1].imshow(vrup / vs_3d, extent=[0, nx * dx, nz * dx, 0], aspect=5, cmap=cm.hot_r)
ax[1].set_xlabel("Strike (km)")
xticks = ax[1].get_xticks()
ax[1].set_xticks(xticks[:-1])
ax[1].set_ylabel("Depth (km)")

plt.savefig(f"plot_vrup_{module}.svg", dpi=600, bbox_inches='tight', pad_inches=0.05)
