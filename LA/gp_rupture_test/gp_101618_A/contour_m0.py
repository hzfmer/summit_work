#! /usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
from sys import argv, exit

if len(argv)!=2:
    print("Usage: %s m0.dat" % (argv[0]))
    exit(-1)

module = os.getcwd().split('/')[-1]
nx = 2960
nz = 180
xoff = 20
dx = 0.1
npt = nx * nz
m0 = np.zeros((nz, nx))
with open (f"{argv[1]}", 'r') as f:
    for line in f:
        dat = line.strip('\n').split()
        xi = int(dat[0]) - 1 - xoff
        zi = int(dat[2]) - 1
        m0_tmp = float(dat[3])
        m0[zi, xi] = m0_tmp
print("Total M0 is: %s\n" % sum(sum(m0)))
print("Magnitude is: ",  2. / 3. * (np.log10(sum(sum(m0))) - 9.1))
fig, ax = plt.subplots(figsize=(9,5))
cax = ax.imshow(m0, extent=[0, nx * dx, nz * dx, 0], aspect = 5,
                        cmap=plt.cm.hot_r)
ax.set_xlabel("Strike (km)")
ax.set_ylabel("Depth (km)")
cb = plt.colorbar(cax, ax=ax, orientation='horizontal', shrink = 0.6,
                          aspect = 18, pad=0.15)
cb.set_label("M0")
plt.savefig(f"plot_m0_{module}.svg", dpi=600, bbox_inches='tight', pad_inches=0.05)
