#! /usr/bin/env python

import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from sys import argv, exit

plt.tight_layout()

if len(argv) != 5:
    print("usage: %s mesh nx ny nz" % argv[0])
    exit()

fname=argv[1]
nx = 9000
ny = 6750
nz = 3072
nx = int(argv[2])
ny = int(argv[3])
nz = int(argv[4])
nvar = 3
nf = 4  # float number bytes
f = open(fname,'rb')
for y in range(100,6750,1000):
    v = np.zeros((nx,nz))
    for z in range(nz):
        f.seek((z*ny*nx + y*nx) * nvar * nf, 0)
        buf = np.frombuffer(f.read(nx * nvar * nf), dtype='f')
        v[:,z] = buf[1::3]
    fig, ax = plt.subplots(figsize=(8,6))
    cax=ax.imshow(v.T, vmin=v.min(), vmax=v.max())
    cb = plt.colorbar(cax, ax=ax, orientation="horizontal")
    cb.set_label('Vs (m/s')
    fig.savefig("slice_%d.png" % y, dpi=600, bbox_inches='tight', pad_inches=0.1)

f.close()
              
