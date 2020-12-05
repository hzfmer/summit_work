#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 24 01:08:57 2018

@author: zhh076

plot plastic strain
"""

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import os

params = {'axes.labelsize': 20,
          'axes.titlesize': 20,
          'font.size': 20,
          'legend.fontsize': 16,
          'xtick.labelsize': 20,
          'ytick.labelsize': 20,
          'xtick.major.size': 10,
          'xtick.minor.size': 5,
          'ytick.major.size': 10,
          'ytick.minor.size': 5}
mpl.rcParams.update(params)

homedir = os.getcwd()
basedir = os.path.dirname(homedir)
cs = homedir.split('_')[-1]
rock = basedir.split('/')[-1][13:]
nx = 3000
ny = 1372
nz = 180
dx = 0.1

x = np.linspace(0, (nx-1) * dx, nx) / 1e3
y = np.linspace(0, (ny-1) * dx, ny) / 1e3


epxx = np.fromfile(os.path.join(homedir, "output_ep/epxx"), 
                    'f').reshape(nz, ny, nx) 
epxx = np.abs(np.flipud(epxx[nz-1, :, :]))
#epxx = np.log10(epxx)
epmax = np.ceil(epxx.max())
print(epmax)
print(np.linspace(0, epmax * 1.1, 5))
#cnt = np.linspace(0, epmax, 5)
fig = plt.figure(figsize=(9.6, 6))
ax = plt.subplot(111)
#plt.imshow(epxx, extent=[0, nx * dx, ny * dx, 0], 
#        vmax = 8e-4, cmap=plt.get_cmap('jet'))
plt.imshow(epxx, extent=[0, nx * dx, ny * dx, 0], 
         norm=LogNorm(vmin=1e-8), cmap='RdBu')
#cb = plt.colorbar(orientation = "vertical", shrink=0.7, 
#              ticks=np.linspace(0, 8e-4, 5), aspect=8)
cb = plt.colorbar(orientation = "vertical", shrink=0.7, 
               ticks=[10**x for x in np.linspace(-8,-1,8)], aspect=8)
#cb.set_label("strain")
plt.xlabel("Strike (km)")
plt.ylabel("Y (km)")
plt.title("Plastic strain")
plt.savefig(os.path.join(basedir,"epxx_{0}_{1}.pdf".format(rock, cs)))





