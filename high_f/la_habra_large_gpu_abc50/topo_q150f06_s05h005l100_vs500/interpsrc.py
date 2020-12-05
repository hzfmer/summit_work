#!/usr/bin/env python3

"""Convert medium domain source indices to large domain
(1) Fetch medium domain grid, containing lon/lat
(2) convert indices of subfault to lon/lat in medium domain
(3) fetch large domain grid
(4) locate the lon/lat in the large domain grid using interpolation
"""

import matplotlib
matplotlib.use('Agg')
from pylab import *
import numpy as np
from scipy.interpolate import NearestNDInterpolator
from sys import stdout
import struct

# hard coded size of domains
nx = 9504
ny = 7020
nx_med = 5940
ny_med= 4860
npt = nx * ny
npt_med = nx_med * ny_med

stdout.write("Reading grid points ...")
stdout.flush()
ll=np.fromfile('../cvm/la_habra_ext_large_cvmsi_20m.grid').reshape(npt,3)
ll_med = np.fromfile('../cvm/la_habra_medium_cvmsi_ijk12_8m_5940.grid').reshape(npt_med, 3)
ll_med = reshape(ll_med[:,0:2], (ny_med,nx_med,2))
stdout.write("ok.\n")

alli = zeros((npt,), dtype=float)
allj = zeros((npt,), dtype=float)

n = 0
for k in range(ny):
    stdout.write("\rDefining interpolation data ({0}/{1})".format(k+1, ny))
    stdout.flush()
    for l in range(nx):
        alli[n] = l + 1
        allj[n] = k + 1
        n += 1

stdout.write("Done.\nSetting up Interpolators ...")
stdout.flush()

IMi = NearestNDInterpolator(ll[:,0:2], alli)
IMj = NearestNDInterpolator(ll[:,0:2], allj)

nsrc = 15625
nt = 5000
fidin = open("../la_habra_small_8m_gpu_dm_abc50/momrate.dat", 'rb')
fidout = open("momrate_LaHabra.20m.large.bin", 'wb')
for ii in range(nsrc):
    stdout.write("\rInterpolating source {0} of {1}... ".format(ii+1, nsrc))
    stdout.flush()
    src_ll = struct.unpack('3I', fidin.read(12))
    idx = src_ll[0]
    idy = src_ll[1]
    lon = ll_med[idy-1, idx-1, 0]
    lat = ll_med[idy-1, idx-1, 1]
    xi = IMi((lon, lat))
    yi = IMj((lon, lat))
    src_tmp = fidin.read(4*6*nt)
    fidout.write(struct.pack('3I',xi,yi,src_ll[2]))
    fidout.write(src_tmp)

fidin.close()
fidout.close()
stdout.write("Done.\n")



