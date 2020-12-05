#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
from pylab import *
import numpy as np
from scipy.interpolate import NearestNDInterpolator
from sys import stdout
import struct

nx = 9000
ny = 6750
nx_sm = 1400
ny_sm = 1400
npt = nx * ny
npt_sm = nx_sm * ny_sm

stdout.write("Reading grid points ...")
stdout.flush()
ll=np.fromfile('../cvm/la_habra_large_cvmsi_20m.grid').reshape(npt,3)
ll_sm = np.fromfile('../cvm/la_habra_small_cvmsi_20m.grid').reshape(npt_sm, 3)
ll_sm = reshape(ll_sm[:,0:2], (ny_sm,nx_sm,2))
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

nsrc = 1560
nt = 5000
fidin = open("momrate_LaHabra.zf.100m.rev.bin", 'rb')
fidout = open("momrate_LaHabra.zf.100m.rev.large.bin", 'wb')
for ii in range(nsrc):
    stdout.write("\rInterpolating source {0} of {1}... ".format(ii+1, nsrc))
    stdout.flush()
    src_ll = struct.unpack('3I', fidin.read(12))
    idx = src_ll[0]
    idy = src_ll[1]
    lon = ll_sm[idy-1, idx-1, 0]
    lat = ll_sm[idy-1, idx-1, 1]
    xi = IMi((lon, lat))
    yi = IMj((lon, lat))
    src_tmp = fidin.read(4*6*nt)
    fidout.write(struct.pack('3I',xi,yi,src_ll[2]))
    fidout.write(src_tmp)

fidin.close()
fidout.close()
stdout.write("Done.\n")



