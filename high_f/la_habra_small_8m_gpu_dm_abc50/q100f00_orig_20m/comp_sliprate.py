#!/usr/bin/env python
"""
Created on Sun Feb 24 2020
@author: Zhifeng Hu <zhh076@ucsd.edu>
    Input:
        *.srf (string): srf source file from CyberShake
        mesh (string): mesh file (la_habra_small_cvmsi_8m_3456.media), including vp, vs, rho

    Output: subfaults.idx, fault information
        sliprate.bin, sliprate time history on each subfault
"""

import numpy as np
from numpy import sin, cos, pi, sqrt
from scipy.signal import resample
from scipy import spatial
from struct import pack
import sys
import glob
from filter_BU import filt_B

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

nt_ref = 5000
dz = 0.020
theta_rot = 39.9

# x-axis: azimuth = 129.9; y-axis: azimuth = 39.9
mx, my = 1400, 1400
grids = np.fromfile('surf.grid', dtype='float64').reshape(my, mx, 3)[:, :, :2] # lon, lat
grids = np.reshape(grids, (-1, 2))
kdtree = spatial.cKDTree(grids)
del grids


f = open(glob.glob('./*.srf')[0],'r')
f.readline()
f.readline()
token = f.readline()
nx = int(token.split()[2])
nz = int(token.split()[3])
npt = nx * nz
mz = int(float(token.split()[5]) / dz) + 1
final_sr = np.zeros((nx, nz, nt_ref), dtype='float32') # (nx, nz, nt) --> (nt, nz, nx) to save

dtop = np.floor(float(f.readline().split()[2]) / dz) * dz 
print(f.readline())

f_sliprate = open("sliprate.bin", 'wb')
f_fault = open("subfaults.idx", 'w')
f_media = open('mesh', 'rb')

#mu = np.zeros((mz, my, mx), dtype='float64')
#with open('mesh', 'rb') as f_media:
#    f_media.seek(int(dtop // dz) * 4 * my * mx * 3, 0)
#    for i in range(mz):
#        try:
#            data = np.frombuffer(f_media.read(4 * my * mx * 3), 
#                             dtype='float32').reshape(my, mx, 3)
#        except:
#            print(int(dtop // dz) + i)
#        mu[i, :, :] = data[:, :, 2] * data[:, :, 1] ** 2
#print(mu.shape)

for j in range(nz):
    sys.stdout.write("\rreading subfault %d of %d" % (j+1, nz))
    sys.stdout.flush()
    for i in range(nx):
        sr = np.zeros((nt_ref,))
        nl1 = f.readline().split()
        nl2 = f.readline().split()
        iy, ix = np.unravel_index(kdtree.query([float(nl1[0]), float(nl1[1])])[1], (my, mx))
        iz = int(float(nl1[2]) / dz + 0.5)
        f_media.seek(iz * 4 * my * mx * 3, 0)
        data = np.frombuffer(f_media.read(4 * my * mx * 3), dtype='float32').reshape(my, mx, 3)
        mu = data[iy, ix, 2] * data[iy, ix, 1] ** 2
        stk = float(nl1[3]) - theta_rot
        dip = float(nl1[4])
        rake = float(nl2[0])
        area = float(nl1[5]) / 1.e4  # cm^2 --> m^2
        tinit= float(nl1[6])
        dt = float(nl1[7])
        nt1 = int(nl2[2])
        nskip1 = int(np.ceil(nt1 / 6.))
        sliprate1 = np.zeros((nt1), dtype=float)
        p1 = 0
        for l in range(nskip1):
            tmp = f.readline()
            nt = len(tmp.split())
            for n in range(nt):
                sliprate1[p1] = float(tmp.split()[n]) / 1.e2  # cm/s2 --> m/s2
                p1 += 1
        t1 = int(np.floor(tinit/dt))
        try:
            sr[t1 : t1+len(sliprate1)] = sliprate1.copy()
        except:
            print(nt1, len(sliprate1), t1)
            print(nl2)
            sys.exit(-1)
        f_fault.write("{0} {1} {2} {3} {4} {5} {6} {7}\n".format(
                    ix, iy, iz, area, mu, stk, dip, rake))
        final_sr[i, j, :] = sr
        # final_sr[i, j, :] = filt_B(sr, 1/dt, highcut=10, causal=False)
    # print(k)
        
final_sr = np.float32(final_sr.T)  # shape = (nt, nz, nx)
final_sr.tofile(f_sliprate, format='float32')
f.close()
f_fault.close()
f_sliprate.close()
f_media.close()

