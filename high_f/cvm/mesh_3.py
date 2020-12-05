#!/usr/bin/env python
'''
Downsample the mesh for AWP-DM, using 3 layers.
'''
import matplotlib
matplotlib.use('Agg')
from pylab import *
import numpy as np
import struct
from sys import stdout

nx = 9000
ny = 6750
nz = 3072

fidin = open('la_habra_large_cvmsi_20m.media','rb')
fidout = open('la_habra_large_cvmsi_20m_3g.media','wb')

for ii in range(58):
    buf = fidin.read(12 * nx * ny * 6)
    fidout.write(buf)
    stdout.write("\rcopying finest layers, k = {0} of {1} ... ".format(ii, 29))
    stdout.flush()

#fidin.seek(12*nx*ny*352,0)
xnum = 3 * nx
xbyte = 4 * xnum
k = 348
while k<nz:
# 348 + 380 + 176
#    if k<352:
#        k += 1
#        jump = 0
#    elif k<1486:
    if k<1488:
        k += 3
        jump = 2
    elif k<nz:
        k += 9
        jump = 8   

#    for j in range(0, ny, jump+1):
#        for i in range(0, nx, jump+1):
#            buf = fidin.read(12)
#            fidin.seek(4 * 3 * jump, 1)
#            fidout.write(buf)
    fidin.seek(4 * 3 * jump * nx * ny, 1)
    xread = xnum / (jump + 1)
    xidx = arange(0, nx, jump+1)
    for j in range(0, ny, jump+1):
        stdout.write("\rDownsampling ({0}/{1}), k = {2} of {3} ...".format(j, ny, k, nz))
        buf = fidin.read(xbyte)
        buf = array(struct.unpack('{0}f'.format(xnum),buf)).reshape(nx,3)
        buf = struct.pack('{0}f'.format(xread),*ravel(buf[xidx,:]))
        fidout.write(buf)
        fidin.seek(jump * xbyte, 1)
            
fidin.close()
fidout.close()
