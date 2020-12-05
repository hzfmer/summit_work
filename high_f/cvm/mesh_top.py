#!/usr/bin/env python

import numpy as np
import sys

#fid_in = open('la_habra_large_cvmsi_het_top_20m.media','rb')
#fidout = open('la_habra_large_cvmsi_het_20m_p0.media','wb')
fid_in = open('la_habra_large_cvmsi_20m.media','rb')
fidout = open('la_habra_large_cvmsi_20m_p0.media','wb')


nx = 9000
ny = 6750
ntop = 380
nf = 4 # size of float
nvar = 3

for ii in range(ntop//20):
    buf = fid_in.read(20*nvar*nf*nx*ny)
    fidout.write(buf)
    sys.stdout.write("\rWriting the {0}/{1} layer.".format(20*ii, ntop))
    sys.stdout.flush()

fid_in.close()
fidout.close()
