#!/usr/bin/env python

import numpy as np
import sys

fid_hom = open('la_habra_large_cvmsi_20m.media','rb')
fid_het = open('la_habra_large_cvmsi_het_top_20m.media','rb')
fidout = open('la_habra_large_cvmsi_het_20m.media','wb')

nx = 9000
ny = 6750
nz = 3072
ntop = 500
nf = 4 # size of float
nvar = 3

for ii in range(ntop//25):
    buf = fid_het.read(25*nvar*nf*nx*ny)
    fidout.write(buf)
    sys.stdout.write("\rWriting the {0}/{1} layer.".format(25*ii, nz))
    sys.stdout.flush()

fid_hom.seek(ntop*nf*nvar*nx*ny)
buf = fid_hom.read(4*nvar*nf*nx*ny)
fidout.write(buf)

for ii in range(107):
    buf = fid_hom.read(24*nvar*nf*nx*ny)
    fidout.write(buf)
    sys.stdout.write("\rWriting the {0}/{1} layer.".format(504+24*ii, nz))
    sys.stdout.flush()
fid_hom.close()
fid_het.close()
fidout.close()
