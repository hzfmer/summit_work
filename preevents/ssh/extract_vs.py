#!/usr/bin/env python

import numpy as np

nx, ny, nz, nvar = 400, 400, 400, 5

mesh_0 = "mesh_TKCH05_het_052119_2v680_bd98_h005l100.bin"
mesh = np.fromfile(mesh_0, dtype='f', count=nvar*nx*ny*nz).reshape(nz,
        ny, nx, nvar)
vs = mesh[:,200,:,1]
vp = mesh[:,200,:,0]
rho = mesh[:,200,:,2]

vs.tofile('vs_west_east.bin')
vp.tofile('vp_west_east.bin')
rho.tofile('rho_west_east.bin')




