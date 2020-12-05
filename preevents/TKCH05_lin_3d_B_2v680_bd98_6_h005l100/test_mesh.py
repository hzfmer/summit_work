#!/usr/env/python

import numpy as np

nx, ny, nz, nvar = 400, 400, 400, 5
mesh = np.fromfile('mesh', dtype='f', count=nvar*nx*ny*nz).reshape(nz,
        ny, nx, nvar)
vs = mesh[:, :, :, 1].swapaxes(0, 2)
print(vs[200, 200, 30:40])
