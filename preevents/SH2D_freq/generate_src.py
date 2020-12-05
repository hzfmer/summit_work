#!/ccs/home/hzfmer/file_back/programs/anaconda3/bin/python

import numpy as np

nz, ny, nx, nvar = 1000, 100, 100, 5
dx = 2
vs = np.fromfile('mesh_0', dtype='float32').reshape(nz, ny, nx, nvar)[:,:,:,1]

f0 = 60
t0 = 4 / f0
tmax = 3
dt = dx / 2 / vs.max()
dt = 0.0005

t = np.arange(np.ceil(tmax / dt)) * dt
src = -2 * (t - t0) * f0 ** 2 * np.exp(-(f0 ** 2) * (t - t0) ** 2)
vx = np.zeros_like(src)
data = np.dstack((vx, src, vx)).squeeze()

np.savetxt('stf.txt_0', data, fmt="%.10f")
# np.savetxt('stf.txt_0', data)
