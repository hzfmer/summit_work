import numpy as np
import sys

mx, my, mz = 3456, 3456, 108
nx, ny, nz = 1152, 1152, 464


vp = 3300
vs = 1600
rho = 2500

vp = 3550
vs = vp / np.sqrt(2)
rho = 2500

vp1 = vp * np.ones((my, mx), dtype='float32')
vs1 = vs * np.ones((my, mx), dtype='float32')
rho1 = rho * np.ones((my, mx), dtype='float32')
with open('homo_0_bin', 'wb') as fout:
    for i in range(mz):
        mesh = np.dstack((vp1, vs1, rho1))
        mesh.astype('float32').tofile(fout)

#sys.exit(-1)

vp2 = vp * np.ones((ny, nx), dtype='float32')
vs2 = vs * np.ones((ny, nx), dtype='float32')
rho2 = rho * np.ones((ny, nx), dtype='float32')
with open('homo_1_bin', 'wb') as fout:
    for i in range(nz):
        mesh = np.dstack((vp2, vs2, rho2))
        mesh.astype('float32').tofile(fout)

