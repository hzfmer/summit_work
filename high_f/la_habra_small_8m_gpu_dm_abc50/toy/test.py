import numpy as np

mx, my, mz = 1080, 1080, 164
nx, ny, nz =  360,  360, 60
px, py, pz =  120,  120, 104
mx, my, mz = 3456, 3456, 108
nx, ny, nz = 1152, 1152, 296
px, py, pz =  120,  120, 104


vp = 4000
vs = 1700
rho = 2500
vp = 3300
vs = 1600
rho = 2500
qs = 0.1 * vs
qp = 2 * qs

vp1 = vp * np.ones((my, mx), dtype='float32')
vs1 = vs * np.ones((my, mx), dtype='float32')
rho1 = rho * np.ones((my, mx), dtype='float32')
qp1 = qp * np.ones((my, mx), dtype='float32')
qs1 = qs * np.ones((my, mx), dtype='float32')

with open('../test0_bin', 'wb') as fout:
    for i in range(mz):
        mesh = np.dstack((vp1, vs1, rho1))
        mesh.astype('float32').tofile(fout)

vp2 = vp * np.ones((ny, nx), dtype='float32')
vs2 = vs * np.ones((ny, nx), dtype='float32')
rho2 = rho * np.ones((ny, nx), dtype='float32')
qp2 = qp * np.ones((ny, nx), dtype='float32')
qs2 = qs * np.ones((ny, nx), dtype='float32')
with open('../test1_bin', 'wb') as fout:
    for i in range(nz):
        mesh = np.dstack((vp2, vs2, rho2))
        mesh.astype('float32').tofile(fout)

vp3 = vp * np.ones((py, px), dtype='float32')
vs3 = vs * np.ones((py, px), dtype='float32')
rho3 = rho * np.ones((py, px), dtype='float32')
qp3 = qp * np.ones((py, px), dtype='float32')
qs3 = qs * np.ones((py, px), dtype='float32')
with open('../test2_bin', 'wb') as fout:
    for i in range(pz):
        mesh = np.dstack((vp3, vs3, rho3))
        mesh.astype('float32').tofile(fout)
