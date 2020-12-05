#!/ccs/home/hzfmer/file_back/programs/anaconda3/bin/python

import numpy as np

def estimateRhofromVs(vs):
    N = len(vs)
    rho = np.zeros(N)
    for i in range(N):
        if vs[i] <= 333:
            rho[i] = 2000 + 190 * (vs[i] - 166.5) / 333
        elif vs[i] < 560:
            rho[i] = 2350
        elif vs[i] < 640:
            rho[i] = 2450
        else:
            rho[i] = 2750
    return rho

def estimateQfromVs(vs):
    N = len(vs)
    q = np.zeros(N)
    for i in range(N):
        if vs[i] <= 1000:
            q[i] = 0.06 * vs[i]
        elif vs[i] < 2000:
            q[i] = 0.04 * vs[i]
        else:
            q[i] = 0.16 * vs[i]
    q = 1e10 * np.ones_like(vs)
    return q

nx, ny, nz = 100, 100, 1000
z = [20]
z.append(nz - sum(z))
vs = [400.0, 1100.0]
rho = estimateRhofromVs(vs)
vp = vs
qs = estimateQfromVs(vs)
qp = 2 * qs

vp = np.vstack([i * np.ones((j, ny, nx)) for (i, j) in zip(vp, z)])
vs = np.vstack([i * np.ones((j, ny, nx)) for (i, j) in zip(vs, z)])
rho = np.vstack([i * np.ones((j, ny, nx)) for (i, j) in zip(rho, z)])
qp = np.vstack([i * np.ones((j, ny, nx)) for (i, j) in zip(qp, z)])
qs = np.vstack([i * np.ones((j, ny, nx)) for (i, j) in zip(qs, z)])

data = np.dstack([x.flatten() for x in (vp, vs, rho, qp, qs)]).squeeze().astype('float32')
# data = np.dstack([x.flatten() for x in (vp, vs, rho)]).squeeze().astype('float32')
data.tofile('mesh_0')
