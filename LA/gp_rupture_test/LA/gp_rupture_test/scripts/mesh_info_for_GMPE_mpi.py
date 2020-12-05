#!/ccs/proj/geo112/hzfmer/summit/opt/anaconda3/bin/python

import numpy as np
from mpi4py import MPI
import os

def comp_z1(vs_profile):
    if vs_profile.min() > 1000:
        return 0.
    else:
        return np.argwhere(vs_profile <= 1000)[-1][0] * dh


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
ncpus = comm.Get_size()
print(f"ncpus: {ncpus}\nrank: {rank}")

# vs30 just use vs[0] + vs[1]?
nx = 6320
ny = 4200
nz = 400
dh = 0.1

# order: (nz, ny, nx, 3)
vs = np.fromfile('../../data_for_Zhifeng/mesh/vs.media',
        dtype='float32').reshape((nz, ny, nx))

for M in [6.35, 7.35, 8.45]:
    try:
        dist = np.fromfile(f'mesh_rjb_M{M}.bin', dtype='float32')
    except:
        print("No mesh_rjb.bin found!\n")
        continue
    z1 = np.zeros((nx * ny, 1), dtype='float32')
    vs30 = np.zeros((nx * ny, 1), dtype='float32')
    print("correct")
    for iy in range(rank, ny, ncpus):
        for ix in range(nx):
            k = iy *nx + ix
            z1[k] = comp_z1(vs[:, iy, ix])
            if vs[0, iy, ix] < 501:
                vs30[k] = 250
            elif vs[0, iy, ix] >= 501:
                vs30[k] = 760
        if rank == 0:
            print(f"\rNow computing row {iy}/{ny}", end=' ', flush=True)
            for p in range(1, ncpus):
                comm.Recv(z1[(p + iy) * nx : (p + iy + 1) * nx], source=p, tag = 102) 
                comm.Recv(vs30[(p + iy) * nx : (p + iy + 1) * nx], source=p, tag = 103) 
        else:
            comm.Send(z1[iy * nx : (iy + 1) * nx], dest=0, tag=102)
            comm.Send(vs30[iy * nx : (iy + 1) * nx], dest=0, tag=103)
        comm.Barrier()
    if rank == 0:
        data = np.column_stack((dist, z1, vs30))
        print(dist[0])
        data.tofile(f'mesh_info_for_GMPE_M{M}.bin', format='float32')
    comm.Barrier()
