#!/ccs/proj/geo112/hzfmer/summit/opt/anaconda3/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue March  19 18:37:05 2019

@author: zhh076
"""

import numpy as np
from numpy import log, exp
from mpi4py import MPI
import sys
import os
sys.path.append("/ccs/home/hzfmer/scratch/LA/gp_rupture_test/LA/gp_rupture_test/scripts")
import BSSA_2014_nga
           

comm = MPI.COMM_WORLD
ncpus = comm.Get_size()
rank = comm.Get_rank()
if rank == 0:
    print(f"ncpus: {ncpus}\nrank: {rank}")
home_dir = "/ccs/home/hzfmer/scratch/LA/gp_rupture_test/LA/gp_rupture_test/"
nx = 6320
ny = 4200
n = nx * ny
# for M in [6.35, 7.35, 8.45]:
for M in [6.45, 8.15]:
    for T in [1, 2, 3, 5]:
        data_gmpe = np.fromfile(f'mesh_info_for_GMPE_M{M}.bin',
                dtype='float32').reshape((-1, 3))
        dist = data_gmpe[:, 0]
        z1 = data_gmpe[:, 1]
        vs30 = data_gmpe[:, 2]
        sa_median_gmpe = np.zeros((ny, nx), dtype='float32')
        sa_sigma_gmpe = np.zeros((ny, nx), dtype='float32')
        for j in range(rank, ny, ncpus):
            for i in range(nx):
                k = j * nx + i 
                sa_median_gmpe[j, i], sa_sigma_gmpe[j, i], _ = BSSA_2014_nga.BSSA_2014_nga(  \
                            M, T, dist[k], 1, 1, z1[k], vs30[k])
            if rank == 0:
                print(f"\rNow computing M={M}, T={T}", end=' ', flush=True)
                print(f"\rNow computing row {j}/{ny}", end=' ', flush=True)
                for p in range(1, ncpus):
                    comm.Recv(sa_median_gmpe[j + p, :], source=p, tag=101)
                    comm.Recv(sa_sigma_gmpe[j + p, :], source=p, tag=102)
            else:
                comm.Send(sa_median_gmpe[j, :], dest=0, tag=101)
                comm.Send(sa_sigma_gmpe[j, :], dest=0, tag=102)
            comm.Barrier()
        if rank == 0:
            print(dist[0])
            output = np.column_stack((dist, np.ravel(sa_median_gmpe), np.ravel(sa_sigma_gmpe)))
            output.tofile(f'sa_GMPE_M{M}_T{T}.bin')
        comm.Barrier()

