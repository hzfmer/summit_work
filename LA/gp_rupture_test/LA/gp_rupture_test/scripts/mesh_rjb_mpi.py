#!/ccs/proj/geo112/hzfmer/summit/opt/anaconda3/bin/python

import numpy as np
from numpy import sin, cos, arcsin, sqrt
from mpi4py import MPI
import os 
import re

def distance(lon1, lat1, lon2, lat2):
    lat1, lon1, lat2, lon2 = np.radians((lat1, lon1, lat2, lon2))
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    d = 0.5 - cos(dlat) / 2 + cos(lat1) * cos(lat2)  * (1 - cos(dlon)) / 2
    return 12742 * arcsin(sqrt(d))

def mindist(npf, lon, lat, tlon, tlat):
    md = 1e13 * np.ones((len(tlon), ))
    for j in range(len(tlon)):
        if j%100 == 0:
            print(f'\r                , Now processing row {j} / {len(tlon)}',
                    end=' ', flush=True)
        for i in range(npf):
            adist = distance(lon[i], lat[i], tlon[j], tlat[j])
            md[j] = min(md[j], adist)
    return md

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()

params = {'6.35': (354, 212), '6.45': (144, 230), '7.35': (980, 240), '8.15': (5116, 220), '8.45': (5116, 220)}
M = re.findall(r"\d+\.\d+", os.getcwd().split('/')[-1])[0]
npf = params[M][0]
nx = 6320
ny = 4200
npt = nx * ny
fdist = np.zeros((ny, nx), dtype='float32')
ll = np.genfromtxt('fault_surf_loc.txt')
tll = np.fromfile('surf.grid', dtype='float32').reshape((ny, nx, 2))
tlon = np.float32(tll[:,:,0])
tlat = np.float32(tll[:,:,1])

for j in range(rank, ny, nprocs):
    fdist[j, :] = mindist(npf, ll[:, 0], ll[:, 1], tlon[j, :], tlat[j, :])
    if rank == 0:
        print(f'\rRow {j} / {ny}', end=' ', flush=True)
        for p in range(1, nprocs):
            comm.Recv(fdist[j+p, :], source=p, tag=100)
    else:
        comm.Send(fdist[j, :], dest=0, tag=100)
    comm.Barrier()

if rank == 0:
    fdist.tofile(f'mesh_rjb_M{M}.bin')
