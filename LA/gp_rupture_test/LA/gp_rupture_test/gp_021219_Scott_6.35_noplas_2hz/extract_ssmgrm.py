#!/ccs/proj/geo112/hzfmer/summit/opt/anaconda3/bin/python

import numpy as np
import os
from filter_BU import filt_B
from read_params import read_params
from mpi4py import MPI


comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


args = read_params()

ntskp, wstep = args['NTISKP'], args['WRITE_STEP']
skip = wstep * ntskp
dt = args['DT']
nt = np.int(args['TMAX'] / dt)
t = np.linspace(0, nt * dt, nt // ntskp)
nx, ny, nz = args['X'], args['Y'], args['Z']
base_x = args['SXRGO'].strip('\'')
base_y = args['SYRGO'].strip('\'')
base_z = args['SZRGO'].strip('\'')
sites = np.genfromtxt('stat.txt', dtype='int', skip_header=1)
print(base_x, type(base_x))

if rank == 0:
    vx = np.zeros((len(sites), nt // ntskp))
    vy = np.zeros((len(sites), nt // ntskp))
    vz = np.zeros((len(sites), nt // ntskp))
else:
    vx = np.empty((len(sites), nt // ntskp))
    vy = np.empty((len(sites), nt // ntskp))
    vz = np.empty((len(sites), nt // ntskp))
comm.Bcast(vx, root=0)
comm.Bcast(vy, root=0)
comm.Bcast(vz, root=0)

for i in range(0, nt // skip):
    idx = np.arange(i * wstep, (i + 1) * wstep)
    if rank == 0:
        print(f"Reading timestep at {(i + 1) * skip}")
        print("filename", base_x + f'{(i + 1) * skip:07d}')
    vx_tmp = np.fromfile(base_x + f'{(i + 1) * skip:07d}', dtype='f').reshape(-1, ny, nx)
    vy_tmp = np.fromfile(base_y + f'{(i + 1) * skip:07d}', dtype='f').reshape(-1, ny, nx)
    vz_tmp = np.fromfile(base_z + f'{(i + 1) * skip:07d}', dtype='f').reshape(-1, ny, nx)
    for j in range(rank, len(sites), size):
        print(f"reading site {j} / {len(sites)} at time {i * skip}")
        sx, sy, _ = sites[j, :]
        vx[j, idx] = vx_tmp[:, sy - 1, sx - 1]
        vy[j, idx] = vy_tmp[:, sy - 1, sx - 1]
        vz[j, idx] = vz_tmp[:, sy - 1, sx - 1]

comm.Barrier()
if rank == 0:
    for p in range(1, size):
        for j in range(p, len(sites), size):
            comm.Recv(vx[j, :], source=p, tag = 3 * j)
            comm.Recv(vy[j, :], source=p, tag = 3 * j + 1)
            comm.Recv(vz[j, :], source=p, tag = 3 * j + 2)
else:
    for j in range(rank, len(sites), size):
        print("Sending site %d\n" % j)
        comm.Send(vx[j, :], dest=0, tag=3 * j)
        comm.Send(vy[j, :], dest=0, tag=3 * j + 1)
        comm.Send(vz[j, :], dest=0, tag=3 * j + 2)
comm.Barrier()

for i in range(rank, len(sites), size):
    sx, sy, sz = sites[i]
    np.savetxt(f'{base_x}{sx:04d}_{sy:04d}_{sz:04d}.dat', vx[i, :],
            newline='\n', fmt='%.6e')
    np.savetxt(f'{base_y}{sx:04d}_{sy:04d}_{sz:04d}.dat', vy[i, :],
            newline='\n', fmt='%.6e')
    np.savetxt(f'{base_z}{sx:04d}_{sy:04d}_{sz:04d}.dat', vz[i, :],
            newline='\n', fmt='%.6e')

