#!/usr/bin/env python

import numpy as np
import sys
from mpi4py import MPI
from numpy import log, exp, sqrt
import re
import os

comm = MPI.COMM_WORLD
ncpus = comm.Get_size()
rank = comm.Get_rank()

if rank == 0:
    print(f'ncpus: {ncpus}')

params = {'6.35': (354, 212, 53), '6.35': (144, 230, 10), '7.35': (980, 240, 10), '8.15': (5116, 220, 10), '8.45': (5116, 220, 10)}
M = re.findall(r"\d+\.\d+", os.getcwd().split('/')[-1])[0]
(nx, nz, inc) = params[M]
nt = 40000
dt = 0.005
m0 = np.zeros((ncpus, 1))

bskip = (3 + 6 * nt) * 4
rskip = bskip - 12
if ncpus != nz // inc:
    print('ERROR! ncpus should be {nz // inc}')
    sys.exit(-1)

with open("momrate.dat", 'rb') as fh:
    for i in range(rank * inc , (rank + 1) * inc):
        for j in range(nx):
            skip = bskip * (i * nx + j) + 12
            fh.seek(skip, 0)
            tmp = np.frombuffer(fh.read(rskip), dtype='f').reshape((-1, 6))
            tmp = np.sum(tmp, axis=0) * dt
            tmp = np.hstack((tmp, tmp[3:]))  # The other 3 comp -> 9 in total
            m0[rank] += sqrt(np.sum(tmp ** 2) / 2)
if rank == 0:
    for p in range(1, ncpus):
        comm.Recv(m0[p], source=p, tag=101)
else:
    comm.Send(m0[rank], dest=0, tag=101)
comm.Barrier()

if rank == 0:
    print(f'Total moment from momrate.dat is {sum(m0)}')
    mag = 2. / 3. * (np.log10(sum(m0)) - 9.1)
    print(f'Magnitude is {mag}')
    with open('momrate_est.txt', 'a') as f:
        f.write(f'Total moment from momrate.dat is {sum(m0)}')
        f.write(f'Magnitude is {mag}')
