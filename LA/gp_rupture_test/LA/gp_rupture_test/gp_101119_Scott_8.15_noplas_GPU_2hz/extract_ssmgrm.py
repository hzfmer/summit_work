#!/ccs/home/hzfmer/file_back/programs/anaconda3/bin/python
import numpy as np
import os
from filter_BU import filt_B
from read_params import read_params


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

vx = np.zeros((len(sites), nt // ntskp))
vy = np.zeros((len(sites), nt // ntskp))
vz = np.zeros((len(sites), nt // ntskp))

for i in range(0, nt // skip):
    idx = np.arange(i * wstep, (i + 1) * wstep)
    print(f"Reading timestep at {(i + 1) * skip}")
    print("filename", base_x + f'{(i + 1) * skip:07d}')
    vx_tmp = np.fromfile(base_x + f'{(i + 1) * skip:07d}', dtype='f').reshape(-1, ny, nx)
    vy_tmp = np.fromfile(base_y + f'{(i + 1) * skip:07d}', dtype='f').reshape(-1, ny, nx)
    vz_tmp = np.fromfile(base_z + f'{(i + 1) * skip:07d}', dtype='f').reshape(-1, ny, nx)
    for j in range(len(sites)):
        print(f"reading site {j} / {len(sites)} at time {i * skip}")
        sx, sy, _ = sites[j, :]
        vx[j, idx] = vx_tmp[:, sy - 1, sx - 1]
        vy[j, idx] = vy_tmp[:, sy - 1, sx - 1]
        vz[j, idx] = vz_tmp[:, sy - 1, sx - 1]

for i in range(len(sites)):
    sx, sy, sz = sites[i]
    np.savetxt(f'{base_x}{sx:04d}_{sy:04d}_{sz:04d}.dat', vx[i, :],
            newline='\n', fmt='%.6e')
    np.savetxt(f'{base_y}{sx:04d}_{sy:04d}_{sz:04d}.dat', vy[i, :],
            newline='\n', fmt='%.6e')
    np.savetxt(f'{base_z}{sx:04d}_{sy:04d}_{sz:04d}.dat', vz[i, :],
            newline='\n', fmt='%.6e')

