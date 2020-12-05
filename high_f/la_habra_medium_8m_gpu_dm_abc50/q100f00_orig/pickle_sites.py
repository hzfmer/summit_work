import numpy as np
import re
import pickle
import collections
import sys

if len(sys.argv) < 4:
    print("Usage: pickle.py mx my stat.idx")
mx = int(sys.argv[1])
my = int(sys.argv[2])
stats = sys.argv[3]

temp = np.genfromtxt(stats, delimiter=" ", dtype="S8, i4, i4")
syn_sites = ((name.decode('UTF-8').replace('_', ''), ix, iy) for name, ix, iy in temp \
             if 0 < int(ix) < mx - 1 and 0 < int(iy) < my - 1)

with open('param.sh', 'r') as fid:
    for line in fid:
        if 'DT' in line:
            dt = float(re.findall(r'DT (\d+\.\d+)', line)[0])
        if 'NTISKP' in line:
            ntiskp = int(re.findall(r'NTISKP (\d+)', line)[0])


with open('stat.txt', 'r') as fid:
    nsite = int(fid.readline())
    vel = collections.defaultdict(dict)
    for line in fid:
        site_name = next(syn_sites)[0]
        line = line.split(' ')
        ix, iy, iz = map(int, line)
        print(f'output_sfc/SX_0_{ix:04d}_{iy:04d}_{iz:04d}.dat')
        vel[site_name]['dt'] = dt * ntiskp
        vel[site_name]['X'] = np.genfromtxt(f'output_sfc/SX_0_{ix:04d}_{iy:04d}_{iz:04d}.dat', dtype='float32')
        vel[site_name]['Y'] = np.genfromtxt(f'output_sfc/SY_0_{ix:04d}_{iy:04d}_{iz:04d}.dat', dtype='float32')
        vel[site_name]['Z'] = np.genfromtxt(f'output_sfc/SZ_0_{ix:04d}_{iy:04d}_{iz:04d}.dat', dtype='float32')

with open('vel_sites.pickle', 'wb') as fid:
    pickle.dump(vel, fid, protocol=pickle.HIGHEST_PROTOCOL)
