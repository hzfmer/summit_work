"""
Input
-----
statlist.idx : str
    Name of file containing site_name, ix, iy
mx : int
    Number of grids along X
my : int
    Number of grids along Y
nd : int [optional]
    Number of layers cropped

Output
------
stat.txt : ASCII file
    Number of sites; (ix, iy, iz) of each site
vel_sites.pickle : pickle file
    Dictionary of three component velocities at each site, site_name as key

"""
import numpy as np
import re
import pickle
import collections
import sys

if len(sys.argv) < 4:
    print("Usage: pickle.py statlist.idx mx my [nd]")
    print("--------")
    print("la_habra_medium_statlist_5940.idx:")
    print("Line 0: number of site")
    print("[ix, iy, iz]")

fh_stat_idx = sys.argv[1]
mx = int(sys.argv[2])
my = int(sys.argv[3])
nd = 250 if len(sys.argv) < 5 else int(sys.argv[4])

stations = np.genfromtxt(fh_stat_idx, delimiter=" ", dtype="S8, i4, i4")

with open('param.sh', 'r') as fid:
    for line in fid:
        if 'DT' in line:
            dt = float(re.findall(r'DT (\d+\.\d+)', line)[0])
        if 'NTISKP' in line:
            ntiskp = int(re.findall(r'NTISKP (\d+)', line)[0])

stat_idx = []
vel = collections.defaultdict(dict)
for name, ix, iy in stations:
    if nd < int(ix) < mx - nd and nd < int(iy) < my - nd:
        iz = 1  # By default, site locates on the surface
        stat_idx.append([ix, iy, iz])
        site_name = name.decode('UTF-8').replace('_', '')
        vel[site_name]['dt'] = dt * ntiskp  # Decimation with ntiskp
        vel[site_name]['X'] = np.genfromtxt(f'output_sfc/SX_0_{ix:04d}_{iy:04d}_{iz:04d}.dat', dtype='float32')
        vel[site_name]['Y'] = np.genfromtxt(f'output_sfc/SY_0_{ix:04d}_{iy:04d}_{iz:04d}.dat', dtype='float32')
        vel[site_name]['Z'] = np.genfromtxt(f'output_sfc/SZ_0_{ix:04d}_{iy:04d}_{iz:04d}.dat', dtype='float32')
        vel[site_name]['t'] = np.arange(len(vel[site_name]['Z'])) * dt * ntiskp  # Decimation with ntiskp


with open('stat.txt', 'w') as fid:
    fid.write(f'{len(stat_idx)}\n')
    for row in stat_idx:
        fid.write(" ".join(map(str, row)) + "\n")

with open('vel_sites.pickle', 'wb') as fid:
    pickle.dump(vel, fid, protocol=pickle.HIGHEST_PROTOCOL)


