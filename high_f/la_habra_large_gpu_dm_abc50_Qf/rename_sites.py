#!/ccs/home/hzfmer/file_back/programs/anaconda3/bin/python
"""Format synthetics at sites for bbp processing by SCEC
Usage:
    python rename_sites.py {model}
    E.g. python rename_sites.py noqf_orig_vs500
"""
import pandas as pd
import shutil
from pathlib import Path
import sys
from scipy.signal import resample
from post_processing.la_habra import *

if len(sys.argv) < 4:
    print(f'Usage: python {sys.argv[0]} fid_syn_sites fid_rec_sites model dt_new')
    print(f'E.g.\n'
          f'    stat_name_idx.txt: site_name ix iy'
          f'    la_habra_large_statlist_070120.txt: site_name lon lat'
          f'    0.02'
          f'    q100f00_orig')
    sys.exit(-1)
else:
    fid_syn_sites = sys.argv[1]
    fid_rec_sites = sys.argv[2]
    dt_new = float(sys.argv[3])
    model = sys.argv[4] if len(sys.argv) == 5 else ""

# get dt and tmax from param.sh
tmax, dt, tskip, _, _ = read_param(model)
dt *= tskip  # decimation
t = np.arange(0, tmax, dt_new)
N = len(t)
print(f'tmax={tmax}, dt={dt}, resampling dt={dt_new}')

# site_name ix iy
sites = pd.read_csv(fid_syn_sites, sep=' ', header=None, index_col=False)
print(f'synthetics include:\n', sites.head(5))

# site_name column only
rec_sites = pd.read_csv(fid_rec_sites, delimiter=" ", comment='#', header=None, index_col=False)[0].tolist()

count = 0
for i in range(len(sites)):
    if sites.iloc[i][0] in rec_sites:
        count += 1
        syn = t.copy()
        for comp in 'XYZ':
            data = np.genfromtxt(Path(model, f'output_sfc/S{comp}_0_{sites.iloc[i][1]:04d}_{sites.iloc[i][2]:04d}_0001.dat'), dtype='float32')
            data = resample(data, N)
            syn = np.vstack((syn, data))

        np.savetxt(Path(model, f'output_sfc/{sites.iloc[i][0]}.dat'), syn.T, fmt='%.7f', delimiter=' ')
print("Shape of output: ", syn.shape)
print("Number of sites in output_sfc renamed: ", count)
