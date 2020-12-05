#!/ccs/home/hzfmer/file_back/programs/anaconda3/bin/python
"""Format synthetics at sites for bbp processing by SCEC
Usage:
    python rename_sites.py {model}
    E.g. python rename_sites.py noqf_orig_vs500
"""
import pandas as pd
import shutil
import sys
from scipy.signal import resample
from post_processing.la_habra import *

if len(sys.argv) == 1:
    model = 'q100f00_orig_vs500'
    rec_list = '15sites.txt'
else:
    model = sys.argv[1]
    rec_list = sys.argv[2]
tmax, dt, tskip, _, _ = read_param(model)
dt *= tskip
print(f'tmax={tmax}, dt={dt}')

dt_new = 0.02
N = len(np.arange(0, tmax, dt_new))
sites = pd.read_csv('la_habra_small_statlist_3456.idx', sep=' ', header=None, index_col=False)
rec_sites = pd.read_csv(rec_list, delimiter=" ", header=None, index_col=False)[0].tolist()
print(sites.head(5))
count = 0
for i in range(len(sites)):
    if sites.iloc[i][0] in rec_sites:
        count += 1
        syn = []
        for comp in 'XYZ':
            data = np.genfromtxt(f'{model}/output_sfc/S{comp}_0_{sites.iloc[i][1]:04d}_{sites.iloc[i][2]:04d}_0001.dat', dtype='float32')
            data = resample(data, N)
            syn.append(data)

            # shutil.copyfile(f'noqf_orig/output_sfc/S{comp}_0_{sites.iloc[i][1]:04d}_{sites.iloc[i][2]:04d}_0001.dat',
            #                 f'noqf_orig/output_sfc/{sites.iloc[i][0]}.dat')
        syn = np.asarray(syn)
        syn = np.vstack((np.arange(0, tmax, dt_new, dtype='float64'), syn))

        np.savetxt(f'{model}/output_sfc/{sites.iloc[i][0]}.dat', syn.T, fmt='%.8f', delimiter=' ')
print("Shape of output: ", syn.shape)
print("Number of sites in output_sfc renamed: ", count)
