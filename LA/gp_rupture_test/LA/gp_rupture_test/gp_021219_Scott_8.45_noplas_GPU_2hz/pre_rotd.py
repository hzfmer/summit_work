#!/ccs/home/hzfmer/file_back/programs/anaconda3/bin/python

import numpy as np
import os
import re
import subprocess
from filter_BU import filt_B

M = re.findall("\d\.\d{2}", os.getcwd())[0]

seis = np.load(f'/gpfs/alpine/geo112/scratch/hzfmer/LA/gp_rupture_test/LA/data_for_Zhifeng/seis{M.split(".")[0]}.npy', allow_pickle=True)[()]
sites = np.genfromtxt("stat.txt", dtype='int', skip_header=1)
site_names = [line.strip('\n') for line in open("site_names.txt")]
idx = np.random.randint(len(sites))
ix, iy, _ = sites[idx]

with open("rotd50_inp.cfg", 'w') as fh:
    fh.write("1\n")  # jInterp
    fh.write("1\n")  # nPair
    fh.write("1\n")  # nHead
    fh.write("accx.txt\n")
    fh.write("accy.txt\n")
    fh.write("rotd50_res.txt\n")

fx = f"output_sfc/SX_0_{ix}_{iy}_0001.dat"
fy = f"output_sfc/SY_0_{ix}_{iy}_0001.dat"
velx = np.loadtxt(fx)
vely = np.loadtxt(fy)
velx, vely = seis[site_names[idx]][1, :] / 100, seis[site_names[idx]][0, :] / 100
nt = 8000
dt = 0.05
velx = filt_B(velx, 1 / dt, 0.01, 2)
vely = filt_B(vely, 1 / dt, 0.01, 2)
accx, accy = (np.diff(x, append=0) / dt for x in (velx, vely))
np.savetxt("accx.txt", accx, fmt="%.8f", newline="\n", header=f"{nt} {dt}", comments='')
np.savetxt("accy.txt", accy, fmt="%.8f", newline="\n", header=f"{nt} {dt}", comments='')

sa_cb = np.fromfile(f'/gpfs/alpine/geo112/scratch/hzfmer/LA/gp_rupture_test/LA/data_for_Zhifeng/sa_cb_site_M{M}_T1.bin', dtype='float32')
print(f"SA-1s at {site_names[idx]}({ix}, {iy}) from CyberShake is: ", sa_cb[idx])

os.system("./rotd50")

with open("rotd50_res.txt") as f_rotd:
    for i, line in enumerate(f_rotd):
        if i == 45:
            sa = float(line.split()[-1]) / 9.8
            break
print(sa)
output = subprocess.run("sed -n '46p' rotd50_res.txt | awk -F' ' '{{print$4}}'", shell=True, stdout=subprocess.PIPE)
print(output.stdout.decode('UTF-8'))
