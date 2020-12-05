#! /usr/bin/env python

import sys 
sys.path.append("/ccs/home/hzfmer")
from pynga import *
from pynga.utils import *
import os
import re

data_dir = "/lustre/atlas/proj-shared/geo112/huzf/LA/gp_rupture_test/Scott/data/"
os.chdir(data_dir)
res_dir = "/lustre/atlas/proj-shared/geo112/huzf/LA/gp_rupture_test/Scott/results/"
files = list(filter(os.path.isfile, os.listdir(data_dir)))
print(files)
models = ['ASK', 'BSSA', 'CB', 'CY']
models = ['ASK', 'CB', 'CY']
models = ['BSSA']
Ftype='SS'

for model_name in models:
    print(model_name)
    for file in files:
        period = float(re.findall(r'\d+', file)[0])
        fout = open(res_dir + model_name + '_' + file, 'w')
        with open(file) as fid:
            for line in fid:
                line2 = line.split()
                Mw = float(line2[1])
                Rjb = float(line2[2])
                Vs30 = float(line2[3])
                median, sigmaT, tau, sigma = NGA14(model_name, Mw, Rjb, Vs30, period, rake=179.9, epislon=0, Ftype=Ftype)
                fout.write("%s %f %f %f %f\n" % (line.strip(), median[0], sigmaT[0], tau[0], sigma[0]))
        fout.close()

