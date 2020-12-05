#!/ccs/proj/geo112/hzfmer/summit/opt/anaconda3/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 5 2019

@author: Zhifeng Hu <zhh076@ucsd.edu>
Input:  *.srf, srf source file from CyberShake
        fault_idx.bin indices of subfaults in the mesh

Output: gp_src.bin, gp-type source for IFAULT=5
"""

import numpy as np
import os
import struct

idx = np.genfromtxt('subfaults.idx')
idx = idx[:, : 3].astype('int')
stf_type = 3.0
res = np.fromfile('gp_src.bin.orig', dtype='float32').reshape((-1, 15))
res = res[:, 3:]

f_out = open('gp_src.bin', 'wb')
fmt = "<3i12f"
for i in range(len(res)):
    f_out.write(struct.pack(fmt, *idx[i, :], *res[i, :]))
f_out.close()

try:
    os.symlink('gp_src.bin', 'gp_src.bin_0')
except FileExistsError:
    os.remove('gp_src.bin_0')
    os.symlink('gp_src.bin', 'gp_src.bin_0')
    pass
