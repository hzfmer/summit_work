#!/usr/bin/env python

import numpy as np

f_in = open('subfaults.idx', 'r')
f_out = open('temp.idx', 'w')
s = f_in.readlines()

nx = 980
nz = 240
for j in range(nz):
    f_out.writelines(reversed(s[nx * j : nx * (j + 1)]))

f_in.close()
f_out.close()
