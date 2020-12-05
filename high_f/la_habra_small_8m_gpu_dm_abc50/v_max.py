#!/usr/bin/env python
'''
Maximum of vp in mesh, assuming the last several layers
'''

import numpy as np
from sys import argv, exit

if len(argv) != 2:
    print(f"Usage:{argv[0]} + velocity mesh")
    exit

if "large" in argv[1]:
    nx = 9000
    ny = 6750
elif "small" in argv[1]:
    nx = 1400
    ny = 1400

nz = 500
nvar = 3
nflt = 4
k = 50 
vp_max = 0.0
vs_max = 0.0
vs_min = 10000.0

f = open(argv[1], 'rb')

for i in range(k):
    buf = f.read(nvar * nflt * nx * ny)
    vs_tmp = np.frombuffer(buf, dtype='f', count=nvar*nx*ny).reshape(-1,3)[:,1]
    vs_min = np.min([vs_min, vs_tmp.min()])
    print("Vs_min in the {0}th layer is {1}\n".format(k, vs_tmp.min()))

f.seek(nvar * nflt * nx * ny * (nz - k), 0)

for i in range(k):
    buf = f.read(nvar * nflt * nx * ny)
    vp_tmp = np.frombuffer(buf, dtype='f', count=nvar*nx*ny).reshape(-1,3)[:,0]
    vs_tmp = np.frombuffer(buf, dtype='f', count=nvar*nx*ny).reshape(-1,3)[:,1]
    vs_max = np.max([vs_max, vs_tmp.max()])
    vp_max = np.max([vp_max, vp_tmp.max()])
    print("Vp_max in the {0}th layer is {1}\n".format(nz-k+i, vp_tmp.max()))
    print("Vs_max in the {0}th layer is {1}\n".format(nz-k+i, vs_tmp.max()))
print("Vp_max in the last {0} layers of {1} is {2}".format(k, argv[1], vp_max))
print("Vs_max in the last {0} layers of {1} is {2}".format(k, argv[1], vs_max))
print("Vs_min in the first {0} layers of {1} is {2}".format(k, argv[1], vs_min))

