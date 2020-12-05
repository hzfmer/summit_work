#!/u/sciteam/hu2/anaconda2/envs/py35/bin/python

'''
Maximum of vp in mesh, assuming the last several layers
'''

import numpy as np
from sys import argv, exit

if len(argv) != 2:
    print("Usage:{0} + velocity mesh".format(argv[0]))
    exit

nx = 4200
ny = 6320
nz = 400
nvar = 3
nflt = 4
k = 0
vp_max = 0.0
vs_max = 0.0
vs_min = 10000.0

f = open(argv[1], 'rb')

#for i in range(k):
#    buf = f.read(nvar * nflt * nx * ny)
#    vs_tmp = np.frombuffer(buf, dtype='f', count=nvar*nx*ny).reshape(-1,nvar)[:,1]
#    vs_min = np.min([vs_min, vs_tmp.min()])
#    print("Vs_min in the {0}th layer is {1}\n".format(k, vs_tmp.min()))

f.seek(nvar * nflt * nx * ny * (nz - k), 0)

for i in range(nz-k):
    buf = f.read(nvar * nflt * nx * ny)
    vp_tmp = np.frombuffer(buf, dtype='f', count=nvar*nx*ny).reshape(-1,nvar)[:,0]
    vs_tmp = np.frombuffer(buf, dtype='f', count=nvar*nx*ny).reshape(-1,nvar)[:,1]
    vs_max = np.max([vs_max, vs_tmp.max()])
    vp_max = np.max([vp_max, vp_tmp.max()])
    vs_min = np.min([vs_min, vs_tmp.min()])
    print("Vp_max in the {0}th layer is {1}\n".format(i, vp_tmp.max()))
    print("Vs_max in the {0}th layer is {1}\n".format(i, vs_tmp.max()))
    print("Vs_min in the {0}th layer is {1}\n".format(i, vs_tmp.min()))
print("Vp_max in the last {0} layers of {1} is {2}".format(nz, argv[1], vp_max))
print("Vs_max in the last {0} layers of {1} is {2}".format(nz, argv[1], vs_max))
print("Vs_min in the first {0} layers of {1} is {2}".format(nz, argv[1], vs_min))

