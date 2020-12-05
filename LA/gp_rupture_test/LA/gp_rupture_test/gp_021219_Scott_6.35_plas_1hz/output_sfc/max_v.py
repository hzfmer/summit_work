#!/u/sciteam/hu2/anaconda3/bin/python

import numpy as np


nt = 1000
nx = 6320
ny = 4200

#sx = np.fromfile('SX', dtype='float32').reshape(
f = open('SX', 'rb')

v_max = 0

for i in range(nt):
    for j in range(nx):
        buf = f.read(ny * 4)
        dat = np.frombuffer(buf, dtype='float32')
        v_max = np.max((v_max, dat.max()))
    print(f'time_step={i}, v_max={v_max}', flush=True)
print(v_max)

