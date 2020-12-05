#!/ccs/home/hzfmer/file_back/programs/anaconda3/bin/python
import numpy as np

tmax = 30
dt = 0.001
nt = 5000
nsrc = 1560
npend = int(tmax / dt) - 5000
zeros = np.zeros((nsrc, npend), dtype='float32')
for comp in ["xx", "yy", "zz", "xz", "yz", "xy"]:
    src = np.fromfile(f'source_{comp}', dtype='float32').reshape(nsrc, -1)
    if len(src[0, :]) == nt:
        src = np.concatenate((src, zeros), axis=1)
        src.tofile(f'source_{comp}')

