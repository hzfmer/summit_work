#!/ccs/home/hzfmer/file_back/programs/anaconda3/bin/python

import numpy as np
import collections
import shutil
import sys

mx, my, mz = 9504, 7020, 384
# The sources are located in the second block
dh = 20
nsrc = 15625
nt_src = 5000
gbuf = 1000

# [iz, ix, iy]
index = np.genfromtxt('fault_idx.txt', dtype='int')

idx = np.zeros((nsrc, 3), dtype='int32')

# For the second block in Z direction, there are (mz - 7) fine grids above, and (mz - 8)dh distance
# For the second block in Y direction, shift one more fine grid further, so yoff is dh smaller
seen = collections.defaultdict(list)
xoff, yoff, zoff = 0, -dh, (mz - 8) * dh
for i in range(nsrc):
    iz, ix, iy = index[i]
    seen[(ix, iy, iz)] += i,
    idx[i] = np.array([ix, iy, iz], dtype='int32')
print(f'Non-duplicate subfaults {len(seen)} / {nsrc}')

# Superimpose duplicate subfaults
shutil.copyfile('momrate.dat', 'source_0')

with open('momrate.dat', 'rb') as fid_in, \
        open("source_0", "r+b") as fid:
    for k, v in seen.items():
        tmp = np.zeros((6 * nt_src,), dtype='float32')
        for i in v:
            fid_in.seek(i * (6 * 4 * nt_src + 12) + 12, 0)
            fid.seek(i * (6 * 4 * nt_src + 12), 0)
            idx[i].astype('int32').tofile(fid)
            tmp += np.frombuffer(fid_in.read(6 * 4 * nt_src), dtype='float32')
        for i in v:
            fid.seek(i * (6 * 4 * nt_src + 12) + 12, 0)
            tmp.tofile(fid)

# Touch source_1 in case not existed
try:
    open('source_1', 'ab', 0).close()
except OSError:
    pass
