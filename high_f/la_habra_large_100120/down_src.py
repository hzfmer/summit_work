#!/ccs/home/hzfmer/file_back/programs/anaconda3/bin/python

import numpy as np
import collections
import shutil
import sys

# The sources are located in the second block
mx, my, mz = 6480, 4968, 380
dh = 24
nsrc = 15625
nt_src = 12500

# [iz, ix, iy]
index = np.genfromtxt('subfaults.idx', dtype='int')[:,:3]
idx = np.zeros((nsrc, 3), dtype='int32')


# For the second block in Z direction, there are (mz - 7) fine grids above, and (mz - 8)dh distance
# For the second block in Y direction, shift one more fine grid further, so yoff is dh smaller
seen = collections.defaultdict(list)
for i in range(nsrc):
    ix, iy, iz = index[i]
    seen[(ix, iy, iz)] += i,
    idx[i] = np.array([ix, iy, iz], dtype='int32')
print(f'Non-duplicate subfaults {len(seen)} / {nsrc}')

# Superimpose duplicate subfaults
shutil.copyfile('momrate.dat', 'source_1')

with open('momrate.dat', 'rb') as fid_in, \
        open("source_1", "wb") as fid:
    for k, v in seen.items():
        print(f"Duplicate subfault {k}, {len(v)}")
        tmp = np.zeros((6 * nt_src,), dtype='float32')
        for i in v:
            fid_in.seek(0, 0)
            fid.seek(i * (6 * 4 * nt_src + 12), 0)
            idx[i].astype('int32').tofile(fid)
            tmp += np.fromfile(fid_in, dtype='float32', count=6*nt_src,
                    offset=i * (6 * 4 * nt_src + 12) + 12)
        for i in v:
            fid.seek(i * (6 * 4 * nt_src + 12) + 12, 0)
            tmp.tofile(fid)

# Touch source_1 in case not existed
try:
    open('source_0', 'ab', 0).close()
    open('source_2', 'ab', 0).close()
except OSError:
    pass
