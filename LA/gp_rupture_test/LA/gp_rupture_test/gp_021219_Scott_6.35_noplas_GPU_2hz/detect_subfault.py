#!/ccs/proj/geo112/hzfmer/summit/opt/anaconda3/bin/python

import numpy as np
import sys
from glob import glob
import os

from posix import symlink
from os.path import basename


flist=glob(f'{os.getcwd()}/fault*')
print(f"Number of nodes with some sources: {len(flist)}")
select = np.random.randint(len(flist), size=20)
for x in select:
    f = os.path.basename(flist[x])
    data = np.fromfile(f, dtype='float32')
    if max(np.abs(data)) < 1e-12:
        print(f"{f}:  Likely all zeros")
    else:
        idx = np.argmax(np.abs(data))
        print(f'{f}  :  {idx}')
        sys.exit(-1)

