import numpy as np

copnt = ['xx', 'yy', 'zz', 'xz', 'yz', 'xy']
nt = 5000
M0 = 1e15
for i in range(len(copnt)):
    data = np.random.rand(nt) * M0
    data.astype('float32').tofile(f'source_0_{copnt[i]}')
