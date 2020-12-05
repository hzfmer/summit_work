#!/ccs/home/hzfmer/file_back/programs/anaconda3/bin/python
import numpy as np

z = 108 * 8
fname = 'source.txt'
data = np.genfromtxt(fname, dtype='int', skip_header=10)
fid = open(fname, 'r')
header = fid.readlines()[:10]
fid.close()

shallow = np.min(data[:, 3])
deep = np.max(data[:, 3])
mid = (shallow + deep) // 2

data[:, 3] = (-0.5 - (mid - data[:, 3]) / (deep - shallow) * 0.5) * z

fout = open('test.txt', 'w')
fout.writelines(header)
np.savetxt(fout, data, fmt='%d')
fout.close()
