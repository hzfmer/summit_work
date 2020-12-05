import numpy as np

data = np.genfromtxt('subfaults.idx.back')
xoff = 3550
yoff = 2050
for i in range(len(data)):
    data[i][0] = data[i][0] - xoff
    data[i][1] = data[i][1] - yoff

fmt = '%d ' * 3 + '%e ' * 2 + '%f ' * 3
np.savetxt('subfaults.idx', data, fmt=fmt) 
