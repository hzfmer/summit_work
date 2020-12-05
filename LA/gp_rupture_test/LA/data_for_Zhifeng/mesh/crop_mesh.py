#!/ccs/proj/geo112/hzfmer/summit/opt/anaconda3/bin/python

import numpy as np

nx, ny, nz, nvar = 6320, 4200, 400, 3
data = np.fromfile('awp.USC.media', dtype='float32').reshape((nz, ny, nx, nvar))

data = data[:, :, :6048, :]
data.tofile('awp.USC.crop.media')
