#!/ccs/proj/geo112/hzfmer/summit/opt/anaconda3/bin/python

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import resample

nsite_x = 8
nsite_y = 4
ix = 4
iy = 3

nx = 400
ny = 400
nz = 400

sx = np.fromfile('output_sfc/SX_0_0750000', dtype='f').reshape((-1, 2, nsite_y, nsite_x))
sy = np.fromfile('output_sfc/SY_0_0750000', dtype='f').reshape((-1, 2, nsite_y, nsite_x))
sz = np.fromfile('output_sfc/SZ_0_0750000', dtype='f').reshape((-1, 2, nsite_y, nsite_x))

vx_surf = sx[:, 0, iy, ix]
vy_surf = sy[:, 0, iy, ix]
vz_surf = sz[:, 0, iy, ix]
v_surf = np.column_stack((vx_surf, vy_surf, vz_surf))
v_surf = resample(v_surf, 750000)
np.savetxt('stf_surf.txt', v_surf, delimiter=' ', fmt='%.10f')

vx_bot = sx[:, 1, iy, ix]
vy_bot = sy[:, 1, iy, ix]
vz_bot = sz[:, 1, iy, ix]
v_bot = np.column_stack((vx_bot, vy_bot, vz_bot))
v_bot = resample(v_bot, 750000)
np.savetxt('stf_bot.txt', v_bot, delimiter=' ', fmt='%.10f')


