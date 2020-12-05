import numpy as np

mx, my = 3456, 3456
nz, mz = 1500, 1400

with open('la_habra_small_cvmsi_8m_3456_nz1400.media', 'wb') as fid:
    for i in range(mz):
        for j in range(my):
            data = np.fromfile('la_habra_small_cvmsi_8m_3456.media', dtype='float32', count=3 * mx)
            data.tofile(fid)

