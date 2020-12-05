import numpy as np

nbit = 4  # float bit

nx, ny, nz, nvar = 19440, 14904, 700, 3
mx, my = nx * 2 // 3 + 1200, ny * 2 // 3
endx, endy = mx + 2592, my + 2160
nlayer = nx * ny * nvar
with open("mesh_0", "wb") as fid:
    for i in range(nz):
        data = np.fromfile("../../cvm/la_habra_large_cvmsi_8m.media", dtype='float32',
            count=nlayer, offset=nlayer * i * nbit).reshape(ny, nx, nvar)
        data[my : endy, mx : endx, :].tofile(fid)
