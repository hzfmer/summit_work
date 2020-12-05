import numpy as np

nx, ny, nz = 100, 100, 100
nvar = 3
std = 0.05
with open("hom.bin", "wb") as fid:
    data = np.arange(nz * ny * nx * nvar).reshape(nz, ny, nx, nvar) + 10
    data.astype('float32').tofile(fid)

with open("ssh.bin", "wb") as fid:
    data_ssh = np.arange(nz * ny * nx).reshape(nz, ny, nx) / 0.1 + 0.05
    data_ssh.astype('float32').tofile(fid)

with open('het_benchmark.bin', "wb") as fid:
    data_out = np.copy(data)
    for i in range(nz):
        data_out[i, :, :, 0] = data[i, :, :, 0] / (1 - data_ssh[i, :, :])
        data_out[i, :, :, 1] = data[i, :, :, 1] / (1 - data_ssh[i, :, :])
        data_out[i, :, :, 2] = data[i, :, :, 2] * (1 + data_ssh[i, :, :])
    data_out.astype('float32').tofile(fid)

