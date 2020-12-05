import numpy as np
import matplotlib.pyplot as plt

def estimateRhofromVs(vs):
    N = len(vs)
    rho = np.zeros(N)
    for i in range(N):
        if vs[i] <= 333:
            rho[i] = 2000 + 190 * (vs[i] - 166.5) / 333
        elif vs[i] < 560:
            rho[i] = 2350
        elif vs[i] < 640:
            rho[i] = 2450
        else:
            rho[i] = 2750
    return rho


def estimateVpfromVs(vs):
    N = len(vs)
    vp = np.zeros(N)
    for i in range(N):
        if vs[i] <= 333:
            vp[i] = 2.5 * vs[i]
        elif vs[i] < 560:
            vp[i] = 2 * vs[i]
        elif vs[i] < 640:
            vp[i] = 1.9 * vs[i] 
        else:
            vp[i] = 1.73 * vs[i]
    return vp


nx, ny, nz = 540, 160, 280

vs_0 = [1000] * (nz // 4) + [1500] * (nz // 4) + [2000] * (nz // 2)
vp_0 = estimateVpfromVs(vs_0)
rho_0 = estimateRhofromVs(vs_0)


vs = np.tile(vs_0, (nx, ny, 1)).T
vp = np.tile(vp_0, (nx, ny, 1)).T
rho = np.tile(rho_0, (nx, ny, 1)).T
print(vp.shape)

data = np.column_stack((x.flatten() for x in (vp, vs, rho)))
data.astype('float32').tofile('mesh')

fig, ax = plt.subplots()
im = ax.imshow(vs[:, ny // 2, :])
fig.colorbar(im)
fig.savefig(f"mesh.png", dpi=600, bbox_inches='tight', pad_inches=0.05)
