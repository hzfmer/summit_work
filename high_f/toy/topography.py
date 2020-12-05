import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

nx, ny = 3500, 3500
mx, my = nx // 54 * 54, ny // 54 * 54
mx, my = 1000, 1000
px, py = (nx - mx) // 2, (ny - my) // 2
pad = 8
grid = np.fromfile('surf.grid', dtype='float64').reshape(ny, nx, 3)
lon = grid[:, :, 0]
lat = grid[:, :, 1]
del grid

print(f"Left lower corner: {lon[0, 0]}, {lat[0, 0]}\n"
      f"Right lower corner: {lon[0, -1]}, {lat[0, -1]}\n"
      f"Left upper corner: {lon[-1, 0]}, {lat[-1, 0]}\n"
      f"Right upper corner: {lon[-1, -1]}, {lat[-1, -1]}")

dem_head = np.genfromtxt('dem.txt', dtype='float')
print(dem_head)

NX, NY = int(dem_head[0]), int(dem_head[1])
left, right, bot, top = dem_head[2:]
dem = np.fromfile('dem.bin', dtype='float32').reshape(NY, NX)

dx = (right - left) / NX
dy = (top - bot) / NY

try:
    raise Exception
    with open('topography.bin', 'rb') as fout:
        buf = fout.read(12)
        mx, my, pad = np.frombuffer(buf, dtype='int32')
        topo = np.frombuffer(fout.read((mx + 2 * pad) * (my + 2 * pad) * 4), dtype='float32').reshape(mx + 2 * pad, my + 2 * pad)

except:
    topo = np.zeros(((mx + 2 * pad), (my + 2 * pad)), dtype='float32')
    idx = np.round((lon - left) / dx)
    idy = np.round((top - lat) / dy)
    for i in range(pad, mx + pad):
        for j in range(pad, my + pad):
            topo[i, j] = dem[int(idy[j + py - pad, i + px - pad]), int(idx[j + py - pad, i + px - pad])]

    f = interpolate.interp2d(pad + np.arange(mx), pad + np.arange(my), topo[pad : mx + pad, pad : my + pad], kind='cubic')
for i in range(mx + 2 * pad):
    for j in range(pad):
        topo[i, j] = f(i, j)
    for j in range(my + pad, my + 2 * pad):
        topo[i, j] = f(i, j)

    for i in range(pad):
        for j in range(pad, my + pad):
            topo[i, j] = f(i, j)
    for i in range(mx + pad, mx + 2 * pad):
        for j in range(pad, my + pad):
            topo[i, j] = f(i, j)


    with open('topography.bin', 'w') as fout:
        header = np.array([mx, my, pad], dtype='int32')
        header.tofile(fout)
        topo.tofile(fout)

fig, ax = plt.subplots(1, 2)
fig.subplots_adjust(wspace=0.5, bottom = 0.2)
vmin = np.min(dem)
vmax = np.max(dem)

dlon, dlat = np.meshgrid(np.linspace(left, right, NX), np.linspace(top, bot, NY))
pc0 = ax[0].pcolormesh(dlon[::5, ::5], dlat[::5, ::5], dem[::5, ::5], vmin=vmin, vmax=vmax, cmap='terrain')
pc1 = ax[1].pcolormesh(lon[py : py + my : 5, px : px + mx : 5], lat[py : py + my : 5, px : px + mx : 5], topo[pad : mx + pad : 5, pad : my + pad : 5].T, vmin=vmin, vmax=vmax, cmap='terrain')
cax = fig.add_axes([0.2, 0.05, 0.6, 0.05])
cb = fig.colorbar(pc1, cax=cax, orientation='horizontal')
cb.ax.set_xlabel('Elevation (m)')
fig.savefig(f"check_DEM.png", dpi=600, bbox_inches='tight', pad_inches=0.05)
