import numpy as np

from numba import jit

@jit(nopython=True)
def vs_to_vp(vs):
    vs = vs / 1000.
    vp = (0.9409 + (2.0947 * vs)
          - (0.8206 * vs ** 2.) 
          + (0.2683 * vs ** 3.) 
          - (0.0251 * vs ** 4.)) * 1000.
    return vp
    
    
@jit(nopython=True)
def vp_to_rho(vp):
    vp = vp / 1000.
    rho = (1.6612 * vp - (0.4721 * vp ** 2)
          + (0.0671 * vp ** 3.) 
          - (0.0043 * vp ** 4.) 
          + (0.000106 * vp ** 5.)) * 1000.
    return rho


@jit(nopython=True)
def to_vs30(v, vs30, k=0.):
    for i in range(my):
        for j in range(mx):
            if vs30[i, j] > 0:
                v[i, j, 1] = (1 - k) * vs30[i, j] + k * v_bot[i, j, 1]
                v[i, j, 0] = vs_to_vp(v[i, j, 1])
                v[i, j, 2] = vp_to_rho(v[i, j, 0])
    return v

import math
@jit(nopython=True)
def get_erf_value(x, xt, y0, y1, x0=0., norm=5):
    """Get scaled erf value
    Input
    -----
    :x (float): location to query
    :x0 (float): start location, here is 0/surface
    :xt (float): width of erf function, here is the tapering depth
    :y0 (float): lower bound of the scaled erf function
    :y1 (float): upper bound of the scaled erf function
    """
    return y0 + (y1 - y0) * (math.erf((x - x0 - xt / 2) * norm / xt) + 1) / 2

@jit(nopython=True)
def find_erf_avg(xt, y0, y1, xq0, xq1, x0=0., N=100, method='harmonic'):
    """Get harmonic average of scaled erf values between (xq0, xq1)
    Input
    -----
    :x0 (float): start location, here is 0/surface
    :xt (float): width of erf function, here is the tapering depth
    :y0 (float): lower bound of the scaled erf function
    :y1 (float): upper bound of the scaled erf function
    :xq0 (float): left bound of query range
    :xq1 (float): right bound of query range
    :N (int): number of points to sample
    :method (str): "harmonic" or "arithmetic"
    """
    xq = np.linspace(xq0, xq1, N)
    res = 0.
    if method.startswith("h"):
        for x in xq:
            res += 1 / get_erf_value(x, xt, y0, y1, x0=x0)
        return N / res
    elif method.startswith("a"):
        for x in xq:
            res += get_erf_value(x, xt, y0, y1, x0=x0)
        return res / N


@jit(nopython=True)
def find_start_erf_value(vs30, v_bot, h, atol=0.2, alpha=0.999, mu=0.5):
    """Search the lower bound of scaled erf function
    :vs30 (float): vs30, constraint of near surface 30 m
    :v_bot (float): vs at the bottom of the tapering layer
    :h (float): thickness of the tapering layer
    :atol (float): absolute tolerance
    :alpha (float): learning_rate
    :mu (float): Nesterov moment optimization parameter
    """
    v_init = vs30
    v_dot = 0.  # average in the shallow 30 m
    v = 0.
    while abs(vs30 - v_dot) > atol:
        v_dot = find_erf_avg(h, v_init, v_bot, 0, 30)
        v_prev = v
        v = mu * v - alpha * (v_dot - vs30)
        v_init += -mu * v_prev + (1 + mu) * v
    return v_init


@jit(nopython=True)
def find_vs_surface_erf(vs30, v_bot, h):
    v = vs30.copy()
    for i in range(my):
        for j in range(mx):
            v[i, j] = find_start_erf_value(vs30[i, j], v_bot[i, j], h)
    return v
            
            
@jit(nopython=True)
def to_vs30_erf(v, v_surf, v_bot, vs30, h, z0, z1, v_orig=None):
    for i in range(my):
        for j in range(mx):
            if vs30[i, j] > 0 and z0[i, j] < h:
                if v_orig is None or v_orig[i, j] > 1000:
                    v[i, j, 0] = find_erf_avg(h, v_surf[i, j, 0], v_bot[i, j, 0], z0[i, j], z1[i, j]) 
                    v[i, j, 1] = find_erf_avg(h, v_surf[i, j, 1], v_bot[i, j, 1], z0[i, j], z1[i, j]) 
                    v[i, j, 2] = find_erf_avg(h, v_surf[i, j, 2], v_bot[i, j, 2], z0[i, j], z1[i, j], method='arithmetic') 
    return v


def read_topo(fname, mx, my):
    """Return topog (my, mx)"""
    pad = 8
    with open(fname, 'rb') as fout:
        mx, my, pad = np.frombuffer(fout.read(12), dtype='int32')
        topo = np.frombuffer(fout.read((mx + 2 * pad) * (my + 2 * pad) * 4), dtype='float32').reshape(mx + 2 * pad, my + 2 * pad).T
        return topo[pad : -pad, pad : -pad]


from osgeo import gdal
raster = gdal.Open(f'../cvm/California_vs30_Wills15_hybrid_7p5c.tif')
geotransform = raster.GetGeoTransform()
originX = geotransform[0]  # Left upper corner
pixelWidth = geotransform[1]
pixelRotateX = geotransform[2]  # Rotation
originY = geotransform[3]
pixelRotateY = geotransform[4]  # Rotation
pixelHeight = geotransform[5]
band = raster.GetRasterBand(1)

arr = np.array(band.ReadAsArray())
arr = arr[:, ::-1]  # reverse latitude
nx, ny = arr.shape
lat = np.arange(ny) * pixelHeight + originY
lat = lat[::-1]  # reverse latitude
lon = np.arange(nx) * pixelWidth + originX

left, right, bot, top = -119.34, -116.77, 33.06, 35.14

idx_lon = np.logical_and(lon >= left, lon <= right)
idx_lat = np.logical_and(lat >= bot, lat <= top)
lon = lon[idx_lon]
lat = lat[idx_lat]
arr = arr[idx_lon, :][:, idx_lat].T
grid_lon, grid_lat = np.meshgrid(lon, lat)
# plt.pcolormesh(lon, lat, arr, cmap='coolwarm', rasterized=True)
# plt.colorbar()
del lon, lat, idx_lon, idx_lat

from scipy import spatial
kdtree = spatial.cKDTree(np.array((grid_lon.flatten(), grid_lat.flatten())).T)

mx, my, mz = 9504, 7020, 3072
try: 
    vs30 = np.fromfile('Vs30_Wills.bin', dtype='float32').reshape(my, mx)
except:
    grids = np.fromfile('surf.grid', dtype='float64').reshape(my, mx, 3)[:, :, :2]

    _, out = kdtree.query(
        np.array(
            (grids[:, :, 0].flatten(),
            grids[:, :, 1].flatten())).T,
        n_jobs=-1)

    vs30 = arr.flatten()[out].reshape(my, mx)
    del arr, grids
    vs30.astype('float32').tofile('Vs30_Wills.bin')

z0 = 30  # vs30
taper = 100  # taper depth
dh = 20
nd = taper // dh  # first layer is repeated
mz = 3072  # number of layers in original mesh
nz = 2900  # number of layers in topo mesh

mesh_base = 'mesh_extlarge_20m_sh005l100.bin'
v_orig = np.fromfile(mesh_base, dtype='float32', count=mx * my * 3).reshape(my, mx, 3)[:, :, 1]
v_bot = np.fromfile(mesh_base, dtype='float32', count=mx * my * 3,
            offset=mx * my * 3 * 4 * nd).reshape(my, mx, 3)
mesh_surf = 'v_surf_sh005l100.bin'
try:
    v_surf = np.fromfile(mesh_surf, dtype='float32').reshape(my, mx, 3)
except Exception as e:
    print(f"v_surf.bin not read: {e}")
    vs_surf = find_vs_surface_erf(vs30, v_bot[:, :, 1], taper)
    vp_surf = vs_to_vp(vs_surf)
    v_surf = np.stack((vp_surf, vs_surf, vp_to_rho(vp_surf)), axis=2)
    del vp_surf, vs_surf
    v_surf.astype('float32').tofile(mesh_surf)

rpt = 1  # repeat surface layer in topo code
topography = read_topo('topography.bin', mx, my) 
H = topography + (nz - 1 - rpt) * dh  # thickness of mesh_topo, one layer is repeated

in_mesh = f'mesh_topo.bin'
in_mesh = 'mesh_extlarge_20m_sh005l100_topo.bin'
out_mesh = f"mesh_topo_1000vs30erf{taper}.bin"
out_mesh = f'mesh_topo_sh005l100_1000vs30erf{taper}.bin'
with open(out_mesh, 'wb') as fid:
    for iz in range(nz):
        v = np.fromfile(in_mesh, dtype='float32', count=mx * my * 3,
                        offset=mx * my * 3 * 4 * iz).reshape(my, mx, 3)
        rk = max(iz - 1, 0) / (nz - 1 - rpt)
        z0 = rk * H
        if np.min(z0) < taper:
            print(iz, np.min(z0))
            rk = max(iz, 1) / (nz - 1 - rpt)
            z1 = rk * H
            v = to_vs30_erf(v, v_surf, v_bot, vs30, taper, z0, z1, v_orig=v_orig)
        v.tofile(fid)

