#!/sw/titan/python/anaconda3-5.1.0/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 11 13:56:23 2018

@author: zhh076
"""

import numpy as np
import pandas
import csv
import utm
import os
import struct
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from math import cos, asin, sqrt
from scipy.interpolate import NearestNDInterpolator
from scipy.interpolate import UnivariateSpline
from scipy.ndimage import zoom
#from scipy.ndimage import gaussian_filter
from jpmesh import Angle, Coordinate, ThirdMesh

plt.style.use('seaborn-talk')
plt.tight_layout()
plt.close('all')

def add_bd_vs(vs, bd):
    # vs at bedrock depth = 1100
    bore_vs = [140,140,430,660,770,640,1100]
    dx = dh * 1000
    print(f"dx in add_bd_vs is: {dx}")
    dim0, dim1, dim2 = vs.shape
    for i in range(dim0):
        for j in range(dim1):
            bore_dp = np.array([5,7,20,30,80,np.max([80, bd[i,j]])])
            for k in range(dim2):      
                if vs[i,j,k]>1100:
                    if k*dx>bd[i,j]:
                        break
                    else:
                        vs[i,j,k] = 1100
                else:
                    if k*dx>=bd[i,j]:
                        vs[i,j,k] = 1100
                    else:
                        depth = np.searchsorted(bore_dp,k*dx)
                        vs[i,j,k] = bore_vs[depth]
    return vs
                
def add_bd_vp(vp, vs, layer, bd):   
    # vs at bedrock depth = 1100
    interp_vp = get_phy(inperp_vp=True)
    vp_bed = float(interp_vp(1100))
    bore_vp = [280,1760,1760,1760,2090,2090,vp_bed]
    dx = dh * 1000
    print(f"dx in add_bd_vp is: {dx}")
    dim0, dim1, dim2 = vs.shape
    for i in range(dim0):
        for j in range(dim1):
            bore_dp = np.array([5,7,20,30,80,np.max([80, bd[i,j]])])
            for k in range(dim2):
                if vs[i,j,k]>1100:
                    if (k + layer) * dx > bd[i,j]:
                        break
                    else:
                        vp[i,j,k] = vp_bed
                else:
                    if (k + layer) * dx >= bd[i,j]:
                        vp[i,j,k] = vp_bed
                    else:
                        depth = np.searchsorted(bore_dp,k*dx)
                        vp[i,j,k] = bore_vp[depth]
    return vp                
                    
            
def get_phy(*vs, **kwargs):
    with open(os.path.join(data_dir,"D_V2_STRUCT_DEEP_PYS.csv")) as f_phy:
        tmp = list(csv.reader(f_phy))
        param = np.asarray(tmp[8:], dtype=int)
    param = np.delete(param, 0, 1)  
    # just pick up default vs-profile
    if len(vs)==0 and len(kwargs)==0:
        return param[:,1]
          
    param = np.vstack({tuple(row) for row in param})  
    param = param[param[:,1].argsort()]         # remove duplicates
    
    svp = param[:, 0] 
    svs = param[:, 1]
    sro = param[:, 2]
    sqp = param[:, 3]
    sqs = param[:, 4]
    
    interp_vp = UnivariateSpline(np.append(140,svs), np.append(280,svp), k=1)
                    # specifically dealing due to borehole data at TKCH05
    interp_ro = UnivariateSpline(svs, sro, k=1)
    interp_qp = UnivariateSpline(svs, sqp, k=1)
    interp_qs = UnivariateSpline(svs, sqs, k=1)
    
    if len(kwargs) > 0:
        for key in kwargs:
            if kwargs[key] is True:
                return interp_vp
    
    vp = np.squeeze(interp_vp(vs))
    ro = np.squeeze(interp_ro(vs))
    qp = np.squeeze(interp_qp(vs))
    qs = np.squeeze(interp_qs(vs))
    return vp, ro, qp, qs
 
def get_vs(vs_depth):
    dx = dh * 1000
    # len(vs_profile) = len(vs_depth) + 1, indicating 
    # the velocity below last layer
    vs_profile = get_phy()
    vs = np.zeros((nz,), dtype=float)
    for i in range(nz):
        # nz-1, make sure the last point is on the surface
        if dx*(nz-1-i)<vs_depth[-1]:
            vs_depth = vs_depth[vs_depth < dx*(nz-1-i)]
            vs_profile = vs_profile[0:(len(vs_depth)+1)]
        vs[i] = vs_profile[-1]
    return np.flipud(vs)

def distance(lat1, lon1, lat2, lon2):
    p = np.pi / 180;
    d = 0.5 - cos((lat2 - lat1) * p) / 2 + cos(lon1 * p) * cos(lon2 * p) \
        * (1 - cos((lon2 - lon1) * p)) / 2
    return 12742 * asin(sqrt(d))   

def get_index(latlon):
    idx = idx_interp(latlon) 
    (id_y, id_x) = np.unravel_index(idx, (ny, nx))
    return [id_y, id_x]
   
    
nz = 800
dh = 5e-3
# threshold of MRVBF 
MR = 1.5

data_dir = '/lustre/atlas/proj-shared/geo112/huzf/preevents/ssh'
res_dir = '/lustre/atlas/proj-shared/geo112/huzf/preevents/ssh'

# MRVBF for each lat/lon      
f = open(os.path.join(data_dir, 'MRVBF_honbetsu.asc') )
north = int(f.readline().split(' ')[-1].split('\n')[0])
south = int(f.readline().split(' ')[-1].split('\n')[0])
east = int(f.readline().split(' ')[-1].split('\n')[0])
west = int(f.readline().split(' ')[-1].split('\n')[0])
ny = int(f.readline().split(' ')[-1].split('\n')[0])
nx = int(f.readline().split(' ')[-1].split('\n')[0])

#mrvbf = np.zeros((ny,nx))    
#for i in range(ny):
#    dat =f.readline()   
#    dat = dat.replace('\n', ' ')
#    mrvbf[i,:] = dat.split(' ')[0:nx]  
f.close()

# reducing the spatial step by 2, 08/13/18  
scale = 2.0
#mrvbf = zoom(mrvbf, scale, order=0)
nx, ny, nz = (int(x * scale) for x in (nx, ny, nz))
dh = dh / scale

lat = np.zeros((ny, nx))
lon = np.zeros((ny, nx))
ver = np.linspace(south, north, ny, dtype='i')
hor = np.linspace(west, east, nx, dtype='i')

for i in range(ny):
    for j in range(nx):
        lat[i,j] = utm.to_latlon(hor[j], ver[i], 54, 'T')[0]
        lon[i,j] = utm.to_latlon(hor[j], ver[i], 54, 'T')[1]
ll = np.asarray([lat.flatten(), lon.flatten()]).T
idx_interp = NearestNDInterpolator(ll,np.arange(nx*ny))

#interp_mrvbf = NearestNDInterpolator(ll, np.ravel(mrvbf))     
#surface_site = [43.127773, 143.615445]
site = [43.121311, 143.618352]
site_idx = get_index(site)
#site_mrvbf = mrvbf[site_idx]

# Artificial bedrock-depth MRVBF relation
# MRVBF = A * (log(x/1.5)) ** B
# assuming site_depth_to_bedrock = 120
#outcrop = mrvbf * (mrvbf>=MR)
#outcrop[outcrop==0] = MR
site_bd = 120
#A = 108
#B = np.log(site_bd/A)/np.log(np.log(site_mrvbf/MR))
#bed_depth = A * (np.log(outcrop / MR)) ** B
bed_depth = site_bd * np.ones((ny, nx))
#bed_depth = gaussian_filter(bed_depth, [3,3])


# Background structure
vs_bg = pandas.read_csv(os.path.join(data_dir,"D_V2_STRUCT_DEEP_LYRD_6443N.csv"), 
                        skiprows=6,index_col=0)
coordinate = Coordinate(lon=Angle.from_degree(lon[i,j]),
                        lat=Angle.from_degree(lat[i,j]))
loc_code = ThirdMesh.from_coordinate(coordinate).code 
tmp = get_vs(vs_bg.loc[loc_code+'N'].values[1:])
vs = np.tile(tmp, [ny, nx, 1])
#vs_database = {}
#vs = np.zeros((ny, nx, nz))
#for i in range(ny):
#    for j in range(nx):
#        coordinate = Coordinate(lon=Angle.from_degree(lon[i,j]),
#                                lat=Angle.from_degree(lat[i,j]))
#        loc_code = ThirdMesh.from_coordinate(coordinate).code 
#        if loc_code not in vs_database:   
#            vs[i,j,:] = get_vs(vs_bg.loc[loc_code+'N'].values[1:])
#            vs_database[loc_code] = vs[i,j,:]
#        else:
#            vs[i,j,:] = vs_database[loc_code]
#         

# plot original VS
# South-North cross section
fig, ax = plt.subplots(2,1,figsize=(5,6), sharex=True,
                       gridspec_kw=dict(width_ratios=[1], height_ratios=[1,5]))
plt.subplots_adjust(left=0.15, right=0.85)
cax1 = ax[1].imshow(np.squeeze(vs[:,site_idx[1],:]).T,cmap='coolwarm_r',
                interpolation='none',extent=[0, 2.4, 4, 0])
cb = plt.colorbar(cax1, ax=ax[1], orientation='vertical', shrink=0.7)
#ax[1].scatter(site_idx[0]*dh, sz*dh,s=300, c='g', marker=(5,2,0))
cb.set_label('Vs (m/s)')
ax[1].set_xlabel('South-North (km)')
ax[1].set_ylabel('Depth (km)')

top_vs = np.squeeze(vs[:, site_idx[1], 0:(nz//8)]).T
cax0 = ax[0].imshow(top_vs, interpolation='none',extent=[0, 2.4, 0.5, 0])
cb = plt.colorbar(cax0, ax=ax[0],ticks=np.linspace(top_vs.min(), top_vs.max(), 3),
                  orientation='vertical', shrink=1.5,aspect=9)
#ax[0].set_ylabel('Depth (km)')
fig.show()
fig.savefig(os.path.join(res_dir,'vs_perpend_orig_TKCH05_1d.svg'), dpi=600,
                    bbox_inches='tight', pad_inches=0.1)

# West-East cross section
fig, ax = plt.subplots(2,1,figsize=(6,6), sharex=True,
                       gridspec_kw=dict(width_ratios=[1], height_ratios=[1,4]))
plt.subplots_adjust(left=0.15, right=0.9)
cax1 = ax[1].imshow(np.squeeze(vs[site_idx[0],:,:]).T,cmap='coolwarm_r',
                interpolation='none',extent=[0, 5, 4, 0])
cb = plt.colorbar(cax1, ax=ax[1], orientation='vertical', shrink=0.7)
#ax[1].scatter(site_idx[1]*dh, sz*dh,s=300, c='g', marker=(5,2,0))
cb.set_label('Vs (m/s)')
ax[1].set_xlabel('West-East (km)')
ax[1].set_ylabel('Depth (km)')

top_vs = np.squeeze(vs[site_idx[0],:,0:(nz//8)]).T
cax0 = ax[0].imshow(top_vs, interpolation='none',extent=[0, 5, 0.5, 0])
cb = plt.colorbar(cax0, ax=ax[0],ticks=np.linspace(top_vs.min(), top_vs.max(), 3),
                  orientation='vertical', shrink=1.5,aspect=9)
#ax[0].set_ylabel('Depth (km)')
fig.show()
fig.savefig(os.path.join(res_dir,'vs_horizon_orig_TKCH05_1d.svg'), dpi=600,
                    bbox_inches='tight', pad_inches=0.1)

## To add more constraints 

vs = add_bd_vs(vs, bed_depth)        
#vp, ro, qp, qs = get_phy(vs)  
#vp = add_bd_vp(vp, vs, bed_depth) 

# Write the mesh file
#f_out = open(os.path.join(data_dir,'mesh_TKCH05_061818.bin'), 'wb')
#for k in range(nz):
#    vp, ro, qp, qs = get_phy(vs[:,:,k])
#    vp = add_bd_vp(vp, vs[:,:,k], k, bed_depth)
#    for i in range(ny-1,-1,-1):      # make sure start from left lower corner
#        for j in range(nx):
#            f_out.write(struct.pack('5f',vp[i,j],vs[i,j,k],ro[i,j],
#                                    qp[i,j],qs[i,j]))
#f_out.close()
          
# flip: make sure start from left lower corner
step = np.linspace(0, nz, 11).astype(int)
for ii in range(10):
    print(f"Saving layer{ii}")
    tmpvs = vs[:,:,step[ii]:step[ii+1]]
    vp, ro, qp, qs = get_phy(tmpvs)
    vp = add_bd_vp(vp, tmpvs, ii * (nz // 10), bed_depth)
    dat = np.ravel(list(map(lambda p: np.rollaxis(p[::-1,:,:].swapaxes(1,2), 1),
                            (vp, tmpvs, ro, qp, qs))))
    dat = dat.reshape(5,-1).T.astype('f')
    dat.tofile(os.path.join(data_dir,f'mesh_TKCH05_1d_081318_{ii}.bin'))

# Plot whatever resulted

## MRVBF
#fig, ax = plt.subplots(figsize=(8,6))
#plt.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.1)
#cax = ax.contourf(lon, lat, np.flipud(mrvbf),cmap='RdBu')
#ax.scatter(site[1], site[0], s=300, c='g', marker='^')
#cb = plt.colorbar(cax, orientation='horizontal', shrink=0.7)
#cb.ax.set_xlabel('MRVBF')
#fig.show()
#fig.savefig(os.path.join(res_dir,'MRVBF_TKCH05_1d.svg'), dpi=600,
#                    bbox_inches='tight', pad_inches=0.1)
#
## MRVBF --- Depth to bedrock
#fig, ax = plt.subplots(figsize=(6,6))
#plt_x = np.append(np.linspace(0.1,MR,60),np.linspace(MR,6,180))
#plt_y = np.append(np.zeros(60),
#                  A * (np.log(np.linspace(MR, 6,180) / MR)) ** B)
#ax.plot(plt_x, plt_y)
#ax.scatter(site_mrvbf, site_bd ,s=300, c='r', marker=(5,2,0))
#ax.set_xticks(np.linspace(0,6,13))
#ax.grid(which='both',linestyle='--')
#ax.set_xlabel('MRVBF')
#ax.set_ylabel('Depth to Bedrock (m)')
#ax.set_title("Let Depth to Bedrock = 120 m at TKCH05\n")
#fig.savefig(os.path.join(res_dir,'MRVBF_bd_TKCH05_1d.svg'), dpi=600,
#                    bbox_inches='tight', pad_inches=0.1)
#
## Depth to bedrock
fig, ax = plt.subplots(figsize=(8,6))
plt.subplots_adjust(left=0.15, right=0.95, top=0.9, bottom=0.1)
cax = ax.contourf(lon, lat, np.flipud(bed_depth),cmap='RdBu')
ax.scatter(site[1], site[0], s=300, c='g', marker='^')
cb = plt.colorbar(cax, orientation='horizontal', shrink=0.7)
cb.ax.set_xlabel('Depth to bedrock (m)')
fig.show()
fig.savefig(os.path.join(res_dir,'bed_depth_TKCH05_1d.svg'), dpi=600,
                    bbox_inches='tight', pad_inches=0.1)
     
# South-North cross section      
fig, ax = plt.subplots(2,1,figsize=(5,6), sharex=True,
                       gridspec_kw=dict(width_ratios=[1], height_ratios=[1,5]))
plt.subplots_adjust(left=0.15, right=0.85)
cax1 = ax[1].imshow(np.squeeze(vs[:,site_idx[1],:]).T,cmap='coolwarm_r',
                interpolation='none',extent=[0, 2.4, 4, 0])
cb = plt.colorbar(cax1, ax=ax[1], orientation='vertical', shrink=0.7)
# ax[1].scatter(site_idx[0]*dh, sz*dh,s=300, c='g', marker=(5,2,0))
cb.set_label('Vs (m/s)')
ax[1].set_xlabel('South-North (km)')
ax[1].set_ylabel('Depth (km)')

top_vs = np.squeeze(vs[:, site_idx[1], 0:(nz//8)]).T
cax0 = ax[0].imshow(top_vs, interpolation='none',extent=[0, 2.4, 0.5, 0])
cb = plt.colorbar(cax0, ax=ax[0],ticks=np.linspace(top_vs.min(), top_vs.max(), 3),
                  orientation='vertical', shrink=1.5,aspect=9)
#ax[0].set_ylabel('Depth (km)')
fig.show()
fig.savefig(os.path.join(res_dir,'vs_perpend_TKCH05_1d.svg'), dpi=600,
                    bbox_inches='tight', pad_inches=0.1)


# West-East cross section
fig, ax = plt.subplots(2,1,figsize=(6,6), sharex=True,
                       gridspec_kw=dict(width_ratios=[1], height_ratios=[1,4]))
plt.subplots_adjust(left=0.15, right=0.9)
cax1 = ax[1].imshow(np.squeeze(vs[site_idx[0],:,:]).T,cmap='coolwarm_r',
                interpolation='none',extent=[0, 5, 4, 0])
cb = plt.colorbar(cax1, ax=ax[1], orientation='vertical', shrink=0.7)
#ax[1].scatter(site_idx[1]*dh, sz*dh,s=300, c='g', marker=(5,2,0))
cb.set_label('Vs (m/s)')
ax[1].set_xlabel('West-East (km)')
ax[1].set_ylabel('Depth (km)')

top_vs = np.squeeze(vs[site_idx[0],:,0:(nz//8)]).T
cax0 = ax[0].imshow(top_vs, interpolation='none',extent=[0, 5, 0.5, 0])
cb = plt.colorbar(cax0, ax=ax[0],ticks=np.linspace(top_vs.min(), top_vs.max(), 3),
                  orientation='vertical', shrink=1.5,aspect=9)
#ax[0].set_ylabel('Depth (km)')
fig.show()
fig.savefig(os.path.join(res_dir,'vs_horizon_TKCH05_1d.svg'), dpi=600,
                    bbox_inches='tight', pad_inches=0.1)

#plt.figure()
#plt.pcolor(lon, lat, np.flipud(mrvbf>=threshold).astype(int),cmap='RdBu')
#cb = plt.colorbar(orientation='horizontal', shrink=0.7)
#plt.scatter(site[1], site[0], s=100, c='g', marker='^')

