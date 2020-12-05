#!/ccs/proj/geo112/hzfmer/summit/opt/anaconda3/bin/python

import numpy as np
from numpy import sin, cos, arcsin, sqrt
import os 
import re

def distance(lon1, lat1, lon2, lat2):
    lat1, lon1, lat2, lon2 = np.radians((lat1, lon1, lat2, lon2))
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    d = 0.5 - cos(dlat) / 2 + cos(lat1) * cos(lat2)  * (1 - cos(dlon)) / 2
    return 12742 * arcsin(sqrt(d))

def mindist(npf, lon, lat, tlon, tlat):
    md = 1e13
    for i in range(len(lon)):
        adist = distance(lon[i], lat[i], tlon, tlat)
        md = min(md, adist)
    return md


params = {'6.35': (354, 212), '6.45': (144, 230), '7.35': (980, 240), '8.15': (5116, 220), '8.45': (5116, 220)}
mag = re.findall(r"\d+\.\d+", os.getcwd().split('/')[-1])[0]
npf = params[mag][0]
M = str(mag).replace('.', '_')
nx = 6320
ny = 4200

sites = np.genfromtxt('stat.txt', dtype='int', skip_header=1)
ll = np.fromfile('../scripts/surf.grid', dtype='float32').reshape(ny, nx, 2)
tll = np.zeros((len(sites), 2))
for i in range(len(sites)):
    tll[i, :] = ll[sites[i, 1] - 1, sites[i, 0] - 1, :] 

ll_fault = np.genfromtxt('fault_surf_loc.txt')
fdist = np.zeros((len(sites), ), dtype='float32')

for j in range(len(sites)):
    fdist[j] = mindist(npf, ll_fault[:, 0], ll_fault[:, 1], tll[j, 0], tll[j, 1])

fdist.tofile(f'site_rjb_M{M}.bin')
