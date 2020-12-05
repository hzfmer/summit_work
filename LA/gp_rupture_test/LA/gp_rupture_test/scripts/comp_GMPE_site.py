#!/ccs/proj/geo112/hzfmer/summit/opt/anaconda3/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue March  19 18:37:05 2019

@author: zhh076
"""

import numpy as np
import sys
import os
sys.path.append("/ccs/home/hzfmer/scratch/LA/gp_rupture_test/LA/gp_rupture_test/scripts")
import BSSA_2014_nga
           

nx = 6320
ny = 4200
n = nx * ny
for mag in [6.35, 6.45, 7.35, 8.15, 8.45]:
    M = str(mag).replace('.', '_')
    data_gmpe = np.fromfile(f'site_info_for_GMPE_M{M}.bin',
            dtype='float32').reshape((-1, 4))
    dist, z1, vs30, vs0 = [data_gmpe[:, i] for i in range(data_gmpe.shape[1])]
    nstat = len(dist)
    for T in [1, 2, 3, 5]:
        sa_median_gmpe, sa_sigma_gmpe, sa_median_gmpe0, sa_sigma_gmpe0 = \
                [np.zeros((nstat, 1), dtype='float32') for _ in range(4)]
        for i in range(nstat):
            sa_median_gmpe[i], sa_sigma_gmpe[i], _ = BSSA_2014_nga.BSSA_2014_nga(  \
                        mag, T, dist[i], 1, 1, z1[i], vs30[i])
            sa_median_gmpe0[i], sa_sigma_gmpe0[i], _ = BSSA_2014_nga.BSSA_2014_nga(  \
                        mag, T, dist[i], 1, 1, z1[i], vs0[i])
        output = np.column_stack((dist, sa_median_gmpe, sa_sigma_gmpe))
        output.tofile(f'sa_GMPE_site_M{M}_T{T}.bin')
        output = np.column_stack((dist, sa_median_gmpe0, sa_sigma_gmpe0))
        output.tofile(f'sa_GMPE_vs0_site_M{M}_T{T}.bin')

