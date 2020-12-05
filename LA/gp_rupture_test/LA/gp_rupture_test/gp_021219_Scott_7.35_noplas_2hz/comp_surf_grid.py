#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np


data = np.loadtxt('../../data_for_Zhifeng/mesh/model_coords_GC_USC')
lon = data[:, 0].reshape(6320, 4200).T
lat = data[:, 1].reshape(6320, 4200).T
out = np.float64(np.vstack([lon.flatten(), lat.flatten()]).T)

out.tofile('surf.grid')
