#!/ccs/home/hzfmer/file_back/programs/anaconda3/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 13:31:19 2019

@author: zhh076

data format:
    X, oriented north
    Y, oriented east
    synthetics: velocity (cm/s)
"""
import numpy as np
# import matplotlib.pyplot as plt
import mysql.connector
import os
import struct
# import shelve
# from collections import OrderedDict

home_dir = '/gpfs/alpine/geo112/scratch/hzfmer/LA/gp_rupture_test/LA/data_for_Zhifeng/'
os.chdir(home_dir)


conn = mysql.connector.connect(host='focal.usc.edu', user='cybershk_ro',
                               password='CyberShake2007',
                               database='CyberShake',
                               use_unicode=True)
cursor = conn.cursor()
cursor.execute(f'select * from CyberShake_Sites')
site_tmp = cursor.fetchall()
site_info = {}
for i in range(len(site_tmp)):
    site_info[site_tmp[i][2]] = [float(site_tmp[i][3]), float(site_tmp[i][4])]

try:
    ff = np.load('null')
    seis6 = np.load('new_data_15_4/seis6_45.npy')[()]
    seis8 = np.load('new_data_15_4/seis8_15.npy')[()]
except FileNotFoundError:
    seis_dir = home_dir + 'new_seismograms_15_4/'
    f_site6 = open('new_data_15_4/site_names6_45.txt', 'w')
    f_site8 = open('new_data_15_4/site_names8_15.txt', 'w')
    site_names6 = []
    site_names8 = []
    seis_list = sorted(os.listdir(seis_dir))
#   seis_shelf = shelve.open(home_dir + 'data_15_12/seis_Scott_SAF_3sce.shlf')

    seis6 = {}
    seis8 = {}
    ll6 = []
    ll8 = []
    for (i, F) in enumerate(seis_list):
        with open(seis_dir + F, "rb") as fp_in:
            header_str = fp_in.read(56)
            site = header_str[8:16].decode('utf-8', 'ignore')
            site = site.split(sep='\x00', maxsplit=1)[0]
#            try:
#                site = header_str[8:16].decode('utf-8').strip('\x00')
#            except:
#                print(header_str[8:16])
#                print(site)
            source_id = struct.unpack('i', header_str[24:28])[0]
            rupture_id = struct.unpack('i', header_str[28:32])[0]
            rup_var_id = struct.unpack('i', header_str[32:36])[0]
            dt = struct.unpack('f', header_str[36:40])
            nt = struct.unpack('i', header_str[40:44])[0]
            comps = struct.unpack('i', header_str[44:48])[0]
#            print(comps)
#            x_flag = y_flag = z_flag = False
#            if (comps & 1)==1:
#                x_flag = True
#            if (comps & 2)==2:
#                y_flag = True
#            if (comps & 4)==4:
#                z_flag = True
            tmp = struct.unpack(f'{nt * 2}f', fp_in.read(nt * 8))
            if source_id == 68:
                site_names8.append(site)
                seis8[site] = np.array(tmp).reshape((2, -1))
                ll8.append([site_info[site][1], site_info[site][0]])
            if source_id == 262:
                site_names6.append(site)
                seis6[site] = np.array(tmp).reshape((2, -1))
                ll6.append([site_info[site][1], site_info[site][0]])

np.save('new_data_15_4/seis6_45.npy', seis6)
np.save('new_data_15_4/seis8_15.npy', seis8)
np.savetxt('new_data_15_4/latlon6_45.txt', np.array(ll6), fmt='%f',
           delimiter=' ', newline='\n')
np.savetxt('new_data_15_4/latlon8_15.txt', np.array(ll8), fmt='%f',
           delimiter=' ', newline='\n')
for st in site_names6:
    f_site6.write(st + '\n')
f_site6.close()
for st in site_names8:
    f_site8.write(st + '\n')
f_site8.close()
