#!/sw/xk6/python_anaconda3/4.2.0/sles11.3_gnu4.9.3/bin/python

import numpy as np
#from scipy.interpolate import NearestNDInterpolator
from sys import stdout
import urllib
import urllib.request
import re
import os

def get_stat(html):
    reg = r'href="(.*).bbp"'
    stat_re = re.compile(reg)
    stat_list = stat_re.findall(html.decode('UTF-8'))
    return stat_list

base_url = 'http://hypocenter.usc.edu/research/highf_data/rwg_20180425/Vel-bbp/'
response = urllib.request.urlopen(base_url)
html = response.read()
stat_list = get_stat(html)

nx = 9000
ny = 6750
npt = nx * ny
n = 0
#alli = np.zeros((npt,), dtype=int)
#allj = np.zeros((npt,), dtype=int)
#ll = np.fromfile('/lustre/atlas/proj-shared/geo112/huzf/high_f/cvm/la_habra_large_cvmsi_20m.grid').reshape(npt,3)
#for k in range(ny):
#    stdout.write("\rDefining interpolating data ({0}/{1})".format(k, ny))
#    stdout.flush()
#    for l in range(nx):
#        alli[n] = l + 1
#        allj[n] = k + 1
#        n += 1
#interpi = NearestNDInterpolator(ll[:,0:2], alli)
#interpj = NearestNDInterpolator(ll[:,0:2], allj)
# 

f_stat = open('la_habra_large_statlist.txt', 'w')
#f_statidx = open('la_habra_large_stat.idx', 'w')
for F in stat_list:
    stdout.write("\rWriting station {0}\n".format(F))
    stdout.flush()
    dat = urllib.request.urlopen(base_url + '/{0}.bbp'.format(F)).readlines()
    f_bbp = open('{0}.bbp'.format(F), 'w')
    for i in range(len(dat)):
        line = dat[i].decode('UTF-8')
        if 'lat' in line:
            lat = re.findall(r'-*\d+.\d+', line)[0]
        elif 'lon' in line:
            lon = re.findall(r'-*\d+.\d+', line)[0]
        if '#' not in line:
            break
    f_stat.write('{0} {1} {2}\n'.format(F, lon, lat))
#    f_statidx.write('{0} {1}'.format(interpi([lon, lat]), interpj([lon, lat])))

    for j in range(i, len(dat)):
        line = dat[j].decode('UTF-8').split('\t')[1:]  # Ignore timestamp, which starts at -6s
        f_bbp.writelines(' '.join(line))

f_bbp.close()
f_stat.close()
