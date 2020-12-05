#! /usr/bin/env python

from glob import glob
from posix import symlink
from os.path import basename
import os

gridnum=0
tdir = os.getcwd()
flist=glob("%s/fault*" % tdir)

nf=len(flist)

for n in range(nf):
    basefile=basename(flist[n]).strip('fault')
    rank=int(basefile.split("_")[0])
    filenum=int(basefile.split("_")[1])
    filename="fault_%d_%07d_%03d" % (gridnum, rank, filenum)
    symlink(flist[n], filename)
