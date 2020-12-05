#! /usr/bin/env python

from glob import glob
from posix import symlink
from os.path import basename
import os

gridnum=0
tdir = os.getcwd()

flist=glob("%s/tpsrc???????" % tdir)

nf=len(flist)

for n in range(nf):
    basefile=basename(flist[n])    
    rank=int(basefile.replace("tpsrc", ""))
    filename="tpsrc_%d_%07d" % (gridnum, rank)
    symlink(flist[n], filename)
