#!/usr/bin/env python

import matplotlib
matplotlib.use("AGG", warn=False)
from pylab import *
import sys
import os
import struct

'''struct seisheader {
  char version[8];
  char site_name[8];
  //in case we think of something later
  char padding[8];
  int source_id;
  int rupture_id;
  int rup_var_id;
  float dt;
  int nt;
  int comps;
  float det_max_freq;
  float stoch_max_freq;
  };'''

class Seismogram:
	
	def __init__(self, filename, label):
		self.filename = filename
		self.label = label
		self.nt = 0
		self.dt = 0.0

	def createTimesteps(self):
		self.timesteps = []
		for i in range(0, self.nt):
			self.timesteps.append(i*self.dt)

	def parseHeader(self, header_str):
		#dt is bytes 36-40, nt is 40-44
		self.dt = struct.unpack("f", header_str[36:40])[0]
		self.nt = struct.unpack("i", header_str[40:44])[0]
		print "Using nt=%d, dt=%f for file %s" % (self.nt, self.dt, self.filename)

	def readData(self):
		fp_in = open(self.filename, "r")
		header_str = fp_in.read(56)
		self.parseHeader(header_str)
		data_str = fp_in.read(4*self.nt)
		self.x_data = struct.unpack("%df" % self.nt, data_str)
		data_str = fp_in.read(4*self.nt)
		self.y_data = struct.unpack("%df" % self.nt, data_str)
		fp_in.close()


if len(sys.argv)<5:
	print "Usage: %s <seismogram 1> <seis 1 plot label> <seismogram 2> <seis 2 plot label> ... <seismogram N> <seis N plot label> <PNG filename> <plot title>" % (sys.argv[0])
	sys.exit(1)

seismograms = []

for i in range(1, len(sys.argv)-2, 2):
	seismograms.append(Seismogram(sys.argv[i], sys.argv[i+1]))

output_filename = sys.argv[-2]
plot_title = sys.argv[-1]

num_seis = len(seismograms)

max_y = 0.0
min_y = 0.0

for seis in seismograms:
	seis.readData()
	seis.createTimesteps()
	max_y = max([max_y, max(seis.x_data), max(seis.y_data)])
	min_y = min([min_y, min(seis.x_data), min(seis.y_data)])

max_y = 1.1*max_y
min_y = 1.1*min_y


clf()
subplot(211, title="X component (cm/s)")
for seis in seismograms:
	plot(seis.timesteps, seis.x_data, label=seis.label)
#xlim(0, seis.timesteps[-1])
xlim(0, 200)
#ylim(-80, 80)
ylim(min_y, max_y)
legend(loc="upper left", prop={'size': 10})
subplot(212, title="Y component (cm/s)")
for seis in seismograms:
        plot(seis.timesteps, seis.y_data, label=seis.label)
#xlim(0, seis.timesteps[-1])
xlim(0, 200)
#ylim(-80, 80)
ylim(min_y, max_y)
legend(loc="upper left", prop={'size': 10})
suptitle(plot_title)
gcf().set_size_inches(14, 7)
savefig(output_filename, format="png")

