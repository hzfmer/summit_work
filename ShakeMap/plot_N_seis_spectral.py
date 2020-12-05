#!/usr/bin/env python

import sys
import os
import numpy
import matplotlib
import struct
from pylab import *

if len(sys.argv)<3 or len(sys.argv)%2!=0:
	print "Usage %s <seis 1> <label 1> <seis 2> <label 2> ... <seis N> <label N> <output file>" % sys.argv[0]
	sys.exit(1)

in_files = []
labels = []

for i in range(1, len(sys.argv)-1, 2):
	in_files.append(sys.argv[i])
	labels.append(sys.argv[i+1])

output_file = sys.argv[-1]
freqs = []
data_x = []
data_y = []

MAX_FREQ = 0.0

for i in range(0, len(in_files)):
	fp_in = open(in_files[i], "rb")
	header = fp_in.read(56)
	dt = struct.unpack("f", header[36:40])[0]
	nt = struct.unpack("i", header[40:44])[0]
	data_str = fp_in.read(nt*4)
	data_vals_x = struct.unpack("%df" % (nt), data_str)
	data_str = fp_in.read(nt*4)
	data_vals_y = struct.unpack("%df" % (nt), data_str)
	arr = numpy.array(data_vals_x)
	N = nt
	print "Using nt=%d, dt=%f" % (nt, dt)
	MAX_FREQ = nt*dt
	fft_x = numpy.fft.fft(data_vals_x)
	#freq = numpy.linspace(0.0, MAX_FREQ, num=MAX_FREQ/dt)
	freq = numpy.linspace(0.0, 1.0/(2.0*dt), num=nt/2)
	fft_y = numpy.fft.fft(data_vals_y)
	max_ind = int(MAX_FREQ/dt)
	freqs.append(freq)
	#data_x.append([abs(x)*dt for x in fft_x[0:max_ind]])
	#data_y.append([abs(x)*dt for x in fft_y[0:max_ind]])
	data_x.append([2.0/nt*abs(x) for x in fft_x[0:nt/2]])
	data_y.append([2.0/nt*abs(x) for x in fft_y[0:nt/2]])
clf()
subplots_adjust(hspace=0.25, left=0.09, right=0.96, top=0.96, bottom=0.08)
subplot(211, title="X component")
for i in range(0, len(freqs)):
	print freqs[i][-1]
	plot(freqs[i], data_x[i], label=labels[i])
loglog()
xlim(0.001, 10.0)
ylim(1e-8, 1.0)
xlabel("Frequency (Hz)")
legend(loc="upper right")
subplot(212, title="Y component")
for i in range(0, len(freqs)):
	plot(freqs[i], data_y[i], label=labels[i])
loglog()
xlim(0.001, 10.0)
ylim(1e-8, 1.0)
xlabel("Frequency (Hz)")
legend(loc="upper right")
gcf().set_size_inches(8, 8)
savefig(output_file, format="png")


