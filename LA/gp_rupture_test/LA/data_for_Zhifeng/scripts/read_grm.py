#!/usr/bin/env python2

import sys
import struct

if len(sys.argv)<2:
   print "Usage: %s <input seismogram>" % sys.argv[0]
   sys.exit(1)

seismogram = sys.argv[1]

with open(seismogram, "rb") as fp_in:
   header_str = fp_in.read(56)
   version = header_str[0:8]
   if version[0:5]!="12.10":
       print "Error: version does not match expected string '12.10', aborting."
       sys.exit(2)
   site = header_str[8:16]
   source_id = struct.unpack('i', header_str[24:28])
   rupture_id = struct.unpack('i', header_str[28:32])
   rup_var_id = struct.unpack('i', header_str[32:36])
   dt = struct.unpack('f', header_str[36:40])
   nt = struct.unpack('i', header_str[40:44])
   comps = struct.unpack('i', header_str[44:48])[0]
   print(comps)
   x_flag = y_flag = z_flag = False
   if (comps & 1)==1:
       x_flag = True
   if (comps & 2)==2:
       y_flag = True
   if (comps & 4)==4:
       z_flag = True
   det_max_freq = struct.unpack('f', header_str[48:52])[0]
   stoch_max_freq = struct.unpack('f', header_str[52:56])[0]
   print "Version = %s" % version
   print "Site = %s" % site
   print "Source ID = %d" % source_id
   print "Rupture ID = %d" % rupture_id
   print "Rupture Variation ID = %d" % rup_var_id
   print "DT = %f" % dt
   print "NT = %d" % nt
   print "X component? %d" % x_flag
   print "Y component? %d" % y_flag
   print "Z component? %d" % z_flag
   print "Maximum deterministic frequency = %f" % det_max_freq
   print "Maximum stochastic frequency = %f" % stoch_max_freq
   fp_in.close()
