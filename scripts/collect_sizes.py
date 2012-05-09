#!/usr/bin/env python

"""
collect all entries in sizes_*.dat and scale by s_min = pi*(j1/k)^2 where j1 is the first zero of bessel function J0
"""

import glob
import sys
from math import pi

if len(sys.argv) < 3:
    print 'usage: ./collect_sizes.py input_dir output_file'
    exit(-1)

j1 = 2.404825557695773 # first zero of J_0 bessel function
alpha = 0.7
scale = alpha**2 / (pi * j1**2)

with open(sys.argv[2], 'w') as outfile:
    first = True
    for f in glob.glob('%s/sizes_*.dat' % sys.argv[1]):
        with open(f) as sizefile:
            if first:
                first = False
            else:
                print >> outfile, ',',
            print >> outfile, ','.join((str(round(int(x) * scale,2)) for x in sizefile.read().split())),
