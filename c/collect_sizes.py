#!/usr/bin/env python

"""
collect all entries in sizes_*.dat and scale by s_min = pi*(j1/k)^2 where j1 is the first zero of bessel function J0
"""

import glob
import sys
from math import pi

if len(sys.argv) < 2:
    print 'usage: ./collect_sizes.py output_file'
    exit(-1)

j1 = 2.404825557695773 # first zero of J_0 bessel function
alpha = 0.5
scale = alpha**2 / (pi * j1**2)

with open(sys.argv[1], 'w') as outfile:
    for f in glob.glob('sizes_*.dat'):
        with open(f) as sizefile:
            print >> outfile, ','.join((str(round(int(x) * scale,2)) for x in sizefile.read().split()))
