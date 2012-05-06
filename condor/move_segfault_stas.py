#!/usr/bin/env python

import os
import glob
import sys

if len(sys.argv) < 2:
    print 'usage: ./move_segfault_stas.py output_directory'
    exit(-1)

outdir = sys.argv[1]

for log_filename in glob.glob('*.log'):
    with open(log_filename) as f:
        for line in f:
            if 'signal 11' in line:
                sta_filename = log_filename[:-3] + 'sta_bin'
                os.rename(sta_filename, '%s/%s' % (outdir, sta_filename))
