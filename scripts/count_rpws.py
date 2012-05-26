#!/usr/bin/env python

import glob
import sys
import re
import os

if len(sys.argv) != 2:
    print 'usage: ./count_rpws.py indir'
    exit(-1)

indir = sys.argv[1]

no_interp_outfile = 'rpw_counts_no_interp.txt'

# remove old output
try:
    os.remove(no_interp_outfile);
except OSError:
    pass

for f in glob.glob('%s/rpw_*.sta_bin' % indir):
    m = re.search('rpw_([\d\.]+).sta_bin', f);
    if not m:
        raise ValueError('Failed to parse filename: %s' % f)
    dx = float(m.group(1));
    cmd = '../c/count -f %s -d %.8f -k 1.0 -n -q >> %s' % (f, dx, no_interp_outfile);
    os.system(cmd)

# sort the output
os.system('cat %s | sort -nr > temp; mv temp %s' % (no_interp_outfile,no_interp_outfile))

