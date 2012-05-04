#!/usr/bin/env python

import os
import re
import sys

if len(sys.argv) < 3:
    print 'usage: count_all.py path_to_sta_bin_files output_file [already_counted_file]'
    exit(-1)

path = sys.argv[1]
escaped_path = re.sub(' ', '\ ', path)
outfile = sys.argv[2]
already_counted_ks = []

if len(sys.argv) >= 4:
    with open(sys.argv[3]) as f:
        for line in f:
            already_counted_ks.append(float(line.split(',')[0]))

alpha = 0.7
for f in os.listdir(path):
    match = re.match('run_(.*)\.sta_bin', f)
    if not match:
        raise ValueError('failed to parse filename: %s' % f)
    k = float(match.groups(1)[0])
    if k in already_counted_ks:
        print 'already counted %f' % k
        continue
    h = alpha / k
    cmd = './count -f %s/%s -l qugrs:1.0:0.4:0.7 -d %f -M 9 -u 20 -k %f >> %s' % (escaped_path, f, h, k, outfile)
    print cmd
    os.system(cmd)
