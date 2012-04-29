#!/usr/bin/env python

import os
import re
import sys

if len(sys.argv) < 2:
    print 'usage: count_all.py path_to_sta_bin_files'
    exit(-1)

path = sys.argv[1]
escaped_path = re.sub(' ', '\ ', path)

alpha = 0.5
for f in os.listdir(path):
    match = re.match('run_(.*)\.sta_bin', f)
    if not match:
        raise ValueError('failed to parse filename: %s' % f)
    k = float(match.groups(1)[0])
    h = alpha / k
    cmd = './count -f %s/%s -l qugrs:1.0:0.4:0.7 -d %f -M 9 -u 20 -k %f -t >> counts.txt' % (escaped_path, f, h, k)
    print cmd
    os.system(cmd)
