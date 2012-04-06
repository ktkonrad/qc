#!/usr/bin/env python

import os
import re

path = '/media/My Passport/sta'
escaped_path = '/media/My\ Passport/sta'
for f in os.listdir(path):
    match = re.match('run_(.*)\.sta_bin', f)
    k = match.groups(1)[0]
    cmd = './count -f %s/%s -l qugrs:1.0:0.4:0.7 -d 0.001 -M 9 -u 20 -k %s >> counts.txt' % (escaped_path, f, k)
    print cmd
    os.system(cmd)
