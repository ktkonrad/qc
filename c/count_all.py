#!/usr/bin/env python

import os
import re

alpha = 0.5

path = '/media/My Passport/sta'
escaped_path = '/media/My\ Passport/sta'
for f in os.listdir(path):
    match = re.match('run_(.*)\.sta_bin', f)
    k = match.groups(1)[0]
    h = alpha/k
    cmd = './count -f %s/%s -l qugrs:1.0:0.4:0.7 -d %f -M 9 -u 20 -k %s >> counts.txt' % (escaped_path, f, h, k)
    print cmd
    os.system(cmd)
