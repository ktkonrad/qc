#!/usr/bin/env python

import re
import sys

if len(sys.argv) < 2:
    print 'usage ./parse_timing_data.py infile'
    exit(-1)

with open(sys.argv[1]) as f:
    i = 0
    for line in f:
        if i == 0:
            match = re.search('k=\[(\d+.\d+),(\d+.\d+)\]', line)
            if not match:
                raise ValueError('failed to parse line: %s' % line)
            k1 = float(match.group(1))
            k2 = float(match.group(2))
            k = (k1 + k2)/2
        elif i == 1:
            match = re.search('solving took (\d+.\d+) seconds', line)
            if not match:
                raise ValueError('failed to parse line: %s' % line)
            solve_time = float(match.group(1))
        elif i == 2:
            match = re.search('evaluating took (\d+.\d+) seconds',line)
            if not match:
                raise ValueError('failed to parse line: %s' % line)
            eval_time = float(match.group(1))
        elif i == 3:
            match = re.search('counting took (\d+.\d+) seconds', line)
            if not match:
                raise ValueError('failed to parse line: %s' % line)
            count_time = float(match.group(1))
        i += 1
        if i == 4:
            print k, solve_time, eval_time, count_time
            i = 0
            
