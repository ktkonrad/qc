#!/usr/bin/env python

import glob
import sys
import re

def process_file(infile, outfile):
    """
    filter out count output from verg_and_count output
    """
    # look for "done."
    for line in infile:
        if line == "done\n":
            break
    for line in infile:
        if re.match("\d", line):
            outfile.write(line)

def collect_output(indir_name, outfile_name):
    with open(outfile_name, 'w') as outfile:
        for infile_name in glob.iglob('%s/*.out' % indir_name):
            with open(infile_name) as infile:
                process_file(infile, outfile)



def main():
    if len(sys.argv) < 3:
        print "usage: ./collect_output.py indir outfile"
        exit(-1)

    collect_output(sys.argv[1], sys.argv[2])


if __name__ == '__main__':
    main()
