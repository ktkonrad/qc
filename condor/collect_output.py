#!/usr/bin/env python

import glob
import sys

def process_file(infile, outfile):
    """
    filter out count output from verg_and_count output
    """
    # look for "done."
    for line in infile:
        if line == "done\n":
            break
    # output the rest of the lines
    for line in infile:
        outfile.write(line)

def collect_output(outfile_name):
    with open(outfile_name, 'w') as outfile:
        for infile_name in glob.iglob('*.out'):
            with open(infile_name) as infile:
                process_file(infile, outfile)



def main():
    if len(sys.argv) < 2:
        print "please supply an output file"
        exit(-1)

    collect_output(sys.argv[1])


if __name__ == '__main__':
    main()
