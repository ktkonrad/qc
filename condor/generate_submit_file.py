#!/usr/bin/env python

import sys
import os

alpha = 0.7 # k*dx

def krange(k_low, k_high, delta_k):
    """
    yield pairs (k1, k2) partitioning the interval [k_low, k_high]
    into subintervals of width 2*delta_k
    """
    k = k_low
    while k < k_high:
        yield (k, k+2*delta_k)
        k += delta_k

def write_prelude(submit_file):
    submit_file.write("""
Universe = vanilla
Requirements = OpSys == "LINUX" && Arch == "X86_64" && \
               Machine != "math-01.grid" && \
               Machine != "math-02.grid" && \
               Machine != "math-03.grid" && \
               Machine != "math-04.grid" && \
               Machine != "math-05.grid" && \
               Machine != "math-06.grid" && \
               Machine != "webwork.dartmouth.edu" && \
               Machine != "stanmore.dartmouth.edu"

should_transfer_files = YES
WhenToTransferOutput = ON_EXIT

transfer_input_files = verg, count

Executable = vc

""")

def write_job(submit_file, k1, k2, billiard):
    k = (k1+k2)/2
    dx = alpha / k
    delta_k = k - k1 # = k2 - k
    if billiard == 'qugrs':
        args = "-n qugrs_%f -l qugrs:1.0:0.4:0.7 -s oyooo:1.5:7:1 -u -4 1 -k %f -V %f -d %f -M 9 -p 30" % (k, k, delta_k, dx)
    elif billiard == 'qust':
        args = "-n qust_%f -l qust:2 -s vepwoo:1.2:40:1.5 -u -k %f -V %f -d %f -M 9 -p 30" % (k, k, delta_k, dx)
    else:
        raise ValueError('unknown billiard: %s' % billiard);
    submit_file.write("Arguments = %s\n" % args)
    submit_file.write("Output = %s_%f.out\n" % (billiard, k))
    submit_file.write("Error = %s_%f.err\n" % (billiard, k))
    submit_file.write("Log = %s_%f.log\n" % (billiard, k))
    submit_file.write("Queue\n\n")

def write_submit_file(k_low, k_high, delta_0, billiard):
    with open("verg_and_count.submit", 'w') as submit_file:
        write_prelude(submit_file)
        for (k1, k2) in krange(k_low, k_high, delta_0):
            write_job(submit_file, k1, k2, billiard)

def write_filtered_submit_file(k_low, k_high, delta_0, billiard):
    with open("verg_and_count.submit", 'w') as submit_file:
        write_prelude(submit_file)
        for (k1, k2) in krange(k_low, k_high, delta_0):
            k = (k1+k2)/2
            if not os.path.exists('run_%f.sta_bin' % k):
                write_job(submit_file, k1, k2, billiard)


def usage():
    print "usage: ./generate_submit_file.py k_low k_high delta_0 billiard [f]"


def main():
    if len(sys.argv) < 5:
        usage()
        exit(-1)
    elif len(sys.argv) > 6:
        usage()
        exit(-1)
    
    k_low = float(sys.argv[1])
    k_high = float(sys.argv[2])
    delta_k = float(sys.argv[3])
    billiard = sys.argv[4]
    
    if len(sys.argv) == 6 and sys.argv[5] == 'f':
        write_filtered_submit_file(k_low, k_high, delta_k, billiard)
    else:
        write_submit_file(k_low, k_high, delta_k, billiard)

if __name__ == "__main__":
    main()
