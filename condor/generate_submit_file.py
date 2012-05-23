#!/usr/bin/env python

import sys
import os

alpha = 0.7 # k*dx

MAX_GAP = 0.05

def krange(k_low, k_high, delta_k):
    """
    yield pairs (k1, k2) partitioning the interval [k_low, k_high]
    into subintervals of width 2*delta_k
    """
    k = k_low
    while k < k_high:
        yield (k, k+2*delta_k)
        k += 2*delta_k

def kmetarange(k_low, k_high, delta_k, meta_width, meta_delta):
    """
    divide [k_low, k_high] into intervals of width meta_width*(k_low/k)**2
    (this width scaling yields a constant number of eigenvalues in each large interval)
    and spacing meta_delta
    divide these intervals into adjacent intervals of width 2*delta_k
    yield pairs (k1, k2) defining the endpoints of these small intervals

    """
    meta_k = k_low
    while meta_k < k_high:
        k = meta_k
        while k < min(meta_k + meta_width*(k_low/meta_k)**2, k_high):
            yield (k, k+2*delta_k)
            k += 2*delta_k
        meta_k += meta_delta

def find_gaps(delta, counts_file):
    """find k gaps in a file and yield k ranges of jobs to fill the gaps"""
    last_k = float('inf') # start last k at infinity to guarantee gap < MAX_GAP on first iteration
    with open(counts_file) as f:
        for line in f:
            k = float(line.split(',')[0])
            gap = k - last_k
            if gap > MAX_GAP:
                start = last_k
                while True:
                    end = min(start + delta, k)
                    yield (start, end)
                    start = end
                    if end == k:
                        break
            last_k = k

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
    delta_k = k - k1 # = k2 - k
    if billiard == 'qugrs':
        args = "-n qugrs_%f -l qugrs:1.0:0.4:0.7 -s oyooo:1.5:7:1 -x 5 -u -4 1 -k %f -V %f -a %f -M 9 -p 30" % (k, k, delta_k, alpha)
    elif billiard == 'qust':
        args = "-n qust_%f -l qust:2 -s vepwoo:1.2:40:1.5 -x 5 -u -k %f -V %f -a %f -M 9 -p 30" % (k, k, delta_k, alpha)
    else:
        raise ValueError('unknown billiard: %s' % billiard);
    submit_file.write("Arguments = %s\n" % args)
    submit_file.write("Output = %s_%f.out\n" % (billiard, k))
    submit_file.write("Error = %s_%f.err\n" % (billiard, k))
    submit_file.write("Log = %s_%f.log\n" % (billiard, k))
    submit_file.write("Queue\n\n")

def write_submit_file(k_low, k_high, delta, billiard):
    with open("verg_and_count.submit", 'w') as submit_file:
        write_prelude(submit_file)
        for (k1, k2) in krange(k_low, k_high, delta):
            write_job(submit_file, k1, k2, billiard)

def write_metarange_submit_file(k_low, k_high, delta, meta_width, meta_delta, billiard):
    with open("verg_and_count.submit", 'w') as submit_file:
        write_prelude(submit_file)
        for (k1, k2) in kmetarange(k_low, k_high, delta, meta_width, meta_delta):
            write_job(submit_file, k1, k2, billiard)

def write_filtered_submit_file(k_low, k_high, delta, billiard):
    with open("verg_and_count.submit", 'w') as submit_file:
        write_prelude(submit_file)
        for (k1, k2) in krange(k_low, k_high, delta):
            k = (k1+k2)/2
            if not os.path.exists('%s_%f.sta_bin' % (billiard, k)):
                write_job(submit_file, k1, k2, billiard)

def write_gaps_submit_file(delta, billiard, counts_file):
    with open("verg_and_count.submit", 'w') as submit_file:
        write_prelude(submit_file)
        for (k1, k2) in find_gaps(delta, counts_file):
            write_job(submit_file, k1, k2, billiard)


def usage():
    print "usage: ./generate_submit_file.py k_low k_high delta_0 billiard [f | g counts_file]"
    print "f filters on existing *.sta_bin files in current directory"
    print "g finds gaps in k in given counts_file"

def main():
    if len(sys.argv) < 5:
        usage()
        exit(-1)

    k_low = float(sys.argv[1])
    k_high = float(sys.argv[2])
    delta_k = float(sys.argv[3])
    billiard = sys.argv[4]

    if len(sys.argv) == 5:
        write_submit_file(k_low, k_high, delta_k, billiard)
    elif len(sys.argv) == 6 and sys.argv[5] == 'f':
        write_filtered_submit_file(k_low, k_high, delta_k, billiard)
    elif len(sys.argv) == 7 and sys.argv[5] == 'g':
        counts_file = sys.argv[6]
        write_gaps_submit_file(delta_k, billiard, counts_file)
    elif len(sys.argv) == 8 and sys.argv[5] == 'm':
        metarange_width = float(sys.argv[6])
        metarange_delta = float(sys.argv[7])
        write_metarange_submit_file(k_low, k_high, delta_k, metarange_width, metarange_delta, billiard)

    else:
        usage()
        exit(-1)

if __name__ == "__main__":
    main()
