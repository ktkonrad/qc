#!/usr/bin/env python

import sys

alpha = 0.5 # k*dx

def krange(k_low, k_high, delta):
    """
    yield pairs (k1, k2) partitioning the interval [k_low, k_high]
    into subintervals of width delta
    """
    k = k_low
    while k < k_high:
        yield (k, k+delta)
        k += delta

def write_submit_file(k_low, k_high, delta_0):
    with open("verg_and_count.submit", 'w') as submit_file:
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

        for (k1, k2) in krange(k_low, k_high, delta_0):
            k = (k1+k2)/2
            dx = alpha / k
            delta_k = k - k1 # = k2 - k
            args = "-n run_%f -l qugrs:1.0:0.4:0.7 -s oyooo:1.5:7:1 -u -4 1 -k %f -V %f -d %f -M 9 -p 20" % (k, k, delta_k, dx)
            submit_file.write("Arguments = %s\n" % args)
            submit_file.write("Output = run_%f.out\n" % k)
            submit_file.write("Error = run_%f.err\n" % k)
            submit_file.write("Log = run_%f.log\n" % k)
            submit_file.write("Queue\n\n")
    
def usage():
    print "usage: ./generate_submit_file.py k_low k_high delta_0"


def main():
    if len(sys.argv) != 4:
        usage()
        exit(-1)
    write_submit_file(*map(float, sys.argv[1:]))

if __name__ == "__main__":
    main()
