#!/usr/bin/env python

import sys

def krange(k_low, k_high, delta_low, delta_high):
    k = k_low + delta_low
    while k <= k_high - delta_high + 1e-6:
        yield k
        k += delta_low + delta_high

def write_submit_file(k_low, k_high, delta_low, delta_high):
    with open("verg_and_count.submit", 'w') as submit_file:
        submit_file.write("""
Universe = vanilla
Requirements = OpSys == "LINUX" && Arch == "X86_64"
should_transfer_files = YES
WhenToTransferOutput = ON_EXIT
notify_user = kyle.t.konrad@gmail.com

Executable = vc

""")

        for k in krange(k_low, k_high, delta_low, delta_high):
            args = "-n run_%f -l qugrs:1.0:0.4:0.7 -s oyooo:1.5:7:1 -u -4 1 -k %f -V %f:%f -d 0.001 -M 9 -p 20" % (k, k, delta_low, delta_high)
            submit_file.write("Arguments = %s\n" % args)
            submit_file.write("Output = run_%f.out\n" % k)
            submit_file.write("Error = run_%f.err\n" % k)
            submit_file.write("Log = run_%f.log\n" % k)
            submit_file.write("Queue\n\n")
    
def usage():
    print "usage: ./generate_submit_file.py k_low k_high delta_low delta_high"


def main():
    if len(sys.argv) < 5:
        usage()
        exit(-1)
    write_submit_file(*map(float, sys.argv[1:]))

if __name__ == "__main__":
    main()
