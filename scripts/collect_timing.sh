#!/bin/sh

if [ $# -lt 2 ]
then
    echo 'usage: ./collet_timing.sh indir outfile'
    exit
fi

grep -h -P 'found|took' $1*.out > $2
./parse_timing_data.py $2