#!/bin/sh

cat qust_*_counts.txt | sort -n -u > temp
mv temp qust_all_counts.txt
cat qugrs_*_counts.txt | sort -n -u > temp
mv temp qugrs_all_counts.txt
rm temp