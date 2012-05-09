#!/bin/sh

cp verg_no_main.cc verg_no_main.cc.backup
cp verg.cc verg_no_main.cc
sed -i -e 's/int main(/int verg_main(/' verg_no_main.cc
sed -i -e '/bas.set_type = NBASES;/ a\
  optind = 0;' verg_no_main.cc
sed -i -e 's/int verb.*/extern int verb;/' verg_no_main.cc
