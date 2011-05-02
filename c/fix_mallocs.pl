#!/usr/bin/perl

while (<>) {
    s/(\w+ \*+)(.*)malloc/\1\2\(\1\)malloc/;
    print $_;
}
