#!/bin/sh
grep executing `grep -l 'signal 11' *.log` | cut -d'<' -f2 | cut -d':' -f1 | sort | uniq