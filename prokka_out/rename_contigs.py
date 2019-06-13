#!/usr/bin/env python2
import sys

for line in open(sys.argv[1]):
	if line[0]==">":
		print "-".join(line.strip().split("_")[:2])
	else:
		print line.strip()


