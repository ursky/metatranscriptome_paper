#!/usr/bin/env python
import sys

contigs={}
# load contigs in bin
for line in open (sys.argv[1]):
	if line[0]==">":
		contigs[line[1:-1]]=None


#curate salmon quant file
for line in open(sys.argv[2]):
	line=line.strip()
	cut=line.split("\t")
	if cut[0]=="Name" or cut[0].split("-")[0] in contigs:
		print line






