#!/usr/bin/env python
import sys

foi={}
for line in open(sys.argv[1]):
	foi[line.strip()]=None


for line in open(sys.argv[2]):
	if line[0]==">":
		contig=line[1:].split()[0]
		if contig in foi:
			p=True
		else:
			p=False
	if p==True:
		print line.strip()
