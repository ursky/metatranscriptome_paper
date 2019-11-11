#!/usr/bin/env python
import sys

interest={}
for line in open(sys.argv[1]):
	interest[line.strip()]=None

for line in open(sys.argv[2]):
	if line[0]==">":
		if line.strip().split()[0][1:] in interest:
			line=line.split()[0]
			p=True
		else:
			p=False
	if p==True and line.strip()!="":
		print line.strip()
