#!/usr/bin/env python
import sys

contigs={}
for line in open ("algae_genome.fa"):
	if line[0]==">":
		contigs[line[1:-1]]=None


for line in open(sys.argv[1]):
	line=line.strip()
	cut=line.split("\t")
	if cut[0]=="Name" or cut[0].split("-")[0] in contigs:
		print line






