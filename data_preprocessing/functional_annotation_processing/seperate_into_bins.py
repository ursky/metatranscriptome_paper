#!/usr/bin/env python
import sys, os

membership={}
file_handles={}
bindir="/scratch/gu/RNA_2017_PROJECT/METAGENOME/BINNING/metawrap_bins"
for f in os.listdir(bindir):
	binname="_".join(f.split(".")[:-1])
	file_handles[binname] = open(binname+".gff", "w")
	for line in open(bindir+"/"+f):
		if line[0]==">":
			membership[line[1:-1]]=binname



for line in open("img_annotation.gff"):
	line=line.strip()
	cut=line.split("\t")
	if cut[0] not in membership:
		continue
	binname=membership[cut[0]]
	file_handles[binname].write(line+"\n")


for f in file_handles:
	file_handles[f].close()
	


