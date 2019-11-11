#!/usr/bin/env python
# ./rename_annotation.py 203202.assembled.names_map 203202.assembled.gff
import sys

names={}
for line in open(sys.argv[1]):
	cut=line.strip().split("\t")
	names[cut[1]]=cut[0]

for line in open(sys.argv[2]):
	cut=line.strip().split("\t")
	contig=cut[0]
	info=cut[8].split(";")

	for i in range(len(info)):
		f=info[i]
		c=f.split("=")
		if c[0]=="locus_tag":
			id=c[1][17:]
			c[1]=names[contig]+'-'+id
			info[i]=c[0]+"="+c[1]

	cut[8]=";".join(info)

	cut[0]=names[contig]
	#cut[0]="_".join(names[contig].split("_")[:4])

	if cut[6]=="1": cut[6]="+"
	if cut[6]=="-1": cut[6]="-"
	if int(cut[0].split("_")[3])<1000: continue

	print "\t".join(cut)
