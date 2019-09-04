#!/usr/bin/env python
import sys

nano={}
for line in open("nanohaloarchaea_contigs.list"):
	nano[line.strip()]=None


for line in open("contig_taxonomy.tab"):
	cut=line.strip().split("\t")
	if cut[0] in nano:
		print cut[0] + "\t" + "Archaea;Nanoarchaeota;unclassified Nanoarchaeota"
		del nano[cut[0]]
	else:
		if len(cut)>1:
			print line.strip()

for contig in nano:
	print contig + "\t"+"Archaea;Nanoarchaeota;unclassified Nanoarchaeota"
