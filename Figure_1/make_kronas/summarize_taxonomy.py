#!/usr/bin/env python2
import sys


# load img taxonomy
taxonomy={}
for line in open("contig_taxonomy.tab"):
	cut=line.strip().split("\t")
	if len(cut)>1:
		taxonomy[cut[0]]="root;"+";".join(cut[1].split(";")[:6])

# load salmon abundances
summary={}
tot_abund=0
for line in open(sys.argv[1]):
	if not line.startswith("NODE"):
		continue
	cut=line.strip().split("\t")
	contig=cut[0]
	abund=float(cut[1]) * int(contig.split("_")[3])
	tot_abund+=abund
	if contig in taxonomy:
		tax=taxonomy[contig]
	else:
		tax="root"
	if tax not in summary:
		summary[tax]=abund
	else:
		summary[tax]+=abund


# collapse rare taxonomies
for i in reversed(range(10)):
	to_remove=[]
	for tax in summary:
		if len(tax.split(";"))!=i:
			continue
		abund = 0
		for tax2 in summary:
			if tax2.startswith(tax):
				abund+=summary[tax2]
		if 100.0*abund/tot_abund < 0.5:
			to_remove.append(tax)
	for tax in to_remove:
		new_tax = ";".join(tax.split(";")[:-1])
		if new_tax not in summary:
			summary[new_tax]=0
		summary[new_tax]+=summary[tax]
		del summary[tax]


for tax in summary:
	print str(summary[tax])+"\t"+"\t".join(tax.split(";"))

