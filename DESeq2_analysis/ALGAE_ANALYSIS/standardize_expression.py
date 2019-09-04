#!/usr/bin/env python
import sys

bins={}
#for line in open("/Users/guritsk1/Desktop/RNA_SEQ_2017/quant_dna/bin_abundance_table.tab"):
for line in open("/Users/guritsk1/Desktop/RNA_SEQ_2017/quant_rna/bin_abundance_table.tab"):
	cut = line.strip().split("\t")
	if cut[0]=="Genomic bins":
		for i in range(1, len(cut)):
			cut[i] = "-".join(cut[i].split("-")[:2])
			if "ALL" in cut[i]:
				cut[i]="ALL"
		head = cut[:]
	else:
		bins[cut[0]]={}
		for i in range(1, len(cut)):
			bins[cut[0]][head[i]] = float(cut[i])
		
for line in open(sys.argv[1]):
	sample = "-".join(sys.argv[1].split("/")[-2].split("-")[:2])
	cut = line.strip().split("\t")
	if cut[0]=="Name":
		print line.strip()
	else:
		tpm=float(cut[4])
		tpm = tpm * bins["bin.algae"]["ALL"] / bins["bin.algae"][sample]
		cut[4]=str(tpm)
		print "\t".join(cut)
