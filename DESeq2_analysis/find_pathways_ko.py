#!/usr/bin/env python
import sys, os

def write_to_file(dic, filename):
	f=open(filename, "w")
	for k in dic:
		f.write(str(dic[k])+"\t"+k+"\n")
	f.close()
		

paths={}
kos={}
for line in open ("img_annotation.master"):
	cut=line.strip().split("\t")
	if cut[4]=="NA": continue
	p=cut[4].split("|")
	ko=cut[1]
	paths[cut[0]]=p
	kos[cut[0]]=ko


over={}
under={}


for line in open ("significance_results.csv"):
	if "baseMean" in line: continue
	cut = line.strip().split(",")
	gene=cut[0][1:-1]
	if cut[6]=="NA": continue
	pval=float(cut[6])
	change=float(cut[2])
	if gene not in paths: continue
	ko=kos[gene]
	for path in paths[gene]:
		path = "\t".join(path.split(";"))
		path = path + "\t" + ko

		if pval<0.1:
			if change>0:
				if path not in over:
					over[path]=1
				else:
					over[path]+=1
			if change<0:
				if path not in under:
					under[path]=1
				else:
					under[path]+=1


write_to_file(over, "morning_genes.tab")
write_to_file(under, "evening_genes.tab")






