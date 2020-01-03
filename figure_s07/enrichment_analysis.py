#!/usr/bin/env python
import sys
import random

level=3

paths=[]
for line in open("all_genes.tab"):
	cut=line.strip().split("\t")
	ct=int(cut[0])
	path="\t".join(cut[1:level])
	for i in range(ct):
		paths.append(path)


total=0
categories={}
for line in open (sys.argv[1]):
	cut=line.strip().split("\t")
	ct=int(cut[0])
	total+=ct
	path="\t".join(cut[1:level])
	if path not in categories:
		categories[path]=ct
	else:
		categories[path]+=ct


for path in categories:
	ct=categories[path]
	positive=0; negative=0
	for i in range(1000):
		test_ct=0
		for j in range(total):
			if random.choice(paths)==path:
				test_ct+=1
		if test_ct>=ct:
			positive+=1
		else:
			negative+=1
		#print positive, negative
	pval=1.0*positive/(positive+negative)
	if pval<0.05:
		print pval, path		

