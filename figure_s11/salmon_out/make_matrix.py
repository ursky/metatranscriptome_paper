#!/usr/bin/env python
import sys

data={}
for d in sys.argv[1:]:
	sample="_".join(d.split("-")[:2])
	for line in open(d+"/quant.sf"):
		if "Name" in line:
			continue
		cut=line.strip().split("\t")
		abund=cut[3]
		gene=cut[0]
		if gene not in data:
			data[gene]={}
		data[gene][sample]=abund

cut=["Gene"]
for sample in data["NODE_1_length_176649_cov_9.362153-1"]:
	cut.append(sample)
print "\t".join(cut)


for gene in data:
	cut=[gene]
	for sample in data[gene]:
		cut.append(data[gene][sample])
	print "\t".join(cut)
