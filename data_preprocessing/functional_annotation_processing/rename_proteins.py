#!/usr/bin/env python
# ./rename_products.py 203202.assembled.names_map 203202.assembled.product_names 
import sys

names={}
for line in open(sys.argv[1]):
	cut=line.strip().split("\t")
	names[cut[1]]=cut[0]

for line in open(sys.argv[2]):
	line = line.strip()
	if line[0]==">":
		name=line[1:]
		contig=name[:17]
		id=name[17:]
		contig = names[contig]
		name=contig+"-"+id
		line = ">" + contig +"|"+name
		l = int(contig.split("_")[3])
	if l>5000:
		print line
