#!/usr/bin/env python
import sys

taxonomy={}
for line in open("contig_taxonomy.tab"):
	cut=line.strip().split("\t")
	if len(cut)<2:
		continue
	taxonomy[cut[0]]=cut[1]

brite2func={}
for line in open("brite2function.tab"):
        cut=line.strip().split("\t")
        if "Drug Dev" in cut[1] or "Human Dis" in cut[1] or "Organismal Sys" in cut[1]: continue
        brite2func[cut[0]]=cut[1]

ko2brite={}
for line in open("ko2brite.tab"):
	cut=line.strip().split("\t")
	brites=cut[1].split(";")
	good_brites=[]
	for i in range(len(brites)):
		brite=brites[i]
		if brite[2:] in brite2func:
			good_brites.append(brite)
	cut[1]="|".join(good_brites)
	if cut[1]!="":
		ko2brite[cut[0]]=cut[1]


print "\t".join(["#Locus", "KEGG ID", "Gene Product", "BRITE Pathways", "Functions", "Contig taxonomy"])

for line in open("img_annotation.products"):
	cut=line.strip().split("\t")
	cut[2]=cut[2].split(",")[0]
	if cut[2].startswith("KO:"):
		if cut[2].split(":")[1] in ko2brite:
			cut.append(ko2brite[cut[2].split(":")[1]])
			brites = ko2brite[cut[2].split(":")[1]].split("|")
			functions=[]
			for brite in brites:
				functions.append(brite2func[brite[2:]])
			cut.append("|".join(functions))
		else:
			cut.append("NA")
			cut.append("NA")
	else:
		cut.append(cut[2])
		cut.append("NA")
	ID=cut[2]; name=cut[1]
	cut[1]=ID; cut[2]=name
	contig = cut[0].split("-")[0]
	if contig in taxonomy:
		cut.append(taxonomy[contig])
	else:
		cut.append("NA")


	print "\t".join(cut)


