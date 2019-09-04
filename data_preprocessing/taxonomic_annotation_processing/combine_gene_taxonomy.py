#!/usr/bin/env python
import sys

def add_to_tree( tree, tax_list, length ):
	if len(tax_list)==0: return tree
	else:
		if tax_list[0] not in tree: 
			tree[tax_list[0]]=[length, {}]
		else: 
			tree[tax_list[0]][0]+=length
		add_to_tree( tree[tax_list[0]][1], tax_list[1:], length )
	return tree


def traverse(tree, taxonomy, weight):
	if len(tree)==0:
		return taxonomy
	total_score=0
	max_score=0
	max_class=""
	for k in tree:
		total_score+=tree[k][0]
		if tree[k][0]>max_score:
			max_score=tree[k][0]
			max_class=k
	if weight!=0: total_score=weight

	if 100*max_score/total_score>50:
		taxonomy.append(max_class) # +"_"+str(100*max_score/total_score))
		taxonomy=traverse(tree[max_class][1], taxonomy, tree[max_class][0])
	else:
		return taxonomy
	return taxonomy


tax_tree={}; contig=""
for line in open("img_annotation.taxonomy"):
	cut=line.strip().split("\t")

	if cut[0].split("-")[0]!=contig and tax_tree!={}:
		if "CHLO" in contig or "MITO" in contig:
			consensus=["Eukaryota","Chlorophyta","Mamiellophyceae","Mamiellales","Mamiellaceae","Micromonas"]
		else:
			consensus=traverse(tax_tree, [], 0)
		print contig + "\t" + ";".join(consensus)
		tax_tree={}

	contig=cut[0].split("-")[0]
	taxcut=cut[4].split(";")
	tax_tree = add_to_tree(tax_tree, taxcut, 1)


consensus=traverse(tax_tree, [], 0)
print contig + "\t" + ";".join(consensus)

