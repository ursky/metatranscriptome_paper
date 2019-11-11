#!/usr/bin/env python
print "loading python packages..."
import sys, getopt, os
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
from matplotlib.cbook import get_sample_data
import seaborn as sns
sns.set_color_codes()
import operator as op
import numpy as np
import math
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from matplotlib.colors import LogNorm
import scipy.stats
import json
import random


def load_abundances(file_name):
	print "loading abundances from "+file_name+"..."
	data={}
	for line in open(file_name):
		if line.startswith("transcript"):
			n=1
			continue
		if line.startswith("Name"):
			n=3
			continue
		cut =line.strip().split("\t")
		if "_length_" not in cut[0]:
			continue

		abund=float(cut[n])
		if "-" in cut[0]:
			gene = cut[0]
			data[gene] = abund
		else:
			contig = cut[0]
			data[contig] = abund
	return data


def load_tabular(filename, col=1):
	data={}
	for line in open(filename):
		cut=line.strip().split("\t")
		data[cut[0]]=cut[col]
	return data


def collapse_paths(pathway_df):
	print "colapsing redundant pathways..."
	df = pd.DataFrame()
	to_collapse=["Xenobiotics biodegradation and metabolism",
		"Metabolism of terpenoids and polyketides",
		"Metabolism of cofactors and vitamins",
		"Glycan biosynthesis and metabolism",
		"Biosynthesis of other secondary metabolites",
		"Replication and repair",
		"Amino acid metabolism",
		"Signal transduction",
		"Metabolism of other amino acids",
		"Transport and catabolism",
		"Cell growth and death",
		"Cellular community - prokaryotes",
		"Nucleotide metabolism"]
	for path,row in pathway_df.iterrows():
		cut = path.split(";")
		if cut[0] in to_collapse:
			cut[1]="All"
		if cut[1]=="Oxidative":
			cut[1]="Oxidative phosphorylation"
		path=";".join(cut)

		if path not in df:
			df[path] = row
		else:
			df[path] += row
	return df.T


def drop_incomplete_paths(pathways, pathway_counts):
	count_list=[]
	for path in pathway_counts:
		count = pathway_counts[path]
		count_list.append(count)
	count_list.sort()
	min_count = max(10, count_list[len(count_list)/2])
	print "min gene count per pathway = "+str(min_count)
	for path in pathway_counts:
		count = pathway_counts[path]
		if count<min_count and "Opsin" not in path:
			for sample in pathways:
				pathways[sample].pop(path)
	return pathways


def get_pathways(filename, dna_df, rna_df, taxonomy, taxon="All"):
	print "making pathway table..."
	dna_data = dna_df.to_dict()
	rna_data = rna_df.to_dict()

	pathways={}
	pathway_counts={}

	for sample in dna_data:
		pathways[sample]={}
	for sample in rna_data:
		pathways[sample]={}

	for line in open(filename):
		if line[0]=="#":
			continue
		cut = line.strip().split("\t")
		if cut[3]=="NA":
			continue
		gene = cut[0]
		gene_name = cut[2]
		contig = gene.split("-")[0]
		if "opsin" in gene_name:
			cut[4]="Metabolism;Energy metabolism;Opsins"
		if cut[4]=="NA":
			continue
		paths = cut[4].split("|")
		for i,path in enumerate(paths):
			paths[i]=";".join(path.split(";")[1:3])

		# select for specific taxa (if used)
		if taxon!="All":
			if contig not in taxonomy:
				continue
			taxa=taxonomy[contig]
			if taxon not in taxa:
				continue
		for path in paths:
			if "Global" in path:
				continue 
			if path not in pathway_counts:
				pathway_counts[path]=1
			else:
				pathway_counts[path]+=1

			for sample in dna_data:
				abund = dna_data[sample][contig]
				if path not in pathways[sample]:
					pathways[sample][path]=abund
				else:
					pathways[sample][path]+=abund
			for sample in rna_data:
				if gene not in rna_data[sample]:
					continue
				abund = rna_data[sample][gene] 
				if path not in pathways[sample]:
					pathways[sample][path]=abund
				else:
					pathways[sample][path]+=abund
	
	pathways = drop_incomplete_paths(pathways, pathway_counts)
	pathways_df = pd.DataFrame.from_dict(pathways)
	return pathways_df


def compute_variance(df, taxon="All"):
	data=df.T.to_dict()
	vars_df = pd.DataFrame(columns=['Taxon', 'Type', 'Variance'])
	for path in data:
		for type in ["DNA", "RNA"]:
			values=[]
			for sample in data[path]:
				if type in sample:
					values.append(data[path][sample])
			values = standardize(values)
			var = np.var(values)
			line = pd.DataFrame({"Taxon":[taxon], "Type":[type], "Variance":var})
			vars_df = vars_df.append(line, ignore_index = True) 
	sns.violinplot(x="Taxon", y="Variance", hue="Type", data=vars_df, split=True, palette={"DNA": "red", "RNA": "blue"})
	plt.show()


def standardize(a_list):
	max_element = max(a_list)
	for i,element in enumerate(a_list):
		a_list[i] = element/max_element
	return a_list


##################   SETUP PLOT   ######################
# main figure layout:
font = {'family': 'arial', 'weight': 'normal', 'size': 8}
plt.rc('font', **font)
sns.set_palette("colorblind")
axis_font=8
label_font=10
heading_font=16


##################   LOADING DATA   ######################
taxonomy = load_tabular("contig_taxonomy.tab")

if "reload" in sys.argv:
	print "loading raw data..."
	dna_data={}
	rna_data={}
	for replicate in range(1,7):
		for time in ["9AM", "9PM"]:
			sample = time+"-"+str(replicate)
			dna_data["DNA-"+sample] = load_abundances(sample+"-DNA.quant")
			rna_data["RNA-"+sample] = load_abundances(sample+"-mRNA.quant")

	dna_df = pd.DataFrame.from_dict(dna_data)
	rna_df = pd.DataFrame.from_dict(rna_data)
	dna_df.to_pickle("df_dna.pkl")
	rna_df.to_pickle("df_rna.pkl")
else:
	dna_df = pd.read_pickle("df_dna.pkl")
	rna_df = pd.read_pickle("df_rna.pkl")


##################   MAKE FUNCTIONAL TABLE   ######################
print "computing pathway table..."
if len(sys.argv)>1:
	pathway_df = get_pathways("img_annotation.master", dna_df, rna_df, taxonomy, taxon=sys.argv[1])
	df = collapse_paths(pathway_df)
	pathway_df.to_pickle("pathway_df.pkl")
else:
	pathway_df = get_pathways("img_annotation.master", dna_df, rna_df, taxonomy, taxon="All")
	pathway_df = pd.read_pickle("pathway_df.pkl")
	df = collapse_paths(pathway_df)


##################   PROCESS HEATMAP MATRIX  ######################
print "processing matrix..."
for index, row in df.iterrows():
	if sum(row)<100:
		df = df.drop([index])

# standardize by sample total coverage
sums = []
for col in df.columns.values:
        sums.append( sum(df[col]) )
df = df.div(df.sum(axis=0), axis=1)
df*=np.mean(sums)



# standardize by maximum value in each row
df = df.div(df.max(axis=1), axis=0)


# log standardize
#df += 1
#df = df.apply(np.log)


# sort pathways by expression
dna_col_list=[]; rna_col_list=[]
for col in list(df):
	if "RNA" in col:
		rna_col_list.append(col)
	else:
		dna_col_list.append(col)
df['DNA'] = df[dna_col_list].mean(axis=1)
df['RNA'] = df[rna_col_list].mean(axis=1)
df['DIF'] = 2*df["DNA"]-df["RNA"]
df = df.sort_values(by=["DIF"])
df = df.drop(["DNA", "RNA", "DIF"], axis=1)





##################   MAKE HEATMAP   ######################
print "making heatmap..."

lut=[]
for sample in df.columns.values:
        #if "DNA" in sample:
        #       df = df.drop([sample], axis=1)
        #       continue
        if "DNA" in sample:
                c="blue"
        elif "RNA" in sample:
                c="red"
        else:
                c="w"
        lut.append(c)

print df
print df.shape
g = sns.clustermap(df, figsize=(5,8), col_colors=lut, col_cluster=True, row_cluster=True, yticklabels=True, xticklabels=True, cmap="magma", method="weighted")


#plt.savefig("figure_4.png", bbox_inches='tight', dpi=300)
if len(sys.argv)>1:
	g.fig.suptitle(sys.argv[1], fontsize=40)
	plt.savefig(sys.argv[1]+".png", bbox_inches='tight', dpi=300)
plt.savefig("figure_std.png", bbox_inches='tight', dpi=300)





