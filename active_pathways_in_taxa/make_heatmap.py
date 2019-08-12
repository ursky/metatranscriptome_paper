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
from PIL import Image
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from matplotlib.colors import LogNorm
from PIL import Image
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


def add_algae(taxonomy, filename):
	for line in open(filename):
		cut=line.strip().split("\t")
		if "Micromonas" in cut[1] or "Chlorophyta" in cut[1] or "Dolichomastix" in cut[1]:
			taxonomy[cut[0]]="Eukaryota;Chlorophyta;Mamiellophyceae;Mamiellales;Mamiellaceae;Dolichomastix"
	return taxonomy


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
		if cut[1] in to_collapse:
			cut[-1]="All"
		if cut[1]=="Oxidative":
			cut[1]="Oxidative phosphorylation"
		path=";".join(cut)

		if path not in df:
			df[path] = row
		else:
			df[path] += row
	return df.T


def get_pathway_ratios(filename, dna_df, rna_df, taxonomy, toi):
	print "making pathway table..."
	dna_data = dna_df.to_dict()
	rna_data = rna_df.to_dict()

	pathways={}
	gene_lists={}

	for line in open(filename):
		if line[0]=="#":
			continue
		cut = line.strip().split("\t")
		gene = cut[0]
		gene_id = cut[1]
		gene_name = cut[2]
		if "potassium uptake" in gene_name:
			cut[4] = "Metabolism;Energy metabolism;Opsins"
		contig = gene.split("-")[0]
		if cut[4]=="NA":
			continue
		paths = cut[4].split("|")
		for i,path in enumerate(paths):
			paths[i]=";".join(path.split(";")[0:3])
			if "Oxidative" in paths[i]:
				paths[i]+=" phosphorylation"

		# select for specific taxa (if used)
		if contig not in taxonomy:
			continue
		taxa = taxonomy[contig].split(";")
		if len(taxa)<2:
			continue
		taxon=taxa[1]
		if taxon not in toi:
			continue
		if taxon=="Euryarchaeota":
			taxon="Halobacteria"
		
		for path in paths:
			skip=False
			key_words=["Global", "Ribosome", "community", "Environmental"]
			for key in key_words:
				if key in path:
					skip=True
			if skip==True:
				continue 
			if path not in pathways:
				pathways[path]={}
				gene_lists[path]={}
			if taxon+"_DNA" not in pathways[path]:
				pathways[path][taxon+"_DNA"]=0
				pathways[path][taxon+"_RNA"]=0
				gene_lists[path][taxon]=[]
			
			

			if "Nanoa" in taxon:
				if "Transcription" in path or "Translation" in path:
					print gene_name
	




			
			for sample in dna_data:
				if contig not in dna_data[sample]:
					continue
				abund = dna_data[sample][contig]
				pathways[path][taxon+"_DNA"]+=abund
			for sample in rna_data:
				if gene not in rna_data[sample]:
					continue
				abund = rna_data[sample][gene] 
				pathways[path][taxon+"_RNA"]+=abund
			if gene_id not in gene_lists[path][taxon]:
				gene_lists[path][taxon].append(gene_id)

	pathways = drop_incomplete_paths(pathways, gene_lists)
	out_df = pd.DataFrame.from_dict(pathways)
	out_df = out_df.fillna(0).T
	return out_df


def get_ratios(df):
	out={}
	data = df.to_dict()
	for rna_sample in data:
		sample = rna_sample.split("_")[0]
		if "DNA" in rna_sample:
			continue
		dna_sample = sample+"_DNA"
		out[sample]={}
		for path in data[rna_sample]:
			rna = data[rna_sample][path]
			dna = data[dna_sample][path]
			if dna==0:
				ratio=0
			else:
				ratio=rna/dna
			out[sample][';'.join(path.split(';')[1:])] = ratio
	return pd.DataFrame.from_dict(out)
	

def drop_incomplete_paths(pathways, gene_lists):
	pathway_names={}
	for line in open("brite2function.tab"):
		cut=line.strip().split("\t")
		if "Oxidative" in cut[1]:
			cut[1]+=" phosphorylation"
		pathway_names[cut[0]] = cut[1]

	pathway_counts={}
	for line in open("ko2brite.tab"):
		cut=line.strip().split("\t")
		for path in cut[1].split(";"):
			if path[2:] not in pathway_names:
				continue
			path_name = pathway_names[path[2:]]
			if path_name not in pathway_counts:
				pathway_counts[path_name]=0
			pathway_counts[path_name]+=1					

	for path in pathways:
		if path not in pathway_counts:
			continue
		max_possible = pathway_counts[path]
		for taxon_sample in pathways[path]:
			taxon = taxon_sample.split("_")[0]
			n = len(gene_lists[path][taxon])
			if n<5 or n<max_possible*0.2:
				pathways[path][taxon_sample]=0
	return pathways



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


def get_top_pathways(df, n):
	# leave only top pathways
	rankings = []
	for index, row in df.iterrows():
	        rankings.append(max(row))
	rankings.sort(reverse=True)

	for index, row in df.iterrows():
	        if max(row)<rankings[n]:
	                df = df.drop([index])
	return df


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
taxonomy = add_algae(taxonomy, "contigs.mags")

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

if "recalc" in sys.argv:
	print "computing pathway table..."
	toi = ["Euryarchaeota", "Actinobacteria", "Bacteroidetes", "Chlorophyta", "Cyanobacteria", "Proteobacteria", "Nanoarchaeota"]
	df = get_pathway_ratios("img_annotation.master", dna_df, rna_df, taxonomy, toi)

	df.to_pickle("pathway_df.pkl")
	
else:
	df = pd.read_pickle("pathway_df.pkl")

##################   PROCESS HEATMAP MATRIX  ######################
print "processing matrix..."
df = collapse_paths(df)
# standardize by sample total
df = df.div(df.sum(axis=0), axis=1)
# compute ratios
df = get_ratios(df)


# standardize by sample ratio total
#df = 10*df.div(df.sum(axis=0), axis=1)

# get top pathways
df = get_top_pathways(df, 40)



##################   MAKE HEATMAP   ######################
print "making heatmap..."

g = sns.clustermap(df, figsize=(3,10), col_cluster=True, row_cluster=True, yticklabels=True, xticklabels=True, cmap="PiYG", method="weighted")
plt.savefig("figure.png", bbox_inches='tight', dpi=300)



