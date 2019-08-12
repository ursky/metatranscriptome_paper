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
from scipy.stats import spearmanr
from matplotlib.colors import LogNorm
from PIL import Image
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
		if cut[4]=="NA":
			continue
		paths = cut[4].split("|")
		for i,path in enumerate(paths):
			paths[i]=";".join(path.split(";")[1:3])

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


def split_dna_rna(df):
	dna={}; rna={}
	data = df.to_dict()
	for k in data:
		if "DNA" in k:
			dna[k]=data[k]
		if  "RNA" in k:
			rna[k]=data[k]
		
	dna_df = pd.DataFrame.from_dict(dna)
	rna_df = pd.DataFrame.from_dict(rna)
	return dna_df, rna_df


def compute_variance(df):
	variances=[]
	paths=[]
	for k in df:
		raw_values = list(df[k])
		values=[]
		for i in raw_values:
			values.append(i/max(raw_values))
		variance= np.var(values)
		variances.append(variance)
		paths.append(k)
	# sort the pathway names by their variances
	sorted_paths = [x for _,x in sorted(zip(variances,paths), reverse=True)]
	return sorted_paths, variances


def standardize_to_mean(df):
	means=[]
	for k in df: 
		means.append(np.mean(df[k]))
	lists=[]
	for index, row in df.iterrows():
		row = row.div(means)
		lists.append(list(row))
	return lists



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
if "recalc" in sys.argv:
	df = get_pathways("img_annotation.master", dna_df, rna_df, taxonomy, taxon="All")
	df.to_pickle("pathway_df.pkl")
else:
	df = pd.read_pickle("pathway_df.pkl")

##################   PROCESS MATRIX  ######################
print "processing matrix..."
for index, row in df.iterrows():
	if sum(row)<100:
		df = df.drop([index])

# standardize by sample total coverage
df = df.div(df.sum(axis=0), axis=1)

# get seperate matrix for pathway abundance in DNA and RNA
dna_df, rna_df = split_dna_rna(df)

# standardize by maximum value in each row
#dna_df = dna_df.div(dna_df.max(axis=1), axis=0).T
#rna_df = rna_df.div(rna_df.max(axis=1), axis=0).T
dna_df=dna_df.T
rna_df=rna_df.T

# compute variances
dna_order, dna_var = compute_variance(dna_df)
rna_order, rna_var = compute_variance(rna_df)

# reorder matrix columns by decending variance
sort_by=dna_order
dna_df = dna_df.reindex(sort_by, axis='columns')
rna_df = rna_df.reindex(sort_by, axis='columns')
dna_var = [x for _,x in sorted(zip(sort_by,dna_var))]
rna_var = [x for _,x in sorted(zip(sort_by,rna_var))]


##################   MAKE PLOTS   ######################
fig = plt.figure(figsize=(6, 6.5), dpi=300)
ax1 = fig.add_axes([0.08, 0.68, 0.9, 0.25])
ax2 = fig.add_axes([0.08, 0.58, 0.9, 0.1])
ax3 = fig.add_axes([0.08, 0.2, 0.9, 0.25])
ax4 = fig.add_axes([0.08, 0.1, 0.9, 0.1])

print "making heatmaps..."
sns.heatmap(dna_df.div(dna_df.max(axis=0)), yticklabels=False, xticklabels=False, cmap="magma", cbar=False, ax=ax1)
sns.heatmap(rna_df.div(rna_df.max(axis=0)), yticklabels=False, xticklabels=False, cmap="magma", cbar=False, ax=ax3)

print "making total plots..."
for series in standardize_to_mean(dna_df):
	ax2.plot(series, c='gray', alpha=0.5, linewidth=1)
for series in standardize_to_mean(rna_df):
	ax4.plot(series, c='gray', alpha=0.5, linewidth=1)


# clean up plots
for ax in [ax2, ax4]:
	ax.set_xlim(0,len(dna_var))
	ax.set_ylim(0.5,1.5)
	ax.set_xticklabels([])
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)

# add labels
ax = fig.add_axes([0,0,1,1])
ax.text(0.52, 0.95, "Pathway Abundance in Metagenome", horizontalalignment='center', fontsize=20)
ax.text(0.52, 0.47, "Pathway Expression in Metatranscriptome", horizontalalignment='center', fontsize=20)
ax.text(0.02, 0.95, "A", horizontalalignment='center', fontsize=24)
ax.text(0.02, 0.47, "B", horizontalalignment='center', fontsize=24)
ax.text(0.52, 0.06, "Pathways", horizontalalignment='center', fontsize=16)
ax.text(0.52, 0.54, "Pathways", horizontalalignment='center', fontsize=16)
ax.text(0.02, 0.32, "Samples", verticalalignment='center', fontsize=16, rotation=90)
ax.text(0.02, 0.8, "Samples", verticalalignment='center', fontsize=16, rotation=90)
#ax.text(0.02, 0.15, "Std(TPM)", verticalalignment='center', fontsize=16, rotation=90)
#ax.text(0.02, 0.63, "Std(CPM)", verticalalignment='center', fontsize=16, rotation=90)

plt.axis("off")



#plt.savefig("figure.png", bbox_inches='tight', dpi=300)
plt.savefig("figure.png", dpi=300)





