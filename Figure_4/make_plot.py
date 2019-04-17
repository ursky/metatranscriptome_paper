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
		contig=cut[0]
		if "_length_" not in contig:
			continue
		length=int(contig.split("_")[3])
		if length<1000:
			continue
		abund=float(cut[n])
		data[contig]={}
		data[contig]["abundance"] = abund
		data[contig]["length"] = length
	return data


def load_tabular(filename, col=1):
	data={}
	for line in open(filename):
		cut=line.strip().split("\t")
		data[cut[0]]=cut[col]
	return data


def get_pathways(filename, data, taxonomy, euk=False):
	pathways={}
	for line in open(filename):
		if line[0]=="#":
			continue
		cut = line.strip().split("\t")
		if cut[3]=="NA":
			continue
		gene = cut[0]
		contig = gene.split("-")[0]
		if contig not in taxonomy:
			continue
		taxa=taxonomy[contig]
		if "Eukaryota" in taxa and euk==False:
			continue
		if "Eukaryota" not in taxa and euk==True:
			continue

		if contig in data:
			abund = data[contig]["abundance"]
		elif gene in data:
			abund = data[gene]["abundance"]
		else:
			continue

		paths=cut[3].split("|")
		path_names=cut[4].split("|")
		for i,path in enumerate(paths):
			name=path_names[i]
			full_path=name
			if full_path not in pathways:
				pathways[full_path]=0
			pathways[full_path]+=abund
	return pathways


def plot_compare(dna_paths, rna_paths, ax, euk=False):
	if euk==True:
		ratio_max=6
		ratio_min=0.3
		dna_min=20
	else:
		ratio_max=0.5
		ratio_min=0.032
		dna_min=100

	for pathway in dna_paths:
		if pathway not in rna_paths:
			continue
		label = "|".join(pathway.split(";")[1:])
		dna = dna_paths[pathway]
		rna = rna_paths[pathway]
		if rna==0:
			continue
		if dna==0:
			continue
		ratio = rna/dna

		ax.set_xscale('log')
		ax.set_yscale('log', basey=2)

		ax.scatter(dna, ratio, c="blue", s=20, alpha=0.3)
		ax.set_xlabel("Pathway abundance in metagenome", fontsize=label_font)
		ax.set_ylabel("Pathway expression relative to abundance (log2 fold)", fontsize=label_font)

		if ((ratio>ratio_max or ratio<ratio_min) and dna>dna_min):
				if label in pathway_rename:
					label=pathway_rename[label]
				label+=" > "
				print label, rna/dna
				ax.text(dna, ratio, label, fontsize=axis_font, alpha=0.8, horizontalalignment="right", verticalalignment='center')
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)







##################   SETUP PLOT   ######################
# main figure layout:
font = {'family': 'arial', 'weight': 'normal', 'size': 8}
plt.rc('font', **font)

sns.set_palette("colorblind")
fig = plt.figure(figsize=(8, 5), dpi=300)

axis_font=8
label_font=10
heading_font=16


##################   LOADING DATA   ######################
print "loading raw data..."
taxonomy = load_tabular("contig_taxonomy.tab")
dna_data = load_abundances("abundance_dna.tpm")
gene_data = load_abundances("abundance_genes.tpm")

pathway_rename = {
	"Signal transduction|FoxO": "FoxO signaling",
	"Energy metabolism|Photosynthesis": "Photosynthesis",
	"Translation|Ribosome": "Ribosome translation",
	"Carbohydrate metabolism|Glyoxylate": "Glyoxylate metabolism",
	"Energy metabolism|Carbon": "Carbon metabolism",
	"Glycan biosynthesis and metabolism|N-Glycan": "N-Glycan synthesis",
	"Xenobiotics biodegradation and metabolism|Steroid": "Steroid metabolism",
	"Biosynthesis of other secondary metabolites|Phenazine": "Phenazine synthesis",
	"Cellular community - prokaryotes|Quorum": "Quorum sensing",
	"Cellular community - prokaryotes|Biofilm": "Biofilm formation" }


##################   DRAW PROKAYOTE PLOT   ######################
print "plotting A..."
ax1 = fig.add_axes([0.1, 0.12, 0.37, 0.78])
dna_pathway_abundances = get_pathways("img_annotation.master", dna_data, taxonomy)
rna_pathway_abundances = get_pathways("img_annotation.master", gene_data, taxonomy)
plot_compare(dna_pathway_abundances, rna_pathway_abundances, ax1)


##################   DRAW EUKARYOTE PLOT   ######################
print "plotting D..."
ax2 = fig.add_axes([0.57, 0.12, 0.37, 0.78])
dna_pathway_abundances = get_pathways("img_annotation.master", dna_data, taxonomy, euk=True)
rna_pathway_abundances = get_pathways("img_annotation.master", gene_data, taxonomy, euk=True)
plot_compare(dna_pathway_abundances, rna_pathway_abundances, ax2, euk=True)


##################   FINISHING PLOT   ######################
ax1.annotate("A", xy=(-0.12, 1.05), xycoords="axes fraction", fontsize=heading_font)
ax2.annotate("B", xy=(-0.12, 1.05), xycoords="axes fraction", fontsize=heading_font)


#plt.subplots_adjust(left=0.1, right=0.95, top=0.9, bottom=0.1)
plt.savefig("figure_4.png", dpi=300)
#plt.savefig("figure_S2.eps", dpi=300)
#plt.show()



