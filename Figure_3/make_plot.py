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

def load_files(dirname, genes=False):
	print "Loading files from "+dirname
	data={}
	for filename in os.listdir(dirname):
		if not filename.startswith("9"):
			continue
		file_data = load_abundances(dirname+"/"+filename)
		if genes==True:
			file_data = combine_genes(file_data)
		for contig in file_data:
			abund = file_data[contig]["abundance"]
			if contig not in data:
				data[contig]={}
				data[contig]["length"]=file_data[contig]["length"]
				data[contig]["abundance"] = []
			data[contig]["abundance"].append(abund)
	return data


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


def get_mag_abundance(data, mags):
	mag_info = {}
	for contig in data:
		if contig not in mags:
			continue
		mag = mags[contig]
		if mag not in mag_info:
			mag_info[mag]={}
			mag_info[mag]["length"]=data[contig]["length"]
			mag_info[mag]["totals"]=data[contig]["abundance"]
		else:
			mag_info[mag]["length"]+=data[contig]["length"]
			mag_info[mag]["totals"]=map(sum, zip(data[contig]["abundance"], mag_info[mag]["totals"]))
	for mag in mag_info:
		mag_info[mag]["abundance"]=[]
		for abund in mag_info[mag]["totals"]:
			mag_info[mag]["abundance"].append(abund/mag_info[mag]["length"])
	return mag_info



def combine_genes(gene_data):
	print "calculating extression across contig based on all genes..."
	temp={}
	data={}
	for gene in gene_data:
		tpm=gene_data[gene]["abundance"]
		contig=gene.split("-")[0]
		if contig not in temp:
			temp[contig]=[]
			data[contig]={}
			data[contig]["length"]=gene_data[gene]["length"]
		temp[contig].append(tpm)
	for contig in data:
		mean=np.mean(temp[contig])
		data[contig]["abundance"]=mean
	return data


def standardize(data):
	out=[]
	minimum = min(data)
	maximum = max(data)
	if maximum==0:
		return data
	else:
		for val in data:
			out.append((val-minimum)/(maximum-minimum))
		return out



def plot_contigs(dna_data, rna_data, n, ax, alpha=0.2):
	print "plotting random contigs..."
	ax.set_xlim(0, 1)
	ax.set_ylim(0, 1)
	contigs=rna_data.keys()
	random.shuffle(contigs)
	ct=0
	for i,contig in enumerate(contigs):
		if ct>=n: break
		dna = dna_data[contig]["abundance"]
		rna = rna_data[contig]["abundance"]
		if max(rna)==0 or dna_data[contig]["length"]<10000:
			continue
		dna = standardize(dna)
		rna = standardize(rna)
		length = dna_data[contig]["length"]
		ax.scatter(dna, rna, s=20, alpha=alpha, c='blue', edgecolors='none')
		ct+=1
	ax.set_xlabel("Contig DNA coverage (CPM)", fontsize=label_font)
	ax.set_ylabel("Mean expression or genes on contigs (TPM)", fontsize=label_font)
	#ax.plot([0.00001, 1000000], [0.000001, 10000000], ls="--", c='k', alpha=0.5)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)

	
def plot_correlations(dna_data, rna_data, n, ax):
	print "calculating correlations..."
	correlations=[]
	contigs=rna_data.keys()
	random.shuffle(contigs)
	ct=0
	for i,contig in enumerate(contigs):
		if ct>=n:
			break
		dna = dna_data[contig]["abundance"]
		rna = rna_data[contig]["abundance"]
		if max(rna)==0  or dna_data[contig]["length"]<10000:
			continue
		corr, pval = scipy.stats.pearsonr(dna, rna)
		correlations.append(corr)
		if "T17" in contig:
			print corr, contig
		ct+=1

	n, bins, patches = ax.hist(correlations, 20, facecolor='blue', alpha=0.5, rwidth=0.95)
	ax.set_xlim(-1,1)
	ax.set_xlabel("Pearson's Coefficient", fontsize=label_font)
	ax.set_ylabel("Number of contigs", fontsize=label_font)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)




##################   SETUP PLOT   ######################
# main figure layout:
font = {'family': 'arial', 'weight': 'normal', 'size': 8}
plt.rc('font', **font)

sns.set_palette("colorblind")
fig = plt.figure(figsize=(6, 6), dpi=300)

axis_font=8
label_font=10
heading_font=16


##################   LOADING DATA   ######################
if "reload" in sys.argv:
	reload=True
else:
	reload=False
if reload==True:
	print "loading raw data..."
	dna_data = load_files("abundances_dna")
	rna_data = load_files("abundances_genes", genes=True)
	with open('dna_data.json', 'w') as fp:
		json.dump(dna_data, fp)
	with open('rna_data.json', 'w') as fp:
		json.dump(rna_data, fp)
else:
	print "loading processed data..."
	with open('dna_data.json', 'r') as fp:
		dna_data = json.load(fp)
	with open('rna_data.json', 'r') as fp:
		rna_data = json.load(fp)

print "estimating MAG coverage..."
mags = load_tabular("contigs.mags")
mag_dna_data = get_mag_abundance(dna_data, mags)
mag_rna_data = get_mag_abundance(rna_data, mags)


##################   DRAW CONTIG PLOTS   ######################
print "plotting A..."
ax1 = fig.add_axes([0.1, 0.57, 0.37, 0.37])
plot_contigs(dna_data, rna_data, 1000, ax1, alpha=0.05)

print "plotting B..."
ax2 = fig.add_axes([0.57, 0.57, 0.37, 0.37])
plot_correlations(dna_data, rna_data, 1000, ax2)


##################   DRAW MAG PLOTS   ######################
print "plotting C..."
ax3 = fig.add_axes([0.1, 0.07, 0.37, 0.37])
plot_contigs(mag_dna_data, mag_rna_data, 1000, ax3)
ax3.set_ylabel("Mean expression or genes in MAGs (TPM)", fontsize=label_font)

print "plotting D..."
ax4 = fig.add_axes([0.57, 0.07, 0.37, 0.37])
plot_correlations(mag_dna_data, mag_rna_data, 1000, ax4)
ax4.set_ylabel("Number of MAGs", fontsize=label_font)


##################   FINISHING PLOT   ######################
ax1.annotate("A", xy=(-0.12, 1.05), xycoords="axes fraction", fontsize=heading_font)
ax2.annotate("B", xy=(-0.12, 1.05), xycoords="axes fraction", fontsize=heading_font)
ax3.annotate("C", xy=(-0.12, 1.05), xycoords="axes fraction", fontsize=heading_font)
ax4.annotate("D", xy=(-0.12, 1.05), xycoords="axes fraction", fontsize=heading_font)


#plt.subplots_adjust(left=0.1, right=0.95, top=0.9, bottom=0.1)
plt.savefig("figure_3.png", dpi=300)
#plt.savefig("figure_S2.eps", dpi=300)
#plt.show()



