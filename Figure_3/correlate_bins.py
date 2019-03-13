#!/usr/bin/env python
import sys, os
import matplotlib.pyplot as plt
import numpy as np
import random
import scipy.stats
import pandas as pd


def load_abundance_table(filename):
	df = pd.read_csv(filename, sep="\t", index_col=0, header=0)
	for col in df.columns:
		if col.startswith("ALL"):
			df = df.drop([col], axis=1)
	df = df.reindex_axis(sorted(df.columns), axis=1)
	data={}
	for mag in list(df.index):
		values = df.loc[mag].tolist()
		data[mag]=values
	return data
		


def plot_contigs(dna_dic, rna_dic, n, ax):
	print "plotting random contigs"
	#ax.set_xscale('log')
	#ax.set_yscale('log')
	ax.set_xlim(-0.05, 1.05)
	ax.set_ylim(-0.05, 1.05)
	for spine in ax.spines.values(): spine.set_visible(False)
	ax.grid(linestyle='--', alpha=0.5)
	contigs=dna_dic.keys()
	random.shuffle(contigs)
	for i,contig in enumerate(contigs):
		if i>n:
			break
		dna = dna_data[contig]
		rna = rna_data[contig]
		dna = standardize(dna)
		rna = standardize(rna)
		ax = draw_scatter(dna, rna, ax)
	


def plot_correlations(dna_dic, rna_dic, n, ax):
	print "calculating correlations"
	correlations=[]
	contigs=dna_dic.keys()
	random.shuffle(contigs)
	for i,contig in enumerate(contigs):
		if i>n:
			break
		dna = dna_data[contig]
		rna = rna_data[contig]
		corr, pval = scipy.stats.pearsonr(dna, rna)
		print contig+"\t"+str(corr)
		if corr<1 and corr>-1:
			correlations.append(corr)
	n, bins, patches = ax.hist(correlations, 10, facecolor='blue', alpha=0.5, rwidth=0.95)
	ax.set_xlim(-1,1)
	ax.set_ylim(-0.05, 7.5)
	ax.set_xlabel("Pearson's coefficient", fontsize=15)
	ax.set_ylabel("Number of MAGs", fontsize=15)
	for spine in ax.spines.values(): spine.set_visible(False)
	ax.grid(linestyle='--', alpha=0.5, zorder=-1)
	

def standardize(data):
	out=[]
	minimum = min(data)
	maximum = max(data)
	average = np.mean(data)
	if maximum==0:
		return data
	else:
		for val in data:
			out.append((val-minimum)/maximum)
			#out.append(val/average)
		return out


def draw_scatter(list1, list2, ax):
	ax.scatter(list1, list2, s=10, alpha=0.3, c='k')
	ax.set_xlabel("Total abundance", fontsize=15)
	ax.set_ylabel("Total activity", fontsize=15)
	return ax





# MAIN
font = {'family': 'arial', 'weight': 'normal', 'size': 12}
plt.rc('font', **font)

dna_data = load_abundance_table("DNA_QUANT/bin_abundance_table.tab")
rna_data = load_abundance_table("RNA_QUANT/bin_abundance_table.tab")


fig = plt.figure(figsize=(10,6))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

plot_contigs(dna_data, rna_data, 100, ax1)
plot_correlations(dna_data, rna_data, 100, ax2)


plt.tight_layout()
plt.savefig("bin_correlations.png", dpi=600)
plt.show()
