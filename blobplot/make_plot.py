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
		if length<3000:
			continue
		abund=float(cut[n])
		data[contig]={}
		data[contig]["abundance"] = abund
		data[contig]["length"] = length
	return data


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


def make_taxonomy_colors(taxonomy, data, depth=4):
	for contig in data:
		if contig not in taxonomy:
			taxon="NA"
			c="white"
		else:
			taxa = taxonomy[contig]
			if len(taxa)<depth+1:
				taxon="white"
				c='white'
			else:
				taxon = taxa[depth]
				if taxon not in scheme:
					c='white'
				else:
					c = scheme[taxon]
		data[contig]["taxonomy"]=taxon
		data[contig]["color"]=c
	return data


def load_taxonomy(file_name):
	data={}
	for line in open(file_name):
		cut=line.strip().split("\t")
		if len(cut)<2:
			continue
		if "Eukaryota" in cut[1]:
			cut[1]="Eukaryota;Chlorophyta;Mamiellophyceae;Mamiellales"
		data[cut[0]]=cut[1].split(";")
	return data


def load_tabular(filename, col=1):
	data={}
	for line in open(filename):
		cut=line.strip().split("\t")
		data[cut[0]]=cut[col]
	return data
		

def representation_by_taxonomy(dna_data, gene_data):
	print "calculating RNA representation in taxonomies..."
	data_out={}
	for gene in gene_data:
		contig = gene.split("-")[0]
		taxon = dna_data[contig]["taxonomy"]
		if taxon=="NA":
			continue
		dna=dna_data[contig]["abundance"]
		rna=gene_data[gene]["abundance"]
		
		if dna==0:
			continue
		ratio=rna/dna
		if taxon not in data_out:
			data_out[taxon]={}
			data_out[taxon]["ratios"]=[]
			data_out[taxon]["abundances"]=[]
			data_out[taxon]["totals"]=0
			data_out[taxon]["taxon"]=taxon
		data_out[taxon]["ratios"].append(ratio)
		data_out[taxon]["abundances"].append(dna)
		data_out[taxon]["totals"]+=rna
	return data_out


def representation_by_mag(dna_data, gene_data, mags, mag_taxonomy):
	print "calculating RNA representation in MAGs..."
	mag_info={}
	for gene in gene_data:
		contig = gene.split("-")[0]
		if contig not in mags:
			continue
		mag=mags[contig]
		taxon = mag_taxonomy[mag].split(";")[1]
		if mag not in mag_info:
			mag_info[mag]={}
			mag_info[mag]["abundances"]=[]
			mag_info[mag]["ratios"]=[]
			mag_info[mag]["taxon"]=taxon
			mag_info[mag]["totals"]=0
		dna=dna_data[contig]["abundance"]
		rna=gene_data[gene]["abundance"]
		ratio = rna/dna

		mag_info[mag]["abundances"].append(dna)
		mag_info[mag]["ratios"].append(ratio)
		mag_info[mag]["totals"]+=rna
	return mag_info


def draw_scatter(dna_data, rna_data, ax, archaea=False):
	x=[]; y=[]; sizes=[]; colors=[]
	for contig in dna_data:
		if archaea==False and dna_data[contig]["color"]==scheme["Euryarchaeota"]:
			continue
		if archaea==True and dna_data[contig]["color"]!=scheme["Euryarchaeota"]:
			continue
		x.append(dna_data[contig]["abundance"])
		if contig in rna_data:
			y.append(rna_data[contig]["abundance"])
		else:
			y.append(0)
		sizes.append(dna_data[contig]["length"]/500)
		colors.append(dna_data[contig]["color"])
	ax.scatter(x, y, s=sizes, alpha=0.1, c=colors)
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlim(1, max(x))
	ax.set_ylim(0.01, max(y))
	ax.set_xlabel("Contig DNA coverage (CPM)", fontsize=label_font)
	ax.set_ylabel("Mean expression of genes on contig (TPM)", fontsize=label_font)
	ax.plot([0.00001, 1000000], [0.0000001, 1000000], ls="--", c='k', alpha=0.3)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	plt.tick_params(labelsize=axis_font)
	plt.grid(alpha=0.3, linestyle='--')
	


def draw_averages (taxon_info, ax):
	xs=[];ys=[]
	for taxa in taxon_info:
		total_abund = taxon_info[taxa]["totals"]
		taxon = taxon_info[taxa]["taxon"]
		if taxon in scheme:
			color=scheme[taxon]
		else:
			continue
		abund = np.mean(taxon_info[taxa]["abundances"])
		ratio = np.mean(taxon_info[taxa]["ratios"])
		ax.scatter(ratio, abund, c=color, s=total_abund/100, alpha=0.6, edgecolors="k", zorder=10)
		xs.append(ratio)
		ys.append(abund)
		if "_" not in taxa:
			ax.text(ratio*1.1, abund*1.02, taxa, fontsize=axis_font, zorder=11)
		if taxa=="T17_Micromonas_55_7":
			ax.text(ratio*1.1, abund*1.02, "Chromosomes", fontsize=axis_font, zorder=11)
		if taxa=="T17_Micromonas_42_105":
			ax.text(ratio*1.1, abund*1.02, "Chloroplast", fontsize=axis_font, zorder=11)
		if taxa=="T17_Micromonas_48_77":
			ax.text(ratio*1.1, abund*1.02, "Mitochondria", fontsize=axis_font, zorder=11)

	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlim(min(xs)/1.5, max(xs)*6)
	ax.set_ylim(min(ys)/1.5, max(ys)*1.8)
	ax.set_xlabel("RNA:DNA representation ratio", fontsize=label_font)
	if "T17_Micromonas_55_7" not in taxon_info:
		ax.set_ylabel("Average DNA coverage of contigs (CPM)", fontsize=label_font)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	plt.tick_params(labelsize=axis_font)
	plt.grid(which='minor', alpha=0.3, linestyle='--')


def length_legend(ax):
	# make blank rectagle
	legend_elements = []
	for taxon in range(12):
		legend_elements.append(Patch(facecolor="w", edgecolor="w", label="    ", linewidth=1))
	ax.legend(handles=legend_elements, loc="lower right", framealpha=1, frameon=True, facecolor='w', ncol=2, columnspacing=1, handlelength=1, prop={'size': label_font})

	x=70
	y=0.18
	ax.text(x-10, y+0.2, "Contig length:", fontsize=label_font, verticalalignment='center', zorder=10)
	for s in [100000, 30000, 10000, 3000]:
		ax.scatter(x, y, s=s/500, alpha=0.5, c='b', zorder=10)
		ax.text(x*1.2, y, str(s/1000)+"kb", fontsize=label_font, verticalalignment='center', zorder=10)
		y/=2

def expression_legend(ax):
	# make blank rectagle
	plt.axis('off')
	legend_elements = []
	for taxon in range(8):
		legend_elements.append(Patch(facecolor="w", edgecolor="w", label="                     ", linewidth=1))
	ax.legend(handles=legend_elements, loc="lower center", framealpha=1, frameon=True, facecolor='w', ncol=2, columnspacing=1, handlelength=1, prop={'size': label_font})

	# make the circle legend
	x=0.4
	y=0.35
	ax.text(x-0.2, y+0.15, "Total transcript contribution:", fontsize=label_font, verticalalignment='center', zorder=10)
	for s in [30000, 10000, 3000, 1000]:
		ax.scatter(x, y, s=s/100, alpha=0.5, c='b', zorder=10)
		ax.text(x, y-0.18, str(s)+"tpm", fontsize=axis_font, horizontalalignment="right", verticalalignment='center', rotation=30, zorder=10)
		x+=0.1
	ax.set_xlim(0,1)
	ax.set_ylim(0,1)


def color_legend(ax):
	plt.axis("off")
	legend_elements = []
	for taxon in scheme:
		legend_elements.append(Patch(facecolor=scheme[taxon], edgecolor='k', label=taxon, linewidth=1))
	ax.legend(handles=legend_elements, loc="lower center", framealpha=1, frameon=True, facecolor='w', ncol=2, columnspacing=1, handlelength=1, prop={'size': label_font})

		
	

##################   SETUP PLOT   ######################
# main figure layout:
font = {'family': 'arial', 'weight': 'normal', 'size': 12}
plt.rc('font', **font)

sns.set_palette("colorblind")
fig = plt.figure(figsize=(6, 9), dpi=300)

axis_font=8
label_font=10
heading_font=16

scheme = {      "Euryarchaeota":        "white",
		"Nanoarchaeota":        "blue",
                "Proteobacteria":       "orange",
                "Bacteroidetes":        "g",
                "Cyanobacteria":        "r",
                "Actinobacteria":       "c",
		"Firmicutes":        	"pink",
                "Chlorophyta":          "black"  }

##################   DRAW BLOBPLOT   ######################
print "drawing DNA vs RNA contig blobplot..."
dna_data = load_abundances("abundance_dna.tpm")
gene_data = load_abundances("abundance_genes.tpm")
rna_data = combine_genes(gene_data)

taxonomy = load_taxonomy("contig_taxonomy.tab")
dna_data = make_taxonomy_colors(taxonomy, dna_data, depth=1)
rna_data = make_taxonomy_colors(taxonomy, rna_data, depth=1)
ax1 = fig.add_axes([0.1, 0.55, 0.85, 0.42])
draw_scatter(dna_data, rna_data, ax1, archaea=True)
draw_scatter(dna_data, rna_data, ax1, archaea=False)
length_legend(ax1)
ax1.annotate("A", xy=(-0.08, 1.01), xycoords="axes fraction", fontsize=heading_font)

##################   TAXONOMY REPRESENTATION   ######################
print "drawing RNA/DNA taxonomy representation plot..."
ax2 = fig.add_axes([0.1, 0.17, 0.4, 0.30])
taxon_representation_info = representation_by_taxonomy(dna_data, gene_data)
draw_averages (taxon_representation_info, ax2)
ax2.annotate("B", xy=(-0.15, 1.01), xycoords="axes fraction", fontsize=heading_font)

##################   MAG REPRESENTATION   ######################
print "drawing RNA/DNA MAG representation plot..."
mags = load_tabular("contigs.mags")
mag_taxonomy = load_tabular("MAGS.info", col=8)
mag_representation_info = representation_by_mag(dna_data, gene_data, mags, mag_taxonomy)
ax3 = fig.add_axes([0.55, 0.17, 0.4, 0.30])
draw_averages (mag_representation_info, ax3)
ax3.annotate("C", xy=(-0.15, 1.01), xycoords="axes fraction", fontsize=heading_font)

##################   BOTTOM LEGENDS   ######################
ax4 = fig.add_axes([0, 0, 0.5, 0.17])
expression_legend(ax4)
ax5 = fig.add_axes([0.45, 0, 0.5, 0.17])
color_legend(ax5)

##################   FINISHING PLOT   ######################

#plt.subplots_adjust(left=0.1, right=0.95, top=0.9, bottom=0.1)
plt.savefig("figure_2.png", dpi=300)
#plt.savefig("figure_S2.eps", dpi=300)
#plt.show()
