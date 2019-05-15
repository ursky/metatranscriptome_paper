#!/usr/bin/env python
import sys, os
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import random
import scipy.stats
import pandas as pd
import seaborn as sns

def load_files(dirname, contigs):
	print "loading files from "+dirname
	data={}
	for filename in os.listdir(dirname):
		sample = "-".join(filename.split(".")[0].split("-")[:2])
		if "RNA" in filename:
			sample="RNA-"+sample
		else:
			sample="DNA-"+sample
		if not filename.endswith("quant.counts"):
			continue
		sample_data = load_abundances(dirname+"/"+filename, contigs)
		data[sample] = sample_data
	return data


def load_abundances(file_name, contigs):
	print file_name
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
		if contig.split("-")[0] not in contigs:
			continue
		abund=float(cut[n])
		data[contig] = abund
	return data


def load_tabs(filename):
	data={}
	for line in open(filename):
		if line[0]=="#":
			continue
		cut=line.strip().split("\t")
		if len(cut)<2:
			continue
		data[cut[0]]=cut[1]
	return data


def combine_genes(coverage):
	print "calculating extression across contig based on all genes"
	data={}
	for sample in coverage:
		print sample
		data[sample]={}
		temp={}
		for gene in coverage[sample]:
			tpm=coverage[sample][gene]
			contig=gene.split("-")[0]
			if contig not in temp:
				temp[contig]=[]
			temp[contig].append(tpm)
		for contig in temp:
			mean=np.mean(temp[contig])
			data[sample][contig]=mean
	return data



def load_viruses(filename, minlen):
	data=[]
	lengths=[]
	for line in open(filename):
		if line[0]==">":
			contig=line[1:-1]
			l = int(contig.split("_")[3])
			if l<minlen:
				continue
			lengths.append(l)
			data.append(contig)

	print "there are "+str(len(data)) + " viral contigs found"
	return data, lengths


def correlate_virus_to_host(virus_hosts, abunds):
	for virus in virus_hosts:
		for host in virus_hosts[virus]:
			#print virus, host
			v_abunds=[]
			h_abunds=[]
			for sample in abunds:
				v_abunds.append(abunds[sample][virus])
				h_abunds.append(abunds[sample][host])
			if sum(h_abunds)==0 or sum(v_abunds)==0:
				continue
			corr, pval = scipy.stats.pearsonr(v_abunds, h_abunds)
			if pval<0.1:
				print virus +"\t" + host + "\t" + str(pval)


def virus_and_host_cooccurrence(virus_hosts, abunds):
	data={}
	for sample1 in abunds:
		if sample1=="ALL": continue
		data[sample1]={}
		for sample2 in abunds:
			if sample2=="ALL": continue
			data[sample1][sample2]=0
			for virus in virus_hosts:
				for host in virus_hosts[virus]:
					if abunds[sample1][virus]>1 and abunds[sample2][host]>1:
						data[sample1][sample2]+=1
	df = pd.DataFrame.from_dict(data)
	print df


def load_virus_host_pairs(filename):
	data={}
	for line in open(filename):
		cut=line.strip().split("\t")
		virus = cut[0]
		host = "_".join(cut[1].split("_")[:6])
		if virus not in data:
			data[virus]={}
		if host not in data[virus]:
			data[virus][host]=1
		else:
			data[virus][host]+=1
	print "there are "+str(len(data))+" viruses that have a crispr hit to a host"
	return data


def virus_dna_vs_rna(viruses, dna_coverage, rna_coverage, hosts, sample, ax):
	dnas=[]
	rnas=[]
	sizes=[]
	colors=[]
	for virus in viruses:
		dna=dna_coverage[sample][virus]+0.01
		rna=rna_coverage[sample][virus]+0.01
		length=int(virus.split("_")[3])
		host=hosts[virus]
		if host=="Euryarchaeota": colors.append("b")
		elif host=="Proteobacteria": colors.append("c")
		elif host=="Cyanobacteria": colors.append("g")
		elif host=="Bacteroidetes": colors.append("gold")
		elif host=="Firmicutes": colors.append("m")
		else: colors.append("k")

		dnas.append(dna)
		rnas.append(rna)
		sizes.append((length-5000)*1.0/100+20)
	ax = draw_scatter(dnas, rnas, sizes, colors, "Virus DNA read depth", "Virus average transcription", sample, ax)


def virus_activity(viruses, rna_coverage, datatype, ax):
	data=[]
	for virus in viruses:
		actives=0
		for sample in dna_coverage:
			if "9" not in sample:
				continue
			rna=rna_coverage[sample][virus]+0.01
			if rna>0.1:
				actives+=1
		data.append(actives)
	n, bins, patches = ax.hist(data, 13)
	for spine in ax.spines.values(): spine.set_visible(False)
	ax.grid(linestyle='--', alpha=0.5, zorder=-100)
	ax.set_xlim(0,13)
        ax.set_xlabel('Occurrence count (number of sample)')
        ax.set_ylabel('Number of viruses')
	if datatype=="RNA":
        	ax.set_title("Number of samples in which viruses are active", fontsize=16)
	else:
		ax.set_title("Number of samples in which viruses are present", fontsize=16)



def draw_scatter(list1, list2, sizes, colors, label1, label2, title, ax):
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlim(min(list1)*0.5, max(list1)*2)
	ax.set_ylim(min(list2)*0.5, max(list2)*2)
	ax.scatter(list1, list2, s=sizes, alpha=0.2, c=colors)
	ax.set_xlabel(label1, fontsize=15)
	ax.set_ylabel(label2, fontsize=15)
	ax.set_title(title, fontsize=20)
	return ax


def draw_histogram(lengths, ax):
	#ax.set_xscale('log')
	n, bins, patches = ax.hist(lengths, 30)
	ax.set_xlabel('Contig length')
	ax.set_ylabel('Number of viral contigs')
	ax.set_title("Viral contig length distribution", fontsize=20)
	ax.set_ylim(0,20)
	ax.yaxis.set_major_locator(ticker.MultipleLocator(2))
	# get rid of the frame
	for spine in ax.spines.values(): spine.set_visible(False)
	ax.grid(linestyle='--', alpha=0.5, zorder=-100)


def compute_ratios(df):
	ratios = {}
	data = df.to_dict()
	for dna_sample in data:
		if "RNA" in dna_sample:
			continue
		sample = "-".join(dna_sample.split("-")[1:])
		ratios[sample]={}
		rna_sample = "RNA-" + sample
		for virus in data[dna_sample]:
			dna = data[dna_sample][virus]
			rna = data[rna_sample][virus]
			if dna==0:
				ratio=0
			else:
				ratio=rna/dna
			ratios[sample][virus]=ratio
	return pd.DataFrame.from_dict(ratios)
			
	

########################### MAIN #####################################
font = {'family': 'arial', 'weight': 'normal', 'size': 12}
plt.rc('font', **font)


############################## LOAD DATA ###############################
viruses, lengths = load_viruses("virus_contigs.fa", 5000)
hosts = load_tabs("virus_taxonomy.tab")

if "reload" in sys.argv:
	dna_coverage = load_files("DNA_QUANT", viruses)
	gene_coverage = load_files("GENE_QUANT", viruses)
	rna_coverage = combine_genes(gene_coverage)
	
	dna_df = pd.DataFrame.from_dict(dna_coverage)
	rna_df = pd.DataFrame.from_dict(rna_coverage)

	dna_df.to_pickle("df_dna.pkl")
	rna_df.to_pickle("df_rna.pkl")
else:
	dna_df = pd.read_pickle("df_dna.pkl")
	rna_df = pd.read_pickle("df_rna.pkl")


df = pd.concat([dna_df, rna_df], axis=1, sort=False)


#################### PROCESS DATA ###########################


#df = compute_ratios(df)

# log standardize
#df+=0.1
#df = df.apply(np.log)

# standardize by sample total coverage
df = df.div(df.sum(axis=0), axis=1)

# standardize by maximum value in each row
df = df.div(df.max(axis=1), axis=0)


##################   MAKE HEATMAP   ######################
print "making heatmap..."


sns.clustermap(df, figsize=(8,12), col_cluster=True, yticklabels=True, xticklabels=True, cmap="magma", method="weighted")


plt.savefig("heatmap.png", dpi=600, bbox_inches="tight")



