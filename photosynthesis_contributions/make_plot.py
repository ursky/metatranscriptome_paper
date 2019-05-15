#!/usr/bin/env python
import sys
import numpy as np
import pandas as pd
import matplotlib
import seaborn as sns; sns.set(color_codes=True)
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
sns.set_palette("colorblind")
plt.rc('font', family='arial')



def load_contigs(filename, out={}):
	for line in open(filename):
		if line[0]==">":
			out[line[1:-1]]=filename.split(".")[0]
	return out


def load_expression(filename):
	out={}
	for line in open(filename):
		cut=line.strip().split("\t")
		if "_length_" in line:
			out[cut[0]] = float(cut[1])
	return out


def get_pathways(filename, taxonomy, expression):
	data={"Halothece":{}, "Euhalothece":{}, "Chloroplast":{}}
	for line in open(filename):
		cut=line.strip().split("\t")
		gene=cut[0]
		contig=gene.split("-")[0]
		protein=cut[2]
		paths = cut[4]
		if "Photosynthesis" not in paths or "photosystem" not in protein or "PsbK" in protein:
			continue
		if contig not in taxonomy:
			continue
		taxa=taxonomy[contig]
		if protein not in data[taxa]:
			data[taxa][protein]=0
		data[taxa][protein]+=expression[gene]
	return data


def make_expression_df(data):
	proteins=[]
	for taxa in data:
		for protein in data[taxa]:
			if protein not in proteins:
				proteins.append(protein)

	header=["Protein"]
	for taxa in data:
		header.append(taxa)
	#print "\t".join(header)
	
	data_dict={}
	for protein in sorted(proteins):
		data_dict[protein]={}
		out=[protein]
		for taxa in data:
			if protein in data[taxa]:
				out.append(str(data[taxa][protein]))
				data_dict[protein][taxa]=data[taxa][protein]
			else:
				out.append("NA")
				data_dict[protein][taxa]=None
		if "NA" in out:
			data_dict.pop(protein, None)
		#print "\t".join(out)
		
	return pd.DataFrame.from_dict(data_dict).T


def add_data_to_df(filename, final_df, df):
	for gene in df.index.values:
		for taxon in df:
			tpm = df.at[gene, taxon]+1
			row = pd.DataFrame([[filename, taxon, gene, tpm]], columns=['sample', 'organism', 'gene', 'tpm'])
			final_df = final_df.append(row, ignore_index=True)
	return final_df


def rename_labels(labels):
	for i,label in enumerate(labels):
		cut=label.split()
		out=[]
		for j,word in enumerate(cut):
			if word=="photosystem":
				if cut[j+1]=="I":
					out.append("PSI")
				else:
					out.append("PSII")
			elif word=="I" or word=="II":
				continue
			elif word=="cytochrome":
				out.append("cytoch")
			elif word=="chlorophyll":
				out.append("Chlo")
			elif word=="apoprotein" or word=="protein" or word=="a":
				continue
			elif word=="alpha":
				out[-2]=out[-2]+r'$\alpha$'
				out=out[:-1]
			elif word=="beta":
				out[-2]=out[-2]+r'$\beta$'
				out=out[:-1]
			else:
				out.append(word)
		labels[i]=" ".join(out)
	return labels


# MAIN

if "reload" in sys.argv:
	contigs = load_contigs("Halothece.fa")
	contigs = load_contigs("Euhalothece.fa", out=contigs)
	contigs = load_contigs("Chloroplast.fa", out=contigs)

	data={}
	final_df = pd.DataFrame(columns=['sample', 'organism', 'gene', 'tpm'])
	
	for time in ["9AM", "9PM"]:
		for i in range(1,7):
			filename=time+"-"+str(i)+"-mRNA.quant.counts"
			print "loading "+filename
			data[filename] = load_expression(filename)
	
			photo_genes = get_pathways("img_annotation.master", contigs, data[filename])
			df = make_expression_df(photo_genes)
			final_df = add_data_to_df(filename, final_df, df)
	final_df.to_pickle("data_df.pkl")
else:
	final_df = pd.read_pickle("data_df.pkl")



# MAKE THE PLOT
font = {'family': 'arial', 'weight': 'normal', 'size': 12}
sns.set_palette("colorblind")
sns.set_style("whitegrid")

fig, ax = plt.subplots(figsize=(10,5))
ax.set_yscale('log')
sns.catplot(x="gene", y="tpm", hue="organism", kind="bar", data=final_df, ax=ax)
sns.swarmplot(x="gene", y="tpm", hue="organism", edgecolor='gray', linewidth=0.3, size=3, data=final_df, split=True, ax=ax)

labels = [item.get_text() for item in ax.get_xticklabels()]
labels = rename_labels(labels)
ax.set_xticklabels(labels, rotation=45)


handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[3:], labels[3:], loc='upper right', ncol=3)


ax.set_ylabel("Total expression (Log10 TPM)")
ax.set_xlabel("Photosynthetic genes")

ax.set_ylim(0.7,3000)
#ax.grid(b=True, which='major', color='w', alpha=0.9, linestyle='-')
ax.yaxis.grid(b=True, which='major', color='k', alpha=0.2, linewidth=1, linestyle='--')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)


fig.tight_layout()
fig.savefig("figure_5.png", dpi=300)






