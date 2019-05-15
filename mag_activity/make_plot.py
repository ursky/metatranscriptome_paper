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

def load_lib_sizes(filename):
	print "loading library sizes..."
	libs={}
	for line in open(filename):
		if line.startswith("#"): continue
		lib=line.strip().split("\t")[0]
		reads=int(line.strip().split("\t")[1])
		libs[lib] = reads
	return libs


def load_data(filename):
	print "loading abundance data..."
	df=pd.read_csv(filename, sep='\t', index_col=0)

	# remove all 0 rows
	df = df[(df.T != 0).any()]

	# standardize columns by total sum in each column
	df = df.div(df.sum(axis=0), axis=1)
	df=1000000*df
	return df


def set_colors_to_timeline(df):
	print "adding colored labels..."
	lut=[]
	for sample in df.columns.values:
		if "2013-04" in sample: lut.append('m')
        	elif "2014-09" in sample: lut.append('r')
        	elif "2015-06" in sample: lut.append('b')
		elif "2015-12" in sample: lut.append('c')
        	elif "2016-02" in sample: lut.append('y')
        	elif "2017-02" in sample: lut.append('g')
		else: lut.append('w')
	return lut


def draw_clustermap(df, lut):
	print "drawing clustermap..."
	sns.set(font_scale=1)
	if lut!=False:
		g = sns.clustermap(df, figsize=(14,8), col_colors=lut, col_cluster=True, yticklabels=True, cmap="magma")
	else:
		g = sns.clustermap(df, figsize=(14,8), col_cluster=True, yticklabels=True, cmap="magma")
	plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90)
	plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)


def fix_sample_naming(df):
	dic = df.to_dict()
	new={}
	for k in dic:
	        if "RNA" in k:
	                name = "-".join(["RNA"]+k.split("-")[:-1])
	        else:
	                name = "DNA-"+k
		new[name]=dic[k]
	df = pd.DataFrame.from_dict(new)
	return df


dna_df = load_data("mag_abundance.tab")
rna_df = load_data("mag_activity.tab")

df = pd.concat([dna_df, rna_df], axis=1, sort=True)
df = df.drop(["ALL-mRNA", "ALL-DNA"], axis=1)
df = fix_sample_naming(df)
print df

# log standardize:
df+=0.01; df=np.log(df)

# standardize rows by maximum value in each row
#df = df.div(df.max(axis=1), axis=0)

lut=False
draw_clustermap(df, lut)

plt.savefig("figure.png", bbox_inches='tight', dpi=300)





