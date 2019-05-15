#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import operator
from mpl_toolkits.axes_grid1 import make_axes_locatable


def load_iso(filename):
	out={}
	for line in open(filename):
		cut=line.strip().split("\t")
		if "Pro_id" in line or len(cut)<4:
			continue
		out[cut[0]]=float(cut[-1])
	return out


####################### START PLOTTING ######################

font = {'family': 'arial', 'weight': 'normal', 'size': 12}
plt.rc('font', **font)
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10, 10))

####################### PANEL A: pI distributions  #########################
print "plotting algae pI distributions..."
increment=0.2
for filename in ["Dolichomastix.iso", "Dunaliella.iso", "Micromonas.iso", "Ostereococcus.iso"]:
	data=[]
	n=0
	for i in np.arange(0.0, 14, increment):
		data.append(0)

	# parse the protein stats file
	for line in open(filename):
		cut =line.strip().split("\t")
		if len(cut)<5 or "Hydrophobicity" in line:
			continue
		pI=float(cut[-1])
		pI_bin = int((1/increment)*pI)
		data[pI_bin]+=1
		n+=1
	print filename, n

	#std data:
	tot=np.sum(data)
	for i in range(len(data)):
		data[i] = 100.0*data[i]/tot

	#plt.scatter(np.arange(0.0, 14, 0.1), data[sample], s=5, c=c)
	ax1.plot(np.arange(0.0, 14, increment), data, linewidth=2, label=filename.split(".")[0])

ax1.set_xticks(np.arange(3, 14, 1))
ax1.set_xlim(3,13)
ax1.set_xlabel("Estimated isoelectric point")
ax1.set_ylabel("Relative gene abundance")
ax1.legend(fontsize=10)
ax1.grid(b=True, which='both', color='0.65', linestyle='--', alpha=0.2)
ax1.annotate("A", xy=(-0.12, 0.96), xycoords="axes fraction", fontsize=20)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

######################## PANEL B: pairwise gene pI comparison ###########################
print "plotting micromonas vs dolichomastix gene pIs..."

halite_iso=load_iso("Dolichomastix.iso")
query_iso=load_iso("Micromonas.iso")

x=[]; y=[]
found={}
for line in open("Micromonas_vs_Dolichomastix.blast"):
	cut=line.strip().split("\t")
	q_iso = query_iso[cut[0]]
	h_iso = halite_iso[cut[1]]
	if cut[1] in found:
		continue
	else:
		found[cut[1]]=None
	pident = float(cut[2])
	pcov = float(cut[3])
	#if pident<50 or pcov<50:
	#	continue
	x.append(h_iso)
	y.append(q_iso)
#	if q_iso>10.5 and h_iso<10.5:
#		print cut[1]

ax2.scatter(x,y, s=10, alpha=0.2, c='k')
ax2.set_xticks(np.arange(3, 14, 1))
ax2.set_yticks(np.arange(3, 14, 1))
ax2.set_xlim(3,13)
ax2.set_ylim(3,13)
ax2.set_xlabel("Dolichomastix estimated gene pI")
ax2.set_ylabel("Micromonas estimated gene pI")
ax2.plot(x,x, '--', c="k")
ax2.grid(which='both', color='0.65', linestyle='--', alpha=0.3, zorder=-1)
ax2.annotate("B", xy=(-0.13, 0.96), xycoords="axes fraction", fontsize=20)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)


##################### SALT ADAPTATION GENE LABELING ###############
print "plotting salt adaptation gene expression..."

goi = {}
for line in open("goi.list"):
	goi[line.strip()]=None

x=[]; y=[]; c=[]
goi_info={}
for line in open("Dolichomastix.master"):
	cut=line.strip().split("\t")
	gene = cut[0]
	tpm = float(cut[1])+0.1
	gc = float(cut[2])
	gene_name = cut[3]

	if "NODE" in line:
		color = "blue"
	if "CHLO" in line:
		color = "magenta"
	if "MITO" in line:
		color = "black"
	if "hypothetical" in gene_name:
		continue

	c.append(color)
	x.append(gc)
	y.append(tpm)

	if gene in goi:
		goi_info[gene]=(gc, tpm, gene_name, color)


ax3.set_yscale('log')
ax3.set_xlabel("GC content (%)")
ax3.set_ylabel("Average gene expression (TPM)")
ax3.scatter(x, y, s=5, alpha=0.05, c=c)

goi_y=[]
for gene in goi_info:
	ax3.scatter(goi_info[gene][0], goi_info[gene][1], s=5, alpha=1, c='k')
	goi_y.append(goi_info[gene][1])


# set up histograms
divider = make_axes_locatable(ax3)
axHistx = divider.append_axes("top", 0.5, pad=0, sharex=ax3)
axHisty = divider.append_axes("right", 0.5, pad=0, sharey=ax3)
axHistx.xaxis.set_tick_params(labelbottom=False)
axHisty.yaxis.set_tick_params(labelleft=False)
axHisty.set_yscale('log')

# draw histograms
xbins = np.arange(25, 75, 2)
ybins = np.logspace(-1, 4, 50)
axHistx.hist(x, bins=xbins)
axHisty.hist(y, bins=ybins, orientation='horizontal')
axHisty.hist(goi_y*100, bins=ybins, orientation='horizontal', color="k")


ax3.set_xlim(25, 75)
ax3.set_ylim(0.1, 10**4)
axHistx.set_xlim(25, 75)
axHisty.set_ylim(0.1, 10**4)

ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
axHistx.axis("off")
axHisty.axis("off")
ax3.annotate("C", xy=(-0.13, 1.1), xycoords="axes fraction", fontsize=20)





####################### PANEL D: organelle  pI distributions  #########################
print "plotting algae pI distributions..."
increment=0.5
for filename in ["Chromosomes.iso", "Mitochondria.iso", "Chloroplast.iso"]:
        data=[]
        n=0
        for i in np.arange(0.0, 14, increment):
                data.append(0)

        # parse the protein stats file
        for line in open(filename):
                cut =line.strip().split("\t")
                if len(cut)<5 or "Hydrophobicity" in line:
                        continue
                pI=float(cut[-1])
                pI_bin = int((1/increment)*pI)
                data[pI_bin]+=1
                n+=1
        print filename, n

        #std data:
        tot=np.sum(data)
        for i in range(len(data)):
                data[i] = 100.0*data[i]/tot

        #plt.scatter(np.arange(0.0, 14, 0.1), data[sample], s=5, c=c)
        ax4.plot(np.arange(0.0, 14, increment), data, linewidth=2, label=filename.split(".")[0])

ax4.set_xticks(np.arange(3, 14, 1))
ax4.set_xlim(3,13)
ax4.set_xlabel("Estimated isoelectric point")
ax4.set_ylabel("Relative gene abundance")
ax4.legend(fontsize=10)
ax4.grid(b=True, which='both', color='0.65', linestyle='--', alpha=0.2)
ax4.annotate("A", xy=(-0.12, 0.96), xycoords="axes fraction", fontsize=20)
ax4.spines['top'].set_visible(False)
ax4.spines['right'].set_visible(False)



##################### EDITING AND SAVING #########################

plt.tight_layout()
plt.savefig("figure.png", dpi=300, bbox_inches="tight")










