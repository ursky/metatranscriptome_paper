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

def clear_ax(ax):
	ax.xaxis.set_visible(False)
	ax.set_yticklabels([])
	ax.set_yticklabels([])
	for spine in ax.spines.values(): spine.set_visible(False)
	ax.yaxis.set_ticks_position('none')
	


def insert_png(filename, ax):
	folder = os.path.dirname(os.path.realpath(__file__))
	im = plt.imread(get_sample_data(folder + "/" + filename))
	ax.imshow(im)
	ax.set_title(filename.split(".")[0], fontsize=title_font)


def manual_png(img, x, y):
	folder = os.path.dirname(os.path.realpath(__file__))
	im = Image.open(folder + "/" + img)

	plot_h = fig.bbox.ymax; plot_w = 0.95*fig.bbox.xmax/2
	h = im.size[1]
	w = im.size[0]

	resize_ratio = min(plot_h/h, plot_w/w)
	h*=resize_ratio
	w*=resize_ratio

	im = im.resize((int(w), int(h)), Image.ANTIALIAS)
	im = np.array(im).astype(np.float) / 255
	plt.figimage(im, x, y*(plot_h-h))
	

def add_labels(labels):
	ax = fig.add_axes([0, 0, 1, 1])
	for i,label in enumerate(labels):
		x = float(i)/2
		x = x+0.05
		y = 0.9
		ax.text(x, y, label, fontsize=heading_font)
	ax.axis("off")
	



##################   SETUP PLOT   ######################

# main figure layout:
font = {'family': 'arial', 'weight': 'normal', 'size': 12}
plt.rc('font', **font)

sns.set_palette("colorblind")
fig = plt.figure(figsize=(8, 3.5), dpi=300)

axis_font=8
label_font=10
heading_font=16


##################   DRAW KRONAS   ######################

print "adding kronagram figures..."
manual_png("krona_dna.png", 0.05 * fig.bbox.xmax/2, 0.90)
manual_png("krona_rna.png", 1.0 * fig.bbox.xmax/2, 0.5)
add_labels(["A", "B"])


##################   FINISHING PLOT   ######################

#plt.subplots_adjust(left=0.1, right=0.95, top=0.9, bottom=0.1)
plt.savefig("figure_1.png", dpi=300)
#plt.savefig("figure_S2.eps", dpi=300)
#plt.show()
