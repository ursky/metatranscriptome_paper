#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import operator


fig, ax = plt.subplots(figsize=(6, 6))

increment=0.2
for filename in sys.argv[1:]:
	data=[]
	n=0
	for i in np.arange(0.0, 14, increment):
		data.append(0)

	# parse the protein stats file
	for line in open(filename):
		cut =line.strip().split("\t")
		if len(cut)<5 or "Hydrophobicity" in line:
			continue
		#MW = float(cut[3])
		#if MW<10000:
		#	continue
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
	ax.plot(np.arange(0.0, 14, increment), data, linewidth=3, label=filename.split(".")[0])




plt.xlim([3, 13])
plt.grid(b=True, which='both', color='0.65', linestyle='-')
ax.set_xlabel("Estimated isoelectric point")
ax.set_ylabel("Relative gene abundance")
plt.legend()
plt.savefig("isoelectric_points.png", dpi=300, bbox_inches="tight")
plt.show()

