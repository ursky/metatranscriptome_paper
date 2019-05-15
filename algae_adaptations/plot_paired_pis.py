#!/usr/bin/env python

import sys
import matplotlib.pyplot as plt
import numpy as np

def load_iso(filename):
	out={}
	for line in open(filename):
		cut=line.strip().split("\t")
		if "Pro_id" in line or len(cut)<4:
			continue
		out[cut[0]]=float(cut[-1])
	return out

halite_iso=load_iso("Dolichomastix.iso")
query_iso=load_iso("Micromonas.iso")

x=[]; y=[]
for line in open("Micromonas_vs_Dolichomastix.blast"):
	cut=line.strip().split("\t")
	q_iso = query_iso[cut[0]]
	h_iso = halite_iso[cut[1]]
	pident = float(cut[2])
	pcov = float(cut[3])
	#if pident<50 or pcov<50:
	#	continue
	x.append(h_iso)
	y.append(q_iso)

fig, ax = plt.subplots()
ax.scatter(x,y, s=10, alpha=0.2)
ax.set_xlim(3,14)
ax.set_ylim(3,14)
plt.xticks(np.arange(3, 15, 1))
plt.yticks(np.arange(3, 15, 1))
ax.set_xlabel("Dolichomastix estimated gene pI")
ax.set_ylabel("Micromonas estimated gene pI")
ax.plot(x,x, ':')
plt.grid()
plt.savefig("figure.png", dpi=300)
#plt.show()
