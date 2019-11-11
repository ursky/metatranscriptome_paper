#!/usr/bin/env python
import sys


alt_map = {'ins':'0'}
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

def reverse_complement(seq):
    for k,v in alt_map.iteritems():
        seq = seq.replace(k,v)
    bases = list(seq)
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    for k,v in alt_map.iteritems():
        bases = bases.replace(v,k)
    return bases


assembly={}; seq=""
for line in open(sys.argv[1]):
        if line[0]==">":
                if seq!="": assembly[contig]=seq
                contig=line.strip()[1:]
                seq=""
        else:
                seq+=line.strip()
assembly[contig]=seq

for line in open(sys.argv[2]):
        cut=line.strip().split('\t')
        st=int(cut[3])
        fi=int(cut[4])
        seq=assembly[cut[0]][st-1:fi]
        if cut[6]=='-': 
		seq=reverse_complement(seq)
        name=cut[8].split("locus_tag=")[-1].split(";")[0]
        print ">"+name
        print seq
