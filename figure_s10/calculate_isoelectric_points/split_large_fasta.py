#!/usr/bin/env python
import sys

max_chars=25000000
chars=0
file_ct=1
too_long=False
for line in open(sys.argv[1]):
	chars+=len(line)
	if chars>max_chars:
		too_long=True

	if line[0]==">" and too_long==True:
		file_ct+=1
		chars=0
		too_long=False
		print file_ct
	
	out = ".".join([ sys.argv[1].split(".")[0], "split", str(file_ct), "fa" ])
	f = open(out, "a")
	f.write(line)


