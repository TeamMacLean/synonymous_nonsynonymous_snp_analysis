#!/usr/bin/env python

import os, sys

from Bio import SeqIO

seqfile=open(sys.argv[1])
seqidset=set()

longest_sequences={}
for record in SeqIO.parse(seqfile, 'fasta'):
	seqid=record.id.replace("_", " ").split()[0]
	if seqid not in seqidset:
		seqidset.add(seqid)
		seqlen=len(str(record.seq))
		#print("Added: ", seqid)
	
	length=len(str(record.seq))
	#print("length ", length)
	if length > seqlen:
		longest_sequences[record.description] = str(record.seq)
		seqlen=length

for key in longest_sequences.keys():
	print(">"+key)
	print(longest_sequences[key])







