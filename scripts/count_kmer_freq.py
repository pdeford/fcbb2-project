#!/usr/bin/env python

import matplotlib.pyplot as plt
import pickle

f = open("data/human_proteome.fasta")

counts = {}

k = 4

seq = ""
for line in f:
	if line.startswith(">"):
		if seq == "": continue
		for i in range(len(seq) - k + 1):
			kmer = seq[i:i+k]
			if kmer not in counts:
				counts[kmer] =  1
			else:
				counts[kmer] += 1
		seq = ""
	else:
		seq += line.strip()

for i in range(len(seq) - k + 1):
	kmer = seq[i:i+k]
	if kmer not in counts:
		counts[kmer] =  1
	else:
		counts[kmer] += 1

plt.hist(counts.values(), bins=30, log=True)
plt.savefig("output/freqdist_k{}.png".format(k))

for pair in sorted(counts.items(), reverse=True, key = lambda x:x[1])[:10]:
	print pair

pickle.dump(counts, open('output/counts_{}.p'.format(k), 'wb'))
