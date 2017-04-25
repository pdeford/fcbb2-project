#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr
import random

random.seed(425)

y = []
sequences = []

f = open("data/bdata.2009.mhci.public.1.txt")
f.readline()
for line in f:
	if line.startswith("human") and "HLA-A-0201" in line:
		fields = line.split()
		allele = fields[1]
		l = int(fields[2])
		sequence = fields[4]
		meas = float(fields[6])

		if l == 10:
			if meas < 150000:
				y.append(meas)
				sequences.append(sequence)

y = np.asarray(y)
sequences = np.asarray(sequences)
l = len(y)
print l
indices = range(l)
l = int(0.9*l)
random.shuffle(indices)
train, test = indices[:l], indices[l:]

func = sum

print "| K  | Full | 80% |"
print "|---:|:----:|:---:|"

plt.figure(figsize=[10,10])
for k in range(1,10):
	plt.subplot(3,3,k)
	plt.title("k = {}".format(k))

	# Trained on full data
	y2 = y[indices]
	kmers = {}
	for n in indices:
		seq = sequences[n]
		count = y[n]
		for i in range(10-k+1):
			kmer = seq[i:i+k]
			if kmer not in kmers:
				kmers[kmer] = count
			else:
				kmers[kmer] += count
	simple = [func([kmers[sequence[i:i+k]] for i in range(10-k+1)]) for sequence in sequences[indices]]

	plt.scatter(y2, simple, edgecolor='orange', facecolor='none', alpha=0.5)
	r1 = pearsonr(y2, simple)[0]


	# Trained on 80% of the data
	kmers = {}
	for n in train:
		seq = sequences[n]
		count = y[n]
		for i in range(10-k+1):
			kmer = seq[i:i+k]
			if kmer not in kmers:
				kmers[kmer] = count
			else:
				kmers[kmer] += count
	simple = []
	y2 = y[test]
	for seq in sequences[test]:
		count = [0]
		for i in range(10-k+1):
			kmer = seq[i:i+k]
			if kmer in kmers:
				count.append(kmers[kmer])
		simple.append(func(count))

	plt.scatter(y2, simple, edgecolor='blue', facecolor='none', alpha=0.5)
	r2 = pearsonr(y2, simple)[0]
	#print k, r1, r2
	print "| {:>2} | {:>1.3f} | {:>1.3f} |".format(k, r1, r2)

	plt.xticks([])
	plt.yticks([])

plt.tight_layout()
plt.savefig("output/kmer_enrichment.png")


"""
##############################################################
# I am seeing that if I train on all of the data, I can 
# predict that data reasonably well. However, if I train on 
# anything less, I am not seeing all of the kmers, and so I 
# can't make any predictions. I should do like that one SELEX
# paper, and use a Markov Model to estimate the kmers that
# aren't seen, to see if that improves my predictions any.
# (Because then I won't just be using 0's for kmers I haven't
# seen before, but instead can use a predicted value from the
# MM).
##############################################################
"""