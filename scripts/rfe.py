#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pickle
from sklearn.feature_selection import RFECV
from sklearn.linear_model import LogisticRegression as logit

from aafeatures import get_matrices


data_dicts = get_matrices()

X, y = [[] for x in range(len(data_dicts)+1)], []
sequences = []


counts = pickle.load(open("output/counts_4.p", "rb"))

flag = True

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
			for i,m in enumerate(data_dicts):
				if flag:
					X[i].append([])
				else:
					X[i].append([0 for x in m['A']])
				for aa in sequence:
					if flag:
						X[i][-1] += m[aa]
					else:
						X[i][-1] = [sum(x) for x in zip(m[aa], X[i][-1])]
			c = []
			for i in range(l-3):
				kmer = sequence[i:i+4]
				if kmer in counts:
					c.append(counts[kmer])
				else:
					c.append(0)
			X[-1].append(c)
			y.append(meas)
			sequences.append(sequence)

X = [np.asarray(x) for x in X]
y = np.asarray(y)
sequences = np.asarray(sequences)

weighted = []

ps = []
t = 0
totals = [t]
for x in X:
	n = x.shape[1]
	t += n
	p = n/10
	totals.append(t)
	ps.append(p)
	weighted.append(np.zeros([p,10]))


ps[-1] = 1
weighted[-1] = np.zeros([1,7]) - t


def get_feature(x):
	for i, t in enumerate(totals[:-1]):
		if t <= x < totals[i+1]:
			p = ps[i]
			position = (x - t)//ps[i]
			feature = (x - t)%ps[i]
			break
	else:
		print x, totals[-1]
	return i, position, feature


idx = np.argsort(y)
y = y[idx]
sequences = sequences[idx]


y[y < 500] = 0
y[y >= 500] = 1



clf1 = logit()

labels = ["B50_rows", "basic_characters", "pKas", "helical", "sigma_properties", 
		  "hydrophobicity", "bin_AA", "netmhc_surface", "count"]
data = dict(zip(labels,X))

X = np.hstack(X)[idx]

step = 10
cv = 3
rfecv = RFECV(estimator=clf1, step=step, cv=cv, scoring='roc_auc')
rfecv.fit(X, y)

print("Optimal number of features : %d / %d" % (rfecv.n_features_, X.shape[1]))


winners = np.arange(X.shape[1])[np.argsort(rfecv.ranking_)]

feature_labels = [
	["A","R","N","D","C","Q","E","G","H","I","L","K",
	 "M","F","P","S","T","W","Y","V","B","J","Z","X","*"],
	["non-polar", "polar", "charged"],
	["pI", "NH2 pKa", "CO2H pKa", "R pKa"],
	["svalue", "wvalue", "-RTln(w)", "ddG"],
	["MW", "residue MW", "CO2H pKa", "NH2 pKa", "R pKa", "pI"],
	[ "pH2",   "pH7" ],
	"ACDEFGHIKLMNPQRSTVWY",
	["A", "B", "C", "D", "E", "F"],
	["kmer"]
]

labels2 = dict(zip(labels, feature_labels))

positions = []
features = []
features2 = []

for i,x in enumerate(winners):
	a,b,c = get_feature(x)
	positions.append(b+1)
	features.append(a)
	features2.append(c)
	weighted[a][c][b] = -i
	#if i < 100:
	#	print "{:<16} {:>2} {}".format(labels[a], b, labels2[labels[a]][c],)

n = rfecv.n_features_

plt.figure()
plt.plot(range(1, len(positions)+1), positions, '.')
plt.xlabel('Ranked feature')
plt.ylabel('Position in amino acid')
plt.ylim([0,11])
plt.axvline(rfecv.n_features_, ls='--', label='Optimum', color='k')
plt.savefig("output/positions.png")
plt.close()

plt.figure()
plt.plot(range(1, len(positions)+1)[:n], positions[:n], '.')
plt.xlabel('Ranked feature')
plt.ylabel('Position in amino acid')
plt.ylim([0,11])
plt.axvline(rfecv.n_features_, ls='--', label='Optimum', color='k')
plt.savefig("output/positions2.png")
plt.close()

plt.figure()
plt.hist(positions[:n], bins=range(12))
plt.savefig('output/positions3.png')
plt.close()

plt.figure()
plt.plot(range(1, len(features)+1), features, '.')
plt.xlabel('Ranked feature')
plt.ylabel('Feature set')
plt.yticks(range(len(labels)), labels)
plt.ylim([-1,len(labels)])
plt.axvline(rfecv.n_features_, ls='--', label='Optimum', color='k')
plt.tight_layout()
plt.savefig("output/features.png")
plt.close()

plt.figure()
plt.plot(range(1, len(features)+1)[:n], features[:n], '.')
plt.xlabel('Ranked feature')
plt.ylabel('Feature set')
plt.yticks(range(len(labels)), labels)
plt.ylim([-1,len(labels)])
plt.axvline(rfecv.n_features_, ls='--', label='Optimum', color='k')
plt.tight_layout()
plt.savefig("output/features2.png")
plt.close()

# Plot number of features VS. cross-validation scores
plt.figure()
plt.xlabel("Number of features selected")
plt.ylabel("Cross validation score (auROC)")
plt.plot(range(1, len(rfecv.grid_scores_)*step + 1, step), rfecv.grid_scores_)
plt.axvline(rfecv.n_features_, ls='--', label='Optimum', color='k')
plt.legend(loc='lower right')
plt.savefig("output/rfe.png")
plt.close()

plt.figure(figsize=[10,10])
i = 0
for array in weighted:
	l1 = labels[i]
	l2 = labels2[l1]
	i += 1
	plt.subplot(3,3,i)
	plt.imshow(array, interpolation='nearest', cmap='Blues', aspect='auto', 
			   vmin=-max(winners), vmax=-min(winners))
	plt.xlabel("Position")
	plt.yticks(range(len(l2)),l2)
	plt.title(l1)
plt.colorbar()
plt.tight_layout()
plt.savefig("output/positionweights.png")
plt.close()