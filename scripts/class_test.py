#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
from sklearn.cross_validation import cross_val_score
from sklearn.cross_validation import KFold
import pickle

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

idx = np.argsort(y)
y = y[idx]
sequences = sequences[idx]


y[y < 500] = 0
y[y >= 500] = 1

from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier as RFC
from sklearn.linear_model import LogisticRegression as logit

kf = KFold(y.shape[0], n_folds=10, )
clf1 = logit()
clf2 = RFC()
clf3 = SVC()

labels = ["B50_rows", "basic_characters", "pKas", "helical", "sigma_properties", "hydrophobicity", "bin_AA", "netmhc_surface", "count"]
clfs = ["logit", "RFC  ", "SVM  "]
data = dict(zip(labels,X))

for l in labels:
	print l
	x = data[l]
	x = x[idx]
	for j,clf in enumerate([clf1, clf2, clf3]):
		perf = cross_val_score(clf, x, y, cv=10, n_jobs=-1, scoring='roc_auc')
		print "  {}  {:0.3f}".format(clfs[j], np.mean(perf))

custom = ["bin_AA", "B50_rows", "sigma_properties"]
print "/".join(custom)
x = np.hstack([data[x] for x in custom])[idx]
for j,clf in enumerate([clf1, clf2, clf3]):
	perf = cross_val_score(clf, x, y, cv=10, n_jobs=-1, scoring='roc_auc')
	print "  {}  {:0.3f}".format(clfs[j], np.mean(perf))

print "All together"
X = np.hstack(X)[idx]
for j,clf in enumerate([clf1, clf2, clf3]):
	perf = cross_val_score(clf, X, y, cv=10, n_jobs=-1, scoring='roc_auc')
	print "  {}  {:0.3f}".format(clfs[j], np.mean(perf))

x = data['bin_AA'][idx]
perf = cross_val_score(clf1, x, y, cv=10, n_jobs=-1, scoring='roc_auc')
clf1.fit(x,y)
#print clf1.coef_.shape
pwm = clf1.coef_.reshape([-1,20])
#print pwm
AAs = "ACDEFGHIKLMNPQRSTVWY"
AA_index = dict(zip(AAs, range(20)))


chemistry_ColorScheme = dict((
	("GSTYC", "green"), #, "polar"),
	("NQ", "purple"), #, "neutral"),
	("KRH", "blue"), #, "basic"),
	("DE", "red"), #, "acidic"),
	("PAWFLIMV", "black"), #, "hydrophobic")
        ))

for key in chemistry_ColorScheme.keys():
	for n in key:
		chemistry_ColorScheme[n] = chemistry_ColorScheme[key]

taylor_ColorScheme = dict((
	('A', '#CCFF00'),
	('C', '#FFFF00'),
	('D', '#FF0000'),
	('E', '#FF0066'),
	('F', '#00FF66'),
	('G', '#FF9900'),
	('H', '#0066FF'),
	('I', '#66FF00'),
	('K', '#6600FF'),
	('L', '#33FF00'),
	('M', '#00FF00'),
	('N', '#CC00FF'),
	('P', '#FFCC00'),
	('Q', '#FF00CC'),
	('R', '#0000FF'),
	('S', '#FF3300'),
	('T', '#FF6600'),
	('V', '#99FF00'),
	('W', '#00CCFF'),
	('Y', '#00FFCC')
))

if flag:
	plt.figure(figsize=(pwm.shape[0]/10,1), frameon=False)
	for i in range(pwm.shape[0]):
		row = pwm[i]
		row = row - min(row)
		row /= sum(row)
		row = zip(AAs, row)
		row.sort(key=lambda x:x[1])
		y = 0
		for AA, p in row:
			plt.text(i + 0.5, y, AA, ha='center', va='baseline', size=int(p*72), 
				color=chemistry_ColorScheme[AA], family='monospace')
			y += p
	plt.xlim(0,pwm.shape[0])
	plt.axis('off')
	plt.savefig('output/logo.png', dpi=600)



# n = X.shape[1]
# m = int(n**0.5)
# p = n//m + 1

# plt.figure(figsize=[30,30])
# for i in range(n):
# 	plt.subplot(m, p, i+1)
# 	plt.hist(X[:,i][y==0], alpha=0.5, color='blue')
# 	plt.hist(X[:,i][y==1], alpha=0.5, color='orange')
# 	plt.xticks([])
# 	plt.yticks([])
# plt.savefig("output/hist_matrix.png")
# plt.close()

# from sklearn.decomposition import PCA
# pca = PCA(n_components=2)
# pca.fit(X)
# X2 = pca.transform(X)

# plt.figure(figsize=[30,30])
# plt.scatter(X2[:,0], X2[:,1], c=y, lw=0, cmap='cool')#cmap='spring')
# plt.savefig('output/class_pca.png')
# plt.close()
