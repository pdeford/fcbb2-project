#!/usr/bin/env python

B50 = {
	"A":{"A": 5, "R":-2, "N":-1, "D":-2, "C":-1, "Q":-1, "E":-1, "G": 0, "H":-2, "I":-1, "L":-2, "K":-1, "M":-1, "F":-3, "P":-1, "S": 1, "T": 0, "W":-3, "Y":-2, "V": 0, "B":-2, "J":-2, "Z":-1, "X":-1, "*":-5},
	"R":{"A":-2, "R": 7, "N":-1, "D":-2, "C":-4, "Q": 1, "E": 0, "G":-3, "H": 0, "I":-4, "L":-3, "K": 3, "M":-2, "F":-3, "P":-3, "S":-1, "T":-1, "W":-3, "Y":-1, "V":-3, "B":-1, "J":-3, "Z": 0, "X":-1, "*":-5},
	"N":{"A":-1, "R":-1, "N": 7, "D": 2, "C":-2, "Q": 0, "E": 0, "G": 0, "H": 1, "I":-3, "L":-4, "K": 0, "M":-2, "F":-4, "P":-2, "S": 1, "T": 0, "W":-4, "Y":-2, "V":-3, "B": 5, "J":-4, "Z": 0, "X":-1, "*":-5},
	"D":{"A":-2, "R":-2, "N": 2, "D": 8, "C":-4, "Q": 0, "E": 2, "G":-1, "H":-1, "I":-4, "L":-4, "K":-1, "M":-4, "F":-5, "P":-1, "S": 0, "T":-1, "W":-5, "Y":-3, "V":-4, "B": 6, "J":-4, "Z": 1, "X":-1, "*":-5},
	"C":{"A":-1, "R":-4, "N":-2, "D":-4, "C":13, "Q":-3, "E":-3, "G":-3, "H":-3, "I":-2, "L":-2, "K":-3, "M":-2, "F":-2, "P":-4, "S":-1, "T":-1, "W":-5, "Y":-3, "V":-1, "B":-3, "J":-2, "Z":-3, "X":-1, "*":-5},
	"Q":{"A":-1, "R": 1, "N": 0, "D": 0, "C":-3, "Q": 7, "E": 2, "G":-2, "H": 1, "I":-3, "L":-2, "K": 2, "M": 0, "F":-4, "P":-1, "S": 0, "T":-1, "W":-1, "Y":-1, "V":-3, "B": 0, "J":-3, "Z": 4, "X":-1, "*":-5},
	"E":{"A":-1, "R": 0, "N": 0, "D": 2, "C":-3, "Q": 2, "E": 6, "G":-3, "H": 0, "I":-4, "L":-3, "K": 1, "M":-2, "F":-3, "P":-1, "S":-1, "T":-1, "W":-3, "Y":-2, "V":-3, "B": 1, "J":-3, "Z": 5, "X":-1, "*":-5},
	"G":{"A": 0, "R":-3, "N": 0, "D":-1, "C":-3, "Q":-2, "E":-3, "G": 8, "H":-2, "I":-4, "L":-4, "K":-2, "M":-3, "F":-4, "P":-2, "S": 0, "T":-2, "W":-3, "Y":-3, "V":-4, "B":-1, "J":-4, "Z":-2, "X":-1, "*":-5},
	"H":{"A":-2, "R": 0, "N": 1, "D":-1, "C":-3, "Q": 1, "E": 0, "G":-2, "H":10, "I":-4, "L":-3, "K": 0, "M":-1, "F":-1, "P":-2, "S":-1, "T":-2, "W":-3, "Y": 2, "V":-4, "B": 0, "J":-3, "Z": 0, "X":-1, "*":-5},
	"I":{"A":-1, "R":-4, "N":-3, "D":-4, "C":-2, "Q":-3, "E":-4, "G":-4, "H":-4, "I": 5, "L": 2, "K":-3, "M": 2, "F": 0, "P":-3, "S":-3, "T":-1, "W":-3, "Y":-1, "V": 4, "B":-4, "J": 4, "Z":-3, "X":-1, "*":-5},
	"L":{"A":-2, "R":-3, "N":-4, "D":-4, "C":-2, "Q":-2, "E":-3, "G":-4, "H":-3, "I": 2, "L": 5, "K":-3, "M": 3, "F": 1, "P":-4, "S":-3, "T":-1, "W":-2, "Y":-1, "V": 1, "B":-4, "J": 4, "Z":-3, "X":-1, "*":-5},
	"K":{"A":-1, "R": 3, "N": 0, "D":-1, "C":-3, "Q": 2, "E": 1, "G":-2, "H": 0, "I":-3, "L":-3, "K": 6, "M":-2, "F":-4, "P":-1, "S": 0, "T":-1, "W":-3, "Y":-2, "V":-3, "B": 0, "J":-3, "Z": 1, "X":-1, "*":-5},
	"M":{"A":-1, "R":-2, "N":-2, "D":-4, "C":-2, "Q": 0, "E":-2, "G":-3, "H":-1, "I": 2, "L": 3, "K":-2, "M": 7, "F": 0, "P":-3, "S":-2, "T":-1, "W":-1, "Y": 0, "V": 1, "B":-3, "J": 2, "Z":-1, "X":-1, "*":-5},
	"F":{"A":-3, "R":-3, "N":-4, "D":-5, "C":-2, "Q":-4, "E":-3, "G":-4, "H":-1, "I": 0, "L": 1, "K":-4, "M": 0, "F": 8, "P":-4, "S":-3, "T":-2, "W": 1, "Y": 4, "V":-1, "B":-4, "J": 1, "Z":-4, "X":-1, "*":-5},
	"P":{"A":-1, "R":-3, "N":-2, "D":-1, "C":-4, "Q":-1, "E":-1, "G":-2, "H":-2, "I":-3, "L":-4, "K":-1, "M":-3, "F":-4, "P":10, "S":-1, "T":-1, "W":-4, "Y":-3, "V":-3, "B":-2, "J":-3, "Z":-1, "X":-1, "*":-5},
	"S":{"A": 1, "R":-1, "N": 1, "D": 0, "C":-1, "Q": 0, "E":-1, "G": 0, "H":-1, "I":-3, "L":-3, "K": 0, "M":-2, "F":-3, "P":-1, "S": 5, "T": 2, "W":-4, "Y":-2, "V":-2, "B": 0, "J":-3, "Z": 0, "X":-1, "*":-5},
	"T":{"A": 0, "R":-1, "N": 0, "D":-1, "C":-1, "Q":-1, "E":-1, "G":-2, "H":-2, "I":-1, "L":-1, "K":-1, "M":-1, "F":-2, "P":-1, "S": 2, "T": 5, "W":-3, "Y":-2, "V": 0, "B": 0, "J":-1, "Z":-1, "X":-1, "*":-5},
	"W":{"A":-3, "R":-3, "N":-4, "D":-5, "C":-5, "Q":-1, "E":-3, "G":-3, "H":-3, "I":-3, "L":-2, "K":-3, "M":-1, "F": 1, "P":-4, "S":-4, "T":-3, "W":15, "Y": 2, "V":-3, "B":-5, "J":-2, "Z":-2, "X":-1, "*":-5},
	"Y":{"A":-2, "R":-1, "N":-2, "D":-3, "C":-3, "Q":-1, "E":-2, "G":-3, "H": 2, "I":-1, "L":-1, "K":-2, "M": 0, "F": 4, "P":-3, "S":-2, "T":-2, "W": 2, "Y": 8, "V":-1, "B":-3, "J":-1, "Z":-2, "X":-1, "*":-5},
	"V":{"A": 0, "R":-3, "N":-3, "D":-4, "C":-1, "Q":-3, "E":-3, "G":-4, "H":-4, "I": 4, "L": 1, "K":-3, "M": 1, "F":-1, "P":-3, "S":-2, "T": 0, "W":-3, "Y":-1, "V": 5, "B":-3, "J": 2, "Z":-3, "X":-1, "*":-5},
	"B":{"A":-2, "R":-1, "N": 5, "D": 6, "C":-3, "Q": 0, "E": 1, "G":-1, "H": 0, "I":-4, "L":-4, "K": 0, "M":-3, "F":-4, "P":-2, "S": 0, "T": 0, "W":-5, "Y":-3, "V":-3, "B": 6, "J":-4, "Z": 1, "X":-1, "*":-5},
	"J":{"A":-2, "R":-3, "N":-4, "D":-4, "C":-2, "Q":-3, "E":-3, "G":-4, "H":-3, "I": 4, "L": 4, "K":-3, "M": 2, "F": 1, "P":-3, "S":-3, "T":-1, "W":-2, "Y":-1, "V": 2, "B":-4, "J": 4, "Z":-3, "X":-1, "*":-5},
	"Z":{"A":-1, "R": 0, "N": 0, "D": 1, "C":-3, "Q": 4, "E": 5, "G":-2, "H": 0, "I":-3, "L":-3, "K": 1, "M":-1, "F":-4, "P":-1, "S": 0, "T":-1, "W":-2, "Y":-2, "V":-3, "B": 1, "J":-3, "Z": 5, "X":-1, "*":-5},
	"X":{"A":-1, "R":-1, "N":-1, "D":-1, "C":-1, "Q":-1, "E":-1, "G":-1, "H":-1, "I":-1, "L":-1, "K":-1, "M":-1, "F":-1, "P":-1, "S":-1, "T":-1, "W":-1, "Y":-1, "V":-1, "B":-1, "J":-1, "Z":-1, "X":-1, "*":-5},
	"*":{"A":-5, "R":-5, "N":-5, "D":-5, "C":-5, "Q":-5, "E":-5, "G":-5, "H":-5, "I":-5, "L":-5, "K":-5, "M":-5, "F":-5, "P":-5, "S":-5, "T":-5, "W":-5, "Y":-5, "V":-5, "B":-5, "J":-5, "Z":-5, "X":-5, "*": 1},
}

AAs = B50.keys()
B50_rows = dict([(AA, [B50[AA][aa] for aa in AAs]) for AA in AAs])

basic_characters = {
	#AA: [non-polar, polar, charged]
	"A":[1, 0, 0],
	"R":[0, 0, 1],
	"N":[0, 1, 0],
	"D":[0, 0, 1],
	"C":[0, 1, 0],
	"Q":[0, 1, 0],
	"E":[0, 0, 1],
	"G":[1, 0, 0],
	"H":[0, 1, 0],
	"I":[1, 0, 0],
	"L":[1, 0, 0],
	"K":[0, 0, 1],
	"M":[1, 0, 0],
	"F":[0, 0, 0],
	"P":[1, 0, 0],
	"S":[0, 1, 0],
	"T":[0, 1, 0],
	"W":[0, 1, 0],
	"Y":[0, 1, 0],
	"V":[1, 0, 0],
}

pKas = {	
	#"AA":[ pI, NH2 pKa, CO2H pKa, R pKa],
	"A":[ 6.02,  9.87, 2.35, 50.00],
	"R":[10.76,  8.99, 1.82, 12.48],
	"N":[ 5.41,  8.72, 2.14, 50.00],
	"D":[ 2.87,  9.90, 1.99,  3.90],
	"C":[ 5.14, 10.70, 1.92,  8.37],
	"Q":[ 5.65,  9.13, 2.17, 50.00],
	"E":[ 3.22,  9.47, 2.10,  4.07],
	"G":[ 5.97,  9.78, 2.35, 50.00],
	"H":[ 7.58,  9.33, 1.80,  6.04],
	"I":[ 6.02,  9.76, 2.32, 50.00],
	"L":[ 5.98,  9.74, 2.33, 50.00],
	"K":[ 9.74,  9.06, 2.16, 10.48],
	"M":[ 5.75,  9.28, 2.13, 50.00],
	"F":[ 5.48,  9.31, 2.20, 50.00],
	"P":[ 6.10, 10.64, 1.95, 50.00],
	"S":[ 5.68,  9.21, 2.19, 50.00],
	"T":[ 6.53,  9.10, 2.09, 50.00],
	"W":[ 5.88,  9.41, 2.46, 50.00],
	"Y":[ 5.65,  9.21, 2.20, 10.46],
	"V":[ 5.97,  9.74, 2.29, 50.00],
}

helical = {
	#AA: [svalue, wvalue, -RTln(w), ddG]
	"A": [1.540,	1.610,	-0.258,	-1.88],
	"R": [1.100,	1.200,	-0.047,	-1.67],
	"L": [0.920,	0.960,	 0.022,	-1.60],
	"K": [0.780,	0.960,	 0.108,	-1.52],
	"E": [(0.63+0.43)/2, (0.45+0.63)/2, (0.433+0.225)/2, (-1.20+-1.37)/2],
	"M": [0.600,	0.630,	 0.251,	-1.37],
	"Q": [0.530,	0.560,	 0.314,	-1.31],
	"I": [0.420,	0.440,	 0.445,	-1.18],
	"Y": [(0.370+0.5)/2,	(0.39+0.53)/2,	(0.344+0.511)/2,	(-1.28+-1.11)/2],
	"H": [0.360,	0.380,	 0.525,	-1.10],
	"S": [0.360,	0.380,	 0.525,	-1.10],
	"C": [0.330,	0.350,	 0.570,	-1.06],
	"N": [0.290,	0.310, 	 0.63,5	-1.00],
	"D": [0.290,	0.310, 	 0.63,5	-1.00],
	"W": [(0.290+0.360)/2,	(0.30+0.38)/2,	(0.525+0.653)/2,	(-1.10+-0.97)/2],
	"F": [0.280,	0.290, 	 0.67,2	-0.95],
	"V": [0.220,	0.230,	 0.797,	-0.83],
	"T": [0.130,	0.140,	 1.070,	-0.56],
	"H": [0.060,	0.060,	 1.530,	-0.10],
	"G": [0.050,	0.050,	 1.620,	 0.00],
	"P": [0.001,	0.001,	 4.000,	 5.00],
}

sigma_properties = {
	#"AA": [ MW, residue MW, CO2H pKa, NH2 pKa, R pKa, pI],
	"A" :  [ 89.10,  71.08, 2.34,  9.69, 50.00,  6.00],
	"R" :  [174.20, 156.19, 2.17,  9.04, 12.48, 10.76],
	"N" :  [132.12, 114.11, 2.02,  8.80, 50.00,  5.41],
	"D" :  [133.11, 115.09, 1.88,  9.60,  3.65,  2.77],
	"C" :  [121.16, 103.15, 1.96, 10.28,  8.18,  5.07],
	"E" :  [147.13, 129.12, 2.19,  9.67,  4.25,  3.22],
	"Q" :  [146.15, 128.13, 2.17,  9.13, 50.00,  5.65],
	"G" :  [ 75.07,  57.05, 2.34,  9.60, 50.00,  5.97],
	"H" :  [155.16, 137.14, 1.82,  9.17,  6.00,  7.59],
	"O" :  [131.13, 113.11, 1.82,  9.65, 50.00, 50.00],
	"I" :  [131.18, 113.16, 2.36,  9.60, 50.00,  6.02],
	"L" :  [131.18, 113.16, 2.36,  9.60, 50.00,  5.98],
	"K" :  [146.19, 128.18, 2.18,  8.95, 10.53,  9.74],
	"M" :  [149.21, 131.20, 2.28,  9.21, 50.00,  5.74],
	"F" :  [165.19, 147.18, 1.83,  9.13, 50.00,  5.48],
	"P" :  [115.13,  97.12, 1.99, 10.60, 50.00,  6.30],
	"U" :  [139.11, 121.09, 50.0, 50.00, 50.00,  5.68],
	"S" :  [105.09,  87.08, 2.21,  9.15, 50.00,  5.68],
	"T" :  [119.12, 101.11, 2.09,  9.10, 50.00,  5.60],
	"W" :  [204.23, 186.22, 2.83,  9.39, 50.00,  5.89],
	"Y" :  [181.19, 163.18, 2.20,  9.11, 10.07,  5.66],
	"V" :  [117.15,  99.13, 2.32,  9.62, 50.00,  5.96],
}

hydrophobicity = {
	#"AA":[ pH2,   pH7 ],
	"L" : [ 100,    97 ],
	"I" : [ 100,    99 ],
	"F" : [  92,   100 ],
	"W" : [  84,    97 ],
	"V" : [  79,    76 ],
	"M" : [  74,    74 ],
	"C" : [  52,    49 ],
	"Y" : [  49,    63 ],
	"A" : [  47,    41 ],
	"T" : [  13,    13 ],
	"E" : [   8,   -31 ],
	"G" : [   0,     0 ],
	"S" : [  -7,    -5 ],
	"Q" : [ -18,   -10 ],
	"D" : [ -18,   -55 ],
	"R" : [ -26,   -14 ],
	"K" : [ -37,   -23 ],
	"N" : [ -41,   -28 ],
	"H" : [ -42,     8 ],
	"P" : [ -46,   -46 ],
}

netmhc_surface = {
	#AA: [A, B, C, D, E, F],
	"G": [1, 0, 0, 0, 0, 0],
	"A": [1, 0, 0, 0, 0, 0],
	"S": [1, 0, 0, 0, 0, 0],
	"C": [0, 1, 0, 0, 0, 0],
	"T": [0, 1, 0, 0, 0, 0],
	"D": [0, 1, 0, 0, 0, 0],
	"V": [0, 1, 0, 0, 0, 0],
	"P": [0, 0, 1, 0, 0, 0],
	"N": [0, 0, 0, 1, 0, 0],
	"L": [0, 0, 0, 1, 0, 0],
	"I": [0, 0, 0, 1, 0, 0],
	"Q": [0, 0, 0, 1, 0, 0],
	"M": [0, 0, 0, 1, 0, 0],
	"E": [0, 0, 0, 1, 0, 0],
	"H": [0, 0, 0, 1, 0, 0],
	"K": [0, 0, 0, 0, 1, 0],
	"F": [0, 0, 0, 0, 1, 0],
	"R": [0, 0, 0, 0, 1, 0],
	"Y": [0, 0, 0, 0, 1, 0],
	"W": [0, 0, 0, 0, 0, 1],
}

AAs = "ACDEFGHIKLMNPQRSTVWY"
AA_index = dict(zip(AAs, range(20)))
bin_AA = {}
for key in AAs:
	bin_AA[key] = [1 if x == key else 0 for x in AAs]


def get_matrices():
	return B50_rows, basic_characters, pKas, helical, sigma_properties, hydrophobicity, bin_AA, netmhc_surface

def get_B50():
	return B50