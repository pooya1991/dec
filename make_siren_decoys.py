import os
import sys
import numpy as np
import scipy.sparse as sp

# Reads in a matrix and reverses the order of the non-zero values
# This is assuming that each column only has a few non-zero values, or else
# this might be inefficient
# Also assumes that the rows appear in ascending order
#matfile = '/net/noble/vol3/user/alexhu/proj/Tsou2016/data/2017-07-05/ms1/B_D140314_SGSDSsample1_R01_MHRM_T0_ms1_sparsebinned_singlescans/0_x.txt'
#fmatfile = '/net/noble/vol3/user/alexhu/proj/Tsou2016/data/2017-07-05/ms1/B_D140314_SGSDSsample1_R01_MHRM_T0_ms1_sparsebinned_singlescans_withdecoys/0_x.txt' # the fmatfile contains both the original and the decoys, where decoys are just appended

xmatdir = sys.argv[1]
if xmatdir[-1] != '/':
    xmatdir = xmatdir + '/'

for matfile in os.listdir(xmatdir):
	if not '_x.txt' in matfile:
		continue
	matfile = xmatdir + matfile
	fmatfile = matfile.replace('_x.txt','_xdecoy.txt')

	fin = open( matfile,'r' )
	dim = fin.readline()
	[M,N] = [ int(v) for v in dim.rstrip().split() ]

	intensities = [ [] for n in range(N) ] # Each list is list of the values for each precursor n
	for line in fin:
		if line[0] == '>':
			continue
		else:
			tokens = line.rstrip().split()
			n = int(tokens[0])
			intensity = tokens[1]	
			intensities[n].append(intensity)
	fin.close()

	# Reverse the intensities
	for n in intensities:
		n.reverse()

	fout = open(fmatfile,'w')
	fin = open(matfile,'r')
	fin.readline()
	fout.write( str(M) + '\t' + str(N*2) + '\n' )
	decoylines = []
	dixes = [0]*N
	for line in fin:
		if '>' == line[0]:
			for dl in decoylines:
				fout.write(dl)
			decoylines = []
			fout.write(line)
		else:
			tokens = line.rstrip().split()
			n = int(tokens[0])
			dn = n + N
			decoylines.append( str(dn) + '\t' + intensities[n][dixes[n]] + '\n' )
			dixes[n] += 1
			fout.write( line )
	fin.close()
	for dl in decoylines:
		fout.write(dl)		
	fout.close()

