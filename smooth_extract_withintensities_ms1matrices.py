import os
import sys
import scipy.signal
import numpy as np

execfile('read_rowmajor_matrix.py')

bfile = sys.argv[1]
minlen = int(sys.argv[2])
afile = sys.argv[3]

print bfile
#sys.exit(1)

outfile = bfile.replace('.txt','_elutionpeaks.txt')
print 'Reading ' + bfile
b = read_rowmajor_matrix_sparse( bfile )
b = sp.csc_matrix(b)
clen = 7
# For this to work, the total length of "curve" can't be longer than minlen

print 'Smoothing'
P = b.shape[1]
Ptargets = P/2
T = b.shape[0]
lens = [0]*P
maxima = [ [] for p in xrange(P) ]
#for p in xrange(15):

smoothedout = open( bfile.replace('.txt','_smoothed.txt'), 'w' )
smoothedout.write( str(P) + '\t' + str(T) + '\n' )

for p in xrange(P):
#for p in xrange(63):
	if p % 2500 == 0:
		print 'Smoothing peptide ' + str(p)
	# Smooth
	prof = b[:,p].todense()
	nz = prof.nonzero()[0]
	if nz.shape[0] < 2:
		continue
	(mint,maxt) = (nz.min(),nz.max()+1)
	if maxt-mint < minlen:
		continue
	lens[p] = maxt-mint
	sub = np.asarray(prof[mint:maxt])[:, 0]

	# Pad the profile with 0s so it can be transformed with the
	# Savitzgy golay filter with 7 points with a 4th order polynomial
	padlen = 0
	newmint = mint
	if lens[p] < clen:
		diff = clen-lens[p]
		padlen = int(np.ceil( float(diff)/2 ))
		pad = np.zeros(padlen)
		sub = np.concatenate((pad, sub, pad))
		newmint = mint - padlen

	smoothed = scipy.signal.savgol_filter(sub, 7 ,4)
	smoothed[ smoothed < 0 ] = 0

	smoothedout.write( '>\t' + str(p) + '\n')
	for t in xrange(smoothed.shape[0]):
		if t + newmint < 0:
			continue
		smoothedout.write( str(t + newmint) + '\t' + str(smoothed[t]) + '\n' )

	smoothedlen = smoothed.shape[0]
	# Find Maxima
	if newmint >= 0 and smoothed[0] >= smoothed[1] and smoothed[0] != 0:
		# If the maximum of the smoothed points is even before the original length,
		# just make it the beginning of the original length
		maxima[p].append((mint, sub[mint - newmint]))

	for t in xrange(1,smoothedlen-1):
		if t + newmint < 0:
			continue
		if smoothed[t] > smoothed[t-1] and smoothed[t] > smoothed[t+1]:
			if t + newmint < mint:  # If the smoothed maximum scan precedes even the
				# original minimum scan, just give the original minimum scan
				maxima[p].append((mint, smoothed[t]))
			else:
				maxima[p].append((t + newmint, smoothed[t]))

	if smoothed[-1] >= smoothed[-2] and smoothed[-1]>0:
		maxima[p].append( (maxt-1,smoothed[-1]) )

smoothedout.close()

ints = []
for m in maxima:
	ints += [ n[1] for n in m ]
ints = np.array(ints)

ends = []

if True:
	numpeaks = 0
	# Add the maxima to the list of annotations
	apeakfile = bfile.replace('.txt','_intensities_annfile.txt')# Like afile, but contains the info about where the elution peak is

	fin = open(afile,'r')
	fout = open(apeakfile,'w')
	header = fin.readline()
	fout.write( header.replace('\n','\tpeak scans\tpeak intensities\tdecoy peak scans\tdecoy peak intensities\n') )	
	p = 0
	for line in fin:
		# Output the targets
		maxes = maxima[p]
		fout.write(line.rstrip() + '\t')
		if len(maxes) > 0:
			maxes.sort(key=lambda m: m[1], reverse=True)
			# Output peak times
			fout.write(str(maxes[0][0]))
			numpeaks += 1
			end = int(lens[p] / 4)
			ends.append(end)
			for (snum, i) in maxes[1:end]:
				fout.write(',' + str(snum))
				numpeaks += 1
			# Output peak intensities and their peak times
			fout.write('\t' + str(maxes[0][1]) + ':' + str(maxes[0][0]))
			for (snum, i) in maxes[1:end]:
				fout.write(',' + str(i) + ':' + str(snum))
		else:
			fout.write('\t')

		# Now output the decoys
		fout.write('\t')
		maxes = maxima[p + Ptargets]
		if len(maxes) > 0:
			maxes.sort(key=lambda m: m[1], reverse=True)
			# Output peak times
			fout.write(str(maxes[0][0]))
			numpeaks += 1
			end = int(lens[p] / 4)
			ends.append(end)
			for (snum, i) in maxes[1:end]:
				fout.write(',' + str(snum))
				numpeaks += 1
			fout.write('\t' + str(maxes[0][1]) + ':' + str(maxes[0][0]))
			for (snum, i) in maxes[1:end]:
				fout.write(',' + str(i) + ':' + str(snum))
		else:
			fout.write('\t')
		fout.write('\n')
		p += 1
	fin.close()
	fout.close()
	print 'Outputfile:', fout.name
	print 'Numpeaks: ' + str(numpeaks)
