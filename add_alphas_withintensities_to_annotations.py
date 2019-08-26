import os
import sys

execfile('read_rowmajor_matrix.py')

xdir = sys.argv[1]
annfile = sys.argv[2]
pathdir = sys.argv[3]

def get_numds( targets, decoys ):	
	numds = []
	targets.sort( reverse = True )
	decoys.sort( reverse = True )

	ti = 0
	di = 0

	nt = len(targets)
	nd = len(decoys)

	while ti < nt and di < nd:
		if targets[ti] > decoys[di]:
			ti += 1
			numds.append(di)
		else:
			di += 1
	return numds
			
# This assumes that the scan number is recorded with the binbounds 
# And that the bins are presented in order
def find_bin_blocks( binbounds ):
	t = binbounds[0][2]
	starts = { t:0 }
	ends = {}
	for i in range(1,len(binbounds)):
		if binbounds[i][2] == -1:
			continue
		if binbounds[i][2] != t:
			ends[t] = i	
			t = binbounds[i][2]
			starts[t] = i
	ends[t] = i

	bounds = { (t,starts[t],ends[t]) for t in starts }
	return bounds

# Give the dimensions of the target blocks for the scans
def get_nbounds( xfile, decoys=False ):
	binbounds = get_binbounds(xfile)
	bounds = find_bin_blocks( binbounds )
	fin = open(xfile,'r')
	dim = fin.readline().rstrip().split()
	dim = ( int(dim[0]), int(dim[1]) )
	tN = dim[1] # number of targets
	if decoys:
		tN = tN/2

	bounds = sorted(list(bounds))

	# Get the first column that is nonzero for a given index
	first_columns = [-1]*dim[0]
	record = False
	rownum = 0
	for line in fin:
		if '>' in line:
			rownum = int( line.rstrip().split()[1]	)
			record = True
			continue
		elif record:
			record = False
			cnum = int( line.rstrip().split()[0] )
			first_columns[rownum] = cnum
	fin.close()
	nstarts = { (t,first_columns[r1]) for (t,r1,r2) in bounds }
	return (tN, nstarts)

# Steps
# 1. Identify the elution peaks from the OLS solution
def read_alphas( f ):
	alphas = []
	fin = open(f)
	for line in fin:
		alphas.append( float( line.rstrip() ) )
	fin.close()
	return alphas

# The matrix in pathfile is assumed to be formatted such that
# the values are reported from left to right, top to bottom.
def first_nonzero_ix_per_index( pathfile ):
	fin = open(pathfile,'r')
	dim = fin.readline().rstrip().split()
	dim = ( int(dim[0]), int(dim[1]) )
	rownum = 0
	ixes = [-1]*dim[0]
	record = False
	for line in fin:
		if '>' == line[0]:
			rownum = int( line.rstrip().split()[1] )
			record = True
			continue
		if record:
			ixes[rownum] = int( line.rstrip().split()[0] )
			record = False
	fin.close()
	return ixes


if True:
	# Get the number of bins per scan
	xfiles = [ f for f in os.listdir(xdir) if 'xdecoy.txt' == f[-10:] ]
	xfiles.sort( key = lambda f: int(f.split('_')[0]) )

	Ms = {}
	Ns = [0]*2500
	totalNs = {}
	for xfile in xfiles:
		t = int( xfile.split('_')[0] )
		print 'Getting nbounds', t
		xfile = xdir + xfile

		#start_time = time.time()
		bb = get_binbounds(xfile)
		#el1 = time.time() - start_time
		

		#start_time = time.time()
		bounds = find_bin_blocks(bb)
		#el2 = time.time() - start_time


		#start_time = time.time()
		(N,ns) = get_nbounds(xfile,decoys=True)
		#el3 = time.time() - start_time
		totalNs[t] = N	
		#start_time = time.time()
		for (t,r1) in ns:
			Ns[t] = r1
		for (t,start,end) in bounds:
			Ms[t] = end-start
		#el4 = time.time() - start_time

		#print 'Times:', (el1, el2, el3, el4)
	
# Read in the elution peaks
if True:	
	anns = read_tab_delimited(annfile)
	numfeatures = len(anns['charge'])

def count_char( s, c ):
	i = 0
	s = str(s)
	for v in s:
		if v == c:
			i += 1
	return i

if True:
	# Peaks will be a 3-element tuple
	# Which is the elution time, the feature index within the scan number, and the index in the annfile
	peaklist = [] # This should be a tuple of ( peak scan, index, annfile index )
	for i in range(0,len(anns['charge'])):
		pts = anns['peak scans'][i]
		ints = anns['peak intensities'][i]
		if pts == '':
			continue

		ptslen = count_char(pts,',')+1
		intslen = count_char(ints,',')+1
		if ptslen != intslen:
			break

		if ptslen == 1:
			pts = [int(pts)]
			ints = [float(ints.split(':')[0])]
		else:
			pts = [int(float(p)) for p in pts.split(',')] 	
			ints = [float(p.split(':')[0]) for p in ints.split(',')]
		start_t = int(anns['start file'][i])
		end_t = int(anns['end file'][i])
		ixes = anns['indices'][i]

		if isinstance(ixes,int):
			ixes = [ixes]
		else:
			if ',' in ixes:
				ixes = [p for p in ixes.split(',')]
			else:
				ixes = [ixes]
		for (pt,intensity) in zip(pts,ints):
			if pt > end_t:
				#print 'Something is wrong'
				#print 'pt:', pt, (start_t,end_t)
				continue
			if pt < start_t:
				print 'Something is wrong'
				print 'pt:', pt, (start_t,end_t)
				continue
			ix = ixes[pt-start_t]
			if isinstance(ix,int):
				peaklist.append( (pt,ix,i,intensity) )
			elif '-' in ix:
				for v in ix.split('-'):
					peaklist.append( (pt,int(v),i,intensity) )
			else:
				peaklist.append( (pt,int(ix),i,intensity) )

	dpeaklist = [] # This should be a tuple of ( peak scan, index )
	for i in range(0,len(anns['charge'])):
		pts = anns['decoy peak scans'][i]
		ints = anns['decoy peak intensities'][i]
		if pts == '':
			continue

		ptslen = count_char(pts,',')+1
		intslen = count_char(ints,',')+1
		if ptslen != intslen:
			break

		if ptslen == 1:
			pts = [int(pts)]
			ints = [float(ints.split(':')[0])]
		else:
			pts = [int(float(p)) for p in pts.split(',')] 	
			ints = [float(p.split(':')[0]) for p in ints.split(',')]
		start_t = int(anns['start file'][i])
		end_t = int(anns['end file'][i])
		ixes = anns['indices'][i]

		if isinstance(ixes,int):
			ixes = [ixes]
		else:
			if ',' in ixes:
				ixes = [ p for p in ixes.split(',')] 	
			else:
				ixes = [ixes]
		for (pt,intensity) in zip(pts,ints):
			if pt > end_t or pt < start_t:
				continue
			ix = ixes[pt-start_t]
			if isinstance(ix,int):
				dpeaklist.append( (pt,ix,i,intensity) )
			elif '-' in ix:
				for v in ix.split('-'):
					dpeaklist.append( (pt,int(v),i,intensity) )
			else:
				dpeaklist.append( (pt,int(ix),i,intensity) )

	print 'Number of target peaks:', len(peaklist)
	print 'Number of decoy dpeaks:', len(dpeaklist)

if True:
	# separate peaks by peak time
	peakdict = {} # value is scan number, element is a tuple (precursor ix, annotation index)
	peaklist.sort()
	for p in peaklist:
		if not p[0] in peakdict:
			peakdict[p[0]] = []
		peakdict[p[0]].append( (p[1],p[2],p[3]) ) # precursor ix, ann ix, intensity
	
	# separate peaks by peak time
	dpeakdict = {}
	dpeaklist.sort()
	for p in dpeaklist:
		if not p[0] in dpeakdict:
			dpeakdict[p[0]] = []
		dpeakdict[p[0]].append( (p[1],p[2],p[3]) ) # precursor ix, ann ix, intensity

if True:
	# Read in information from the paths		

	# Both of these are a list of tuples, where each tuple is an elution peak time inferred
	# by Siren. The tuple contains: (alpha,index,time)
	targetalphas = []
	decoyalphas = []
	fs = [ f for f in os.listdir(pathdir) if not '.txt' in f ]
	for f in sorted( fs, key = lambda k: int(k) ):
	#for f in []:
		print 'f:',f
		# if not os.path.exists(pathdir+f+'/log.txt'):
		# 	continue
		for pf in sorted(os.listdir( pathdir + f )):
			if not 'path' in pf:
				continue
			# Get the maximum alpha for which each feature is non-zero
			tokens = pf.split('_')
			t = int( tokens[1] )
			print t			
			if not t in Ms:
				continue

			startn = int( tokens[2].split('t')[0] ) 	
			af = 'alphas_'+str(t)+'.txt' 
			alphas = np.array(read_alphas( pathdir + f + '/' + af ))
			alphas *= 2*Ms[t]	
	
			nzixes = first_nonzero_ix_per_index( pathdir + f + '/' + pf )
			nzalphas = [ alphas[n] for n in nzixes ]
			lnz = len(nzalphas)
			initial_n = int( pf.split('_')[2].split('t')[0] )
			peakalphas = []
			dpeakalphas = []

			# What are the LARS indices given the annotaFalseindices?
			if t in peakdict:
				peaks = peakdict[t] # Here are the annotated indices
				# The tuple is alpha, index, and peak scan
				for (p,i,intensity) in peaks:
					if p+Ns[t]-initial_n < lnz:
						peakalphas.append( (nzalphas[p+Ns[t]-initial_n],i,t,intensity) )
					else:
						print 'Something is wrong with the target indices', pf
						peakalphas.append( (0.0,i,t,intensity) )			

			if t in dpeakdict:
				dpeaks = dpeakdict[t] # Here are the annotated indices
				# The tuple is alpha, index, and peak scan
				for (p,i,intensity) in dpeaks:
					if p+Ns[t]-initial_n+totalNs[int(f)] < lnz:
						dpeakalphas.append( (nzalphas[p+Ns[t]-initial_n+totalNs[int(f)]],i,t,intensity) )
					else:
						print 'Something is wrong with the decoy indices', pf
						dpeakalphas.append( (0.0,i,t,intensity) )

			targetalphas += peakalphas
			decoyalphas += dpeakalphas	


def getintlist_commasep( s ):
	if s == '':
		return []
	if ',' in s:
		s = s.split(',')
		return [ int(float(v)) for v in s ]
	else:
		return [int(float(s))]

def tuple_to_string(t, numvals=0):
    s = str(t[0])
    if numvals == 0:
        numvals = len(s)
    for v in t[1:numvals]:
        s += ':' + str(v)
    return s

# Read in the annfile
# And output the same thing but with the alpha information
if True:
	print 'Arranging and outputting'
	# Arrange the alphas by ann index
	talphalist = [ [] for i in range(numfeatures) ]
	# The tuples in targetalphas are alpha, index, and peak scan
	for a in targetalphas:
		talphalist[a[1]].append( (a[0],a[2],a[3]) ) # alpha, peak scan, intensity
	dalphalist = [ [] for i in range(numfeatures) ]
	for a in decoyalphas:
		dalphalist[a[1]].append( (a[0],a[2],a[3] )) # alpha, peak scan, intensity
		if pt > end_t or pt < start_t:
			print 'Something is wrong'
			print 'pt:', pt, (start_t,end_t)
			continue

	alphaannfile = annfile.replace('.txt','_alphas_withintensities.txt')
	fout = open(alphaannfile,'w')
	fin = open(annfile,'r')
	header = fin.readline().rstrip().lstrip()
	htokens = header.split('\t')
	tscanix = htokens.index('peak scans')
	dscanix = htokens.index('decoy peak scans')
	fout.write(header + '\ttarget alphas\tdecoy alphas\n')
	i = 0
	for line in fin:
		tokens = line[:-1].split('\t')
		fout.write( line[:-1] + '\t' )	
		written = 0
		if len(talphalist[i]) > 0:
			# There might be duplicates for the same peak time
			# Choose the one with the highest alpha
			talphalist[i].sort( key = lambda a: (-a[0],a[1]) )
			alphas = [];
			times = set([])
			for a in talphalist[i]:
				if a[1] in times:
					continue
				alphas.append( a )
				times.add(a[1]) 
			fout.write(tuple_to_string(alphas[0], 2))
			#fout.write( str(alphas[0][0])+':'+str(alphas[0][1]) )
			written += 1
			for a in alphas[1:]:
				#fout.write( ','+str(a[0])+':'+str(a[1]) )
				fout.write(',' + tuple_to_string(a, 2))
				written += 1	
		# UGH when it doesn't work!

		existing = 0
		if tokens[tscanix] != '':
			if ',' in tokens[tscanix]:
				existing = len(tokens[tscanix].split(','))
			else:
				existing = 1
		if not existing == written:
			print 'Something is wrong'
			print 'Line number:', i
			print 'existing:', existing, 'written:', written
			print line	
			continue

		fout.write('\t')	
		if len(dalphalist[i]) > 0:
			dalphalist[i].sort( key = lambda a: (-a[0],a[1]) )
			alphas = [];
			times = set([])
			for a in dalphalist[i]:
				if a[1] in times:
					continue
				alphas.append( a )
				times.add(a[1]) 
			#fout.write( str(alphas[0][0])+':'+str(alphas[0][1]) )
			fout.write(tuple_to_string(alphas[0], 2))
			for a in alphas[1:]:
				#fout.write( ','+str(a[0])+':'+str(a[1]) )
				fout.write(',' + tuple_to_string(a, 2))
		fout.write('\n')
		i += 1			
	fin.close()
	fout.close()








