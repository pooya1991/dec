mport os
import sys

execfile('read_rowmajor_matrix.py')

xdir = sys.argv[1]
bdir = sys.argv[2]
annfile = sys.argv[3]

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
	nstarts = { t:first_columns[r1] for (t,r1,r2) in bounds }
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
	xfiles = [ f for f in os.listdir(xdir) if 'xdecoy.txt' in f ]
	xfiles.sort( key = lambda f: int(f.split('_')[0]) )

	Ms = {}
	Ns = {}
	totalNs = {}
	for xfile in xfiles:
		t = int( xfile.split('_')[0] )
		xfile = xdir + xfile
		bb = get_binbounds(xfile)
		bounds = find_bin_blocks(bb)
		print 'Getting nbounds', t
		(N,ns) = get_nbounds(xfile,decoys=True)
		Ns.update(ns)
		totalNs[t] = N	
		for (t,start,end) in bounds:
			Ms[t] = end-start
	print Ms		


	
# Read in the elution peaks
anns = read_tab_delimited(annfile)

peaks = [] # This should be a tuple of ( peak scan, index )
for i in range(0,len(anns['charge'])):
	pts = anns['peak scans'][i]
	if pts == '':
		continue
	if isinstance(pts,int) or isinstance(pts,float):
		pts = [int(pts)]
	else:
		pts = [int(float(p)) for p in pts.split(',')] 	
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
	for pt in pts:
		if pt >= end_t:
			continue
		ix = ixes[pt-start_t]
		if isinstance(ix,int):
			peaks.append( (pt,ix) )
		elif '-' in ix:
			for v in ix.split('-'):
				peaks.append( (pt,int(v)) )
		else:
			peaks.append( (pt,int(ix)) )

dpeaks = [] # This should be a tuple of ( peak scan, index )
for i in range(0,len(anns['charge'])):
	pts = anns['decoy peak scans'][i]
	if pts == '':
		continue
	if isinstance(pts,int) or isinstance(pts,float):
		pts = [int(pts)]
	else:
		pts = [int(float(p)) for p in pts.split(',')] 	
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
	for pt in pts:
		if pt >= end_t:
			continue
		ix = ixes[pt-start_t]
		if isinstance(ix,int):
			dpeaks.append( (pt,ix) )
		elif '-' in ix:
			for v in ix.split('-'):
				dpeaks.append( (pt,int(v)) )
		else:
			dpeaks.append( (pt,int(ix)) )

print 'Number of target peaks:', len(peaks)
print 'Number of decoy dpeaks:', len(dpeaks)

# separate peaks by peak time
peakdict = {}
peaks.sort()
for p in peaks:
	if not p[0] in peakdict:
		peakdict[p[0]] = []
	peakdict[p[0]].append( p[1] )
	
# separate peaks by peak time
dpeakdict = {}
dpeaks.sort()
for p in dpeaks:
	if not p[0] in dpeakdict:
		dpeakdict[p[0]] = []
	dpeakdict[p[0]].append( p[1] )

if True:
	# Read in information from the paths		

	targetalphas = []
	decoyalphas = []

	for f in sorted( os.listdir(pathdir), key = lambda k: int(k) ):
		if not os.path.exists(pathdir+f+'/log.txt'):
			continue
		print f
		for pf in os.listdir( pathdir + f ):
			if not 'path' in pf:
				continue

			# Get the maximum alpha for which each feature is non-zero
			tokens = pf.split('_')
			t = int( tokens[1] )
			
			if not t in Ms:
				continue

			startn = int( tokens[2].split('t')[0] ) 	
			af = 'alphas_'+str(t)+'.txt' 
			alphas = np.array(read_alphas( pathdir + f + '/' + af ))
			alphas *= 2*Ms[t]	
	
			nzixes = first_nonzero_ix_per_index( pathdir + f + '/' + pf )
			nzalphas = [ alphas[n] for n in nzixes ]
			initial_n = int( pf.split('_')[2].split('t')[0] )

			# What are the LARS indices given the annotated indices?
			if not t in peakdict:
				continue
			peaks = peakdict[t] # Here are the annotated indices
			peakalphas = [ nzalphas[p+Ns[t]-initial_n] for p in peaks ]

			if t in dpeakdict:			
				dpeaks = dpeakdict[t] # Here are the annotated indices
				dpeakalphas = [ nzalphas[p+Ns[t]-initial_n+totalNs[int(f)]] for p in peaks ]
						
			targetalphas += peakalphas
			decoyalphas += dpeakalphas	



if True:
	# Examine target decoy alphas

	numds = get_numds( targetalphas,decoyalphas )

	fdr = np.array(numds,dtype=float)/np.array(len(numds))
	qval = np.copy(fdr)
	for i in range(len(qval)-1,0,-1):
		if qval[i] < qval[i-1]:
			qval[i-1] = qval[i]

	plt.clf()
	plt.plot( qval, range(len(qval)) )
	plt.xlim( (0,0.1) )
	plt.xlabel('false discovery proportion')
	plt.ylabel('number of accepted elution peaks')
	plt.title('Number of Inferred Elution Peaks')
	#plt.savefig('target_decoy_alphas.png')
	plt.show()











