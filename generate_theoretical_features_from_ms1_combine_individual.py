execfile('read_rowmajor_matrix.py')
execfile('peakclusterforbinning.py')
execfile('make_averagine_models.py')

mass_accuracy = 6e-06
maxcharge = 6

def get_ix(n):
	return -1*n-1
def get_n(ix):
	return (ix+1)*-1
def get_footprint( x1, n ):
	nz = x1[:,n].nonzero()[0].tolist()
	nz = nz + [0]*max(0,4-len(nz))
	return (nz[0],nz[1],nz[2],nz[3])	

# This outputs a series of x and y matrices in sparse format
# where it's al
# The Y matrix is one long vector
# The X matrix is one big matrix,
class output_vector:

	def __init__( self, outfilesuffix ):

		print '\n\nOPENING ' + outfilesuffix + '\n\n'		

		self.M = 0
		self.N = 0
		self.suffix = outfilesuffix
		self.xout = open( outfilesuffix + '_x.txt', 'w' )
		self.yout = open( outfilesuffix + '_y.txt', 'w' )
		self.annout = open( outfilesuffix + '_annfile.txt', 'w' )
		self.annout.write('t\tmin m/z\tmax m/z\tcharge\n')


	# X starts out as 
	# Here, the Y doesn't even contain the binbounds or the row number
	# Because it's 1 dimensional! :) 
	def output_xy( self, x, features, cids, binbounds,t ):
		nrow = x.shape[0]
		ncol = x.shape[1]
		(I,J) = x.nonzero()
		previ = -1
		for n in xrange(0,I.shape[0]):
			i = I[n]
			if previ != i:
				self.xout.write( '>\t' + str(i+self.M) + '\t' + str(binbounds[i][0]) + '\t' + str(binbounds[i][1]) + '\t' + str(t) + '\n' )
			previ = i
			j = J[n]
			self.xout.write( str(j+self.N) + '\t' + str( x[i,j] ) + '\n' )
		self.N += x.shape[1]

		for (n,charge,mz1,mz2) in features:
			self.annout.write( str(t) + '\t' + str(mz1) + '\t' + str(mz2) + '\t' + str(charge) + '\n' )

	
		n = 0
		for cid in cids:
			pc = peakclusters[n] # this has both observed peaks and theoretical peaks	
			if pc.idnum >= 0: # it's an observed peak
				npeaks = len( pc.peaks )
				pc.peaks.sort()
				intensity = 0
				for i in xrange(0,npeaks):
					if i == (npeaks-1) or pc.peaks[i][0] != pc.peaks[i+1][0]:
						self.yout.write( str( self.M + cid ) + '\t' + str( pc.peaks[i][1] ) + '\n' )
						if pc.peaks[i][0] < t:
							pc.p()
							print pc.peaks
						intensity = 0		
			n += 1
	
		self.M += x.shape[0]	
	
	# copose the xfile and write out its correct dimensions
	def finish( self ):

		print '\n\nFINISHING AND CLOSING!!!'
		print 'M:', self.M, 'N:', self.N, '\n\n'
		

		tempfile = self.suffix + '_tempx.txt'
		xname = self.xout.name
		self.xout.close()
		fout = open( tempfile, 'w' )
		fout.write( str(self.M) + '\t' + str(self.N) + '\n' )
		fin = open(xname,'r')
		for line in fin:
			fout.write( line )
		fout.close()
		fin.close()
		os.system( 'mv ' + tempfile + ' ' + xname ) 

		tempfile = self.suffix + '_tempy.txt'
		yname = self.yout.name
		self.yout.close()
		fout = open( tempfile, 'w' )
		fout.write( str(self.M) + '\t1\tS\n' ) # The S means that the y matrix is in sparse format!
		fin = open(yname,'r')
		for line in fin:
			fout.write( line )
		fout.close()
		fin.close()
		os.system( 'mv ' + tempfile + ' ' + yname ) 

		self.annout.close()


ms1filename = sys.argv[1]


# Each scan is a tuple (retention time, [ list of tuples (mz,intensity) ] )
# Read in the ms1 file into a list of peaks
if True:
	print 'Reading the ms1file'
	ms1file = open(ms1filename,'r')
	scannum = -1
	peaklist = []
	rtime = 0	
	numpeaks = 0
	for line in ms1file:
		if ((line[0] == "H") or
			(line[0] == "Z") or
			(line[0] == "D")):
			continue
		if (line[0] == "I"):
			if 'RTime' in line:
				rtime = float(line.rstrip().split()[-1])
				peaklist.append( (rtime,[]) )
			continue

		if line[0] == 'S':
			if scannum % 200 == 0:
				print 'Reading scan ' + str(scannum)
			scannum += 1
			continue
		if len(line) > 1:
			tokens = line.rstrip().split()
			mz = float(tokens[0])
			intensity = float(tokens[1])
			#intensity = int(tokens[1])
			numpeaks += 1
			peaklist[-1][1].append( (mz,intensity) )
	ms1file.close() 
	print 'Numpeaks: ' + str(numpeaks)

if True:
	tol = 0.005
	print 'Deleting unduplicated peaks'
	justpeaks = [ peaklist[t][1] for t in range(len(peaklist)) ]	
	justpeaks =  delete_unduplicated_peaklist( justpeaks, tol )

# Here, see how many clusters you get if you combine consecutive scans into blocks rather
# and cluster the combined peaks rather than cluster the peaks scan by scan
if True:
	outdir = ms1filename.replace('.ms1','_ms1_sparsebinned_singlescans_tannotated/')
	if not os.path.exists(outdir):
		os.system( 'mkdir ' + outdir )

	
spans = []
all_binbounds = {}
ratios = []
mspans = []
maxfeatures = 15000
maxbins = maxfeatures
t = 0

toprint = False
while t < len(justpeaks):

	numfeats = 0
	
	xs = []
	ys = []

	ov = output_vector( outdir+str(t)  )	
	while numfeats < maxfeatures and ov.M < maxbins:

		if t >= len(justpeaks):
			break
		print 'Numfeats:', numfeats


		blockpeaks = [ (justpeaks[t][k][0],justpeaks[t][k][1],t) for k in xrange(0,len(justpeaks[t])) ]
		if len(blockpeaks) == 0:
			t += 1
			continue
	
		tol = 0.01
		(clusts, cixes ) = peaklist_to_clusters( blockpeaks, mass_accuracy)

		if toprint:
			print 'num clusters: ' + str(len(clusts)) + ' ' + str( float(len(clusts))/len(blockpeaks) )

		features = [] # each element is a tuple. The tuple is (index of peakcluster correspond to the monoisotope, charge)
		N = len(clusts) 

		tol = 0.01
		if toprint:	
			print 'Hypothesizing isotope features'
		# For each one of these clusters, generate a monoisotopic precursor ion of each charge
		# Only create the feature if both the monoisotopic and the 2nd isotopic peak is is present
		#for n in xrange(10):
		for n in xrange(N):
			#print 'cluster ' + str(n)
			p = clusts[n]
			#p.p()		
			n2 = n+1

			for c in xrange(maxcharge,0,-1):
				# There must be a peakcluster that matches isop
				# where isop is the second isotopic peak
				isop = peakcluster(p.minmz + 1.003355/c, p.maxmz + 1.003355/c, 0, 0)
				while n2 < N:
					compare = compare_clusters(clusts[n2], isop, mass_accuracy)
					#print '-n2: ' + str(n2) + ' comp: ' + str(compare)
					if compare > 0:
						n2 += 1
					if compare < 0:	
						break
					# There is a match! Make this feature and add it to a list somehow
					if compare == 0:
						features.append((n, c, p.minmz * (1 - mass_accuracy), p.maxmz * (1 + mass_accuracy)))
						#print ('-added ' + str(c) )
						break

		if toprint:	
			print str(len(features)) + ' features!' 


		isopeaks = run_computems1_on_features( features, outdir )

		justobsclusts = len(clusts)
		if toprint:	
			print 'Before:' + str(len(clusts))
		peakclusters = clusts
		fn = 0
		for (n,charge,mz1,mz2) in features:
			mz1s = [ mz1+i*1.003355/charge for i in range(0,7) ]
			mz2s = [ mz2+i*1.003355/charge for i in range(0,7) ]

			ipeaks = isopeaks[fn]
			#for i in range(0,len(mz1s)):
			for i in range(0,min(len(ipeaks),len(mz1s))):
				#def __init__( self, minmz, maxmz, startT, endT, intensity ):
				#peakclusters.append( peakcluster(mz1s[i], mz2s[i],t, 0 ) )  
				peakclusters.append( peakcluster(mz1s[i], mz2s[i],t, ipeaks[i][1] ) )  
				peakclusters[-1].idnum = get_ix(fn)
			fn += 1
		if toprint:	
			print 'After:' + str(len(clusts))

		(binbounds, cids) = assign_bins(peakclusters)
		spans = [ bb[1] - bb[0] for bb in binbounds ]
		if toprint:	
			print 'Lenbins: ' + str(len(binbounds)) + ' # observed peak clusters: ' + str(justobsclusts)
			ratio = float(len(binbounds))/justobsclusts 
			ratios.append(ratio)
			print 'Ratio of bins to observed clusters: ' + str( ratio )
			print 'Max span: ' + str( max(spans) )
			mspans.append( max(spans) )	

		# Here, construct x into rowmajor matrices	
		N = len(features)
		numfeats += N
		M = len(binbounds)
		
		# If you're making the x matrix:
		if True:
			x = sp.lil_matrix( (M,N) )
			n = 0
			ts = []
			if toprint:
				print 'Filling the matrix'
			for cid in cids:
				pc = peakclusters[n] # this has both observed peaks and theoretical peaks	
				if pc.idnum < 0: # it's a theoretical peak
					colnum = get_n( pc.idnum )
					x[cid,colnum] = pc.peaks[0][1]
					#print str((cid,colnum)) + ': ' + str(pc.peaks[0][1])
				n += 1
		
			x = x.todense()
		
			ns = range(0,x.shape[1])	 
			footprints = [ get_footprint( x, n ) for n in ns ]
			ns.sort( key = lambda n: footprints[n] )
			
			if len(ns) <= 0:
				t += 1
				continue	
			to_preserve = [ns[0]]
			for i in range(1,len(ns)):
				if footprints[ns[i]] != footprints[ns[i-1]]:
					to_preserve.append(ns[i])

			x = x[:,to_preserve]
	
			ov.output_xy( x, [features[n] for n in to_preserve], cids, binbounds,t )
		t += 1
	ov.finish()		

print t

