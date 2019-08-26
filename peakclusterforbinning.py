massHMono = 1.0078246

class precursor:
	numps = 0
	def __init__( self, mz1, mz2, charge, t, index ):
		self.minmz = mz1 # over all times
		self.maxmz = mz2 # over all times
		self.charge = charge
		self.idnum = precursor.numps
		self.startscan = t
		self.endscan = t
		self.indices = [index] # Each index is the index of a precursor at the corresponding time, indicated as the range between self.startscan and self.endscan. Negative indices indicate the positive index, but it represents the same scan as the previous scan.

		self.minmass = (mz1-massHMono)*charge 
		self.maxmass = (mz2-massHMono)*charge
		self.to_output = True
		precursor.numps += 1

	def incorporate( self, p ):

		#print 'idnums:', self.idnum, p.idnum

		#print 'Incorporating'
		if self.idnum == p.idnum:
		#	print '...not'
			return

		#print (self.minmz,self.maxmz), self.startscan, self.endscan, self.indices
		#print (p.minmz,p.maxmz), p.startscan, p.endscan, p.indices
	
		# Combine their indicesssss
		# Pad the indices with Nones until they span the same distance
		if self.startscan < p.startscan:
			p.indices = [None]*(p.startscan-self.startscan)+p.indices  			
		elif self.startscan > p.startscan:
			self.indices = [None]*(self.startscan-p.startscan)+self.indices  			
		if self.endscan > p.endscan:
			p.indices = p.indices + [None]*(self.endscan-p.endscan)  			
		elif self.endscan < p.endscan:
			self.indices = self.indices + [None]*(p.endscan-self.endscan)  			

		p.idnum = self.idnum # the incorporated one changes its ID for future reference so it doesn't get incorporated again.
		self.minmz = min( self.minmz, p.minmz )
		self.maxmz = max( self.maxmz, p.maxmz )
		self.startscan = min( self.startscan, p.startscan )
		self.endscan = max( self.endscan, p.endscan )

		if False:
		#if self.idnum == 878:
		#if None in self.indices or None in p.indices:
			print 'Padded indices'
			print self.indices
			print p.indices
			self.p()

		# Combine the padded indices
		newindices = [ [] for i in range(self.endscan-self.startscan+1) ]
		i = 0
		j = 0
		I = len(self.indices)
		J = len(p.indices)			
	
		#print 'Newlen:', len(newindices)	
		for t in range(len(newindices)):
			#print str(t) + ' i,j:', (i,j)
			z = True
			if self.indices[i] != None:
				newindices[t].append( self.indices[i] )
				z = False
				i += 1
			
				while i<I and self.indices[i] != None and self.indices[i] < 0:
					newindices[t].append( self.indices[i] )
					i += 1		
			else:
				i += 1

			if p.indices[j] != None:
				if z:
					newindices[t].append( p.indices[j] )
				else:
					# This means the indices overlap based on scan number
					newindices[t].append( -p.indices[j] )
				j += 1
				while j<J and p.indices[j] != None and p.indices[j] < 0:
					newindices[t].append( p.indices[j] )
					j += 1		
			else:
				j += 1

		self.indices = []
		for t in range(0,len(newindices)):
			if len( newindices[t] ) == 0:
				self.indices.append(None)
			else:
				self.indices += newindices[t]
	
		p.to_output = False
	


	def p(self):
		print 'ID: ' + str(self.idnum) + ', charge: ' +  str(self.charge) + ', m/z: ' + str(self.minmz) + '-' + str(self.maxmz)+ ', scannums: ' + str( (self.startscan,self.endscan) )

# Here, precursors is a list of precursors
#Below is the header for the annotation files
#t       min m/z max m/z charge

def segregate_annotations_by_time( anns ):
	tindices = {} # to contain a list of indices for the row in anns that corresponds to 
	# a particular time
	keys = anns.keys()
	ts = anns['t']
	
	for i in xrange(len(ts)):
		t = ts[i]
		if not t in tindices:
			tindices[t] = [i]
		else:
			tindices[t].append(i)

	times = { t:{} for t in tindices }
	for k in keys:
		for t in times:
			times[t][k] = [ anns[k][i] for i in tindices[t] ]	
	return times


# Here, precursors is a list of precursors
def segregate_dictionary_by_time( anns ):
	times = {}
	keys = precursors.keys()


	for p in precursors:
		if p.startscan in times:
			times[p.startscan].append(p)
		else:
			times[p.startscan] = [p]
	return times
			
# a minmz and a maxmz associated with a set of peaks 
# that are contiguous across scans and connected by a path
# of some m/z tolerance.
# The existence of peakclusters is mostly agnostic to
# time/scannumber though.
class peakcluster:
	numpcs = 0
	#def __init__( self, mz , t,  intensity ):
	#	self.minmz = mz # over all times
	#	self.maxmz = mz # over all times
	#	self.idnum = peakcluster.numpcs
	#	self.peaks = [ (t,intensity) ]
	#	peakcluster.numpcs += 1

	def __init__( self, mz1, mz2, t, intensity ):
		self.minmz = mz1 # over all times
		self.maxmz = mz2 # over all times
		self.idnum = peakcluster.numpcs
		self.peaks = [ (t,intensity) ]
		peakcluster.numpcs += 1

	# This is agnostic to time
	def add( self, mz, tol):
		if (mz+tol) >= self.minmz:
			if (mz-tol) <= self.maxmz:
				self.maxmz = max(maxmz,mz)
				self.minmz = min(minmz,mz)
		#		self.mzlist.add(mz)
				return True
		return False

	def p(self):
		print str(self.idnum) + ' ('+str(self.minmz)+','+str(self.maxmz)+') '

	def pfile(self, fout, minlen = 2 ):
		fout.write(str(self.minmz)+'\t'+str(self.maxmz)+'\n')

# peaklist_to_peakclusters
# Each element in peaklist is a tuple: (mz, intensity, scannum) 
def peaklist_to_clusters( peaklist, mass_accuracy):
	peaklist.sort(key=lambda p: p[0])
	clusters = []
	mz = peaklist[0][0]
	intensity = peaklist[0][1]
	realt = peaklist[0][2]
	clusters.append( peakcluster( mz, mz, realt, intensity ) ) 

	cids = [0]*len(peaklist)
	ix = 1
	for (mz,intensity,realt) in peaklist[1:]:
		if clusters[-1].maxmz >= (mz * (1 - mass_accuracy)):
			clusters[-1].maxmz = mz
			cids[ix] = cids[ix-1]
			clusters[-1].peaks.append((realt, intensity))
		else:
			clusters.append(peakcluster(mz, mz, realt, intensity))	# FIXME
			cids[ix] = cids[ix-1]+1
		ix += 1
	return (clusters, cids)

def compare_clusters( a, b, mass_accuracy):
	if a.maxmz * (1 + mass_accuracy) < b.minmz:
		return 1 # b is greater than a
	if b.maxmz * (1 + mass_accuracy) < a.minmz:
		return -1 # b is less than a
	return 0 # they overlap

# peaklist_to_peakclusters
# try to meld clusters to see if they overlap. They have to overlap, not just
# be within a certain tolerance of each other
def assign_bins( clist ):
	clist.sort(key=lambda c: c.minmz )
	binbounds = []
	N = len(clist)
	cids = [0]*N
	startmz = clist[0].minmz
	endmz = clist[0].maxmz
	for n in xrange(1,N):
		# They're connected!
		if clist[n].minmz < endmz:
			cids[n] = cids[n-1]
			endmz = max(endmz,clist[n].maxmz)
		else:
			# They're not connected!
			binbounds.append( (startmz,endmz) )
			cids[n] = cids[n-1] + 1
			startmz = clist[n].minmz
			endmz = clist[n].maxmz		
	binbounds.append( (startmz,endmz) )	
	return binbounds, cids

# Try deleting peaks that don't exist within a certain tolerance in at least two consecutive scans
def delete_unduplicated_peaks( scanvals, tol=0.005 ): 
	svs = [ set([]) for i in xrange(len(scanvals)) ]
	for t in range(0,len(svs)-1):
		if t % 100 == 0:
			print 't: ' + str(t)
		if t == 0:
			a = sorted(list(scanvals[t]))
		else:
			a = b
		b = sorted(list(scanvals[t+1]))

		i = 0
		j = 0
		I = len(a)
		J = len(b)

		while i<I and j<J:
			diff = a[i]-b[j]
			if abs(diff) <= tol:
				svs[t].add(a[i])	
				svs[t+1].add(b[j])	
			if diff > 0:
				j += 1
			else:
				i += 1
	return svs	 

# Try deleting peaks that don't exist within a certain tolerance in at least two consecutive scans
# Peaklist is a list of lists of tuples, where each tuple is (mz, intensity)
def delete_unduplicated_peaklist( peaklists, tol=0.005 ): 
	svs = [ [] for i in xrange(len(peaklists)) ]
	for t in range(0,len(svs)-1):
		if t % 100 == 1:
			print 't: ' + str(t) + ' ' + str(len(svs[t-1])), len(a), len(b)		
		if t == 0:
			a = peaklists[t]
		else:
			a = b
		b = peaklists[t+1]

		i = 0
		j = 0
		I = len(a)
		J = len(b)

		while i<I and j<J:
			diff = a[i][0]-b[j][0]
			if abs(diff) <= tol:
				svs[t].append(a[i])	
				svs[t+1].append(b[j])	
			if diff > 0:
				j += 1
			else:
				i += 1
                svs[t] = list(set(svs[t]))
                svs[t].sort()
	return svs	


