import os
import sys
from scipy import sparse as sp

execfile('read_rowmajor_matrix.py')

if len(sys.argv) == 4:
	d = sys.argv[1] # Directory containing the annotation files
	alignmentfile = sys.argv[2] # File containing the joint alignment
	bfiledir = sys.argv[3] # Directory containing all of the learned Bs. The B file extracted is the /bB.txt file
else:

	sys.stdout.err('Need three arguments: 1. Directory containing annotation files, 2. Alignment file, 3. Directory containing learned abudnances per spectrum.')

annfile = d + 'INITIALT_annfile.txt'
	
if True:
	print 'Reading alignment file'
	alignments = read_tab_delimited(alignmentfile)

	# Go through the alignments and for each scan number, extract the precursor indices they correspond to
	mint = min( alignments['start file'] )
	endt = max( alignments['end file'] )

if True:
	indices_by_time = { t:[] for t in range(mint,endt+1) } # The key is the scans' start time, the element is a list of tuples ( combined_index, sub_index )
	
	P = len( alignments['end file'] )
	for p in xrange(0,P):
		indices = str(alignments['indices'][p])
		if ',' in indices:
			indices = alignments['indices'][p].split(',')
		else:
			indices = [ indices ]
		t = alignments['start file'][p]
		for i in range(0,len(indices)):
			if not '-' in indices[i]:	
				indices_by_time[t].append( (p,int(indices[i])) )
			else:
				ixes = indices[i].split('-')
				for ix in ixes:
					indices_by_time[t].append( (p,int(ix)) )
			t += 1
		
T = endt + 1 
N = len( alignments['end file'] )

print 'Outputting values'

fs = [ f for f in os.listdir(bfiledir) if os.path.exists( bfiledir+f+'/bB.txt') ]
fs.sort( key = lambda f: float(f) )

#fs = fs[:10]
notexists = []
#b = sp.dok_matrix( (N,T) )
boutdir = bfiledir+'combined_B_sp0.txt'

print boutdir
fout = open( boutdir,'w')
fout.write( str(T) + '\t' + str(N*2) + '\n' )
bypassed = 0
# Combined T should be output in as TxN
for i in range(0,len(fs)):
	bp = 0
	f = fs[i]
	f =  bfiledir + f + '/bB.txt'
	if not os.path.exists(f):
		continue
	ti = int(f.split('/')[-2])
	if i == len(fs)-1:
		tf = int(T)
	else:
		tf = int( fs[i+1] )

	anns = read_tab_delimited( annfile.replace('INITIALT',str(ti)) )
 
	#if True: #ti % 100 == 0:
	#	print f
	#print 'Reading'
	b = read_rowmajor_matrix( f ) # When first read in, it's N x 1
	b = b.flatten() # Now it's just N

	N = b.shape[0]
	NT = N/2# number of targets
	bt = b[:NT] # Target matrix
	bd = b[NT:] # Decoy matrix

	# Build a dictionary that gives you the combined index (ci) given the specific index (si)
	ci_by_si = {}
	for t in range(ti,tf):
		if t in indices_by_time:
			ci_by_si[t] = { si:ci for (ci,si) in indices_by_time[t] }
	
	prevt = -1
	si = 0
	decoylines = []
	for (realt,bi) in zip( anns['t'], range(bt.shape[0])):
		#print realt,n
		if realt != prevt:
			for dl in decoylines:
				fout.write(dl)
			decoylines = []
			fout.write( '>\t' + str(realt) + '\n' )
			si = 0
			prevt = realt
		if not si in ci_by_si[realt]:
			bypassed += 1
			bp += 1
			si += 1
			continue
		intensity = bt[bi]
		dintensity = bd[bi]
		if intensity > 0:
			fout.write( str(ci_by_si[realt][si]) + '\t' + str(intensity) + '\n' )	
		if dintensity > 0:
			decoylines.append( str(ci_by_si[realt][si]+P) + '\t' + str(dintensity) + '\n' )		
		si += 1

	for dl in decoylines:
		fout.write(dl)
	#print 'Bypassed ratio:', float(bp)/b.shape[0]
	#break
print 'Num bypassed:', bypassed

fout.close()
