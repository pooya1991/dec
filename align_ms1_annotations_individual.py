import os
import sys

execfile('read_rowmajor_matrix.py')
execfile('peakclusterforbinning.py')
def sort_dict_charge_minmz( d ):
	order = range( len(d['charge']) )	
	order.sort( key = lambda i: (d['t'][i],d['charge'][i],d['min m/z'][i]) )
	for k in d:
		d[k] = [ d[k][o] for o in order ]
	return order

def output_precursor( fout, p ):
	fout.write(  str(p.minmz) + '\t' + str(p.maxmz) + '\t' +  str(p.charge) + '\t' + str(p.startscan) + '\t' + str(p.endscan) + '\t' )
	indices = p.indices
	fout.write( str(indices[0]) )
	for i in indices[1:]:
		if i == None:
			print 'NONE!'
			p.p()

		if i != None and i < 0:
			#print p.indices
			fout.write(str(i)) # This means that multiple precursors within that scan contribute to the overall one
		else:
			fout.write(','+str(i))
	fout.write('\n')

d = sys.argv[1]

massHMono = 1.0078246 

annotationfiles = [ d+f for f in os.listdir(d) if 'annfile' in f and not '.swp' in f and not 'joint' in f ]
annotationfiles.sort( key = lambda af: int( af.split('/')[-1].split('_')[0] ) )
print annotationfiles
af = annotationfiles[0]
anns = read_tab_delimited( af )

tanns = segregate_annotations_by_time( anns )
ts = sorted(tanns.keys())
t = ts[0]
print ts
anns = tanns[t]
N = len(anns['charge'])
order = sort_dict_charge_minmz(anns)
#def __init__( self, mz1, mz2, charge, t, index ):
precursors = [ precursor( anns['min m/z'][n], anns['max m/z'][n], anns['charge'][n], t, order[n]) for n in xrange(0,N) ]

numfinished = 0
lens = []
minlen = 4

tindex = 1
afindex = 1

outfile = d + 'joint_annfile_minlen'+str(minlen)+'.txt'
fout = open(outfile,'w')
fout.write( 'start m/z\tend m/z\tcharge\tstart file\tend file\tindices\n')

total = len(precursors)
unaligned = 0
ttol = 1
while tindex < len(ts):
	t = ts[tindex]
	prevt = t-ttol
	anns = tanns[t]
	total += len(tanns[t]['charge'])
	tindex += 1

	if tindex % 50 == 0:
		print tindex
                print 'Num features evaluated:', numfinished
                print 'Features that do not persist over ' + str(minlen) + ' scans'
		print str(unaligned)+'/'+str(total) + ',' + str( float(unaligned)/total )

	# If we need to read more annotation files, let's do it!
	if tindex == len(ts) and afindex < len(annotationfiles):
		af = annotationfiles[afindex]
		afindex += 1
		tanns = segregate_annotations_by_time( read_tab_delimited(af) )
		ts += sorted(tanns.keys())

	precursors.sort( key = lambda p: (p.charge,p.minmz) )
	order = sort_dict_charge_minmz(anns)

	J = len(anns['charge'])
	I = len(precursors)	
	ixes_to_continue = set([])

	i = 0
	j = 0
	while i < I and j < J:
		acharge = anns['charge'][j]
		pci = precursors[i]
		minmz = anns['min m/z'][j]
		maxmz = anns['max m/z'][j]
		# They don't align
		if ( minmz > pci.maxmz and acharge == pci.charge) or acharge > pci.charge:
			if not i in ixes_to_continue: # If precursor[i] is finished!
				#finished.append( pci )
				flen = pci.endscan - pci.startscan
				#if pci.endscan < prevt:
				#if True:
				if flen >= minlen-1:
					if pci.to_output:
						output_precursor( fout, pci ) #
						lens.append(flen)	
					numfinished += 1
				else:
					unaligned += len(pci.indices)
					#if pci.startscan == pci.endscan:
						#unaligned += 1
			#	else:
					# Preserve it if it is within ttol of t, even if it isn't
					# extended in t.
					#ixes_to_continue.add( i )
			i += 1
			continue
		# They don't align
		elif ( maxmz < pci.minmz and acharge == pci.charge) or acharge < pci.charge: # anns[x][j] does not yet exist
			ixes_to_continue.add( len(precursors) )
			precursors.append( precursor( anns['min m/z'][j], anns['max m/z'][j], anns['charge'][j], t, order[j]  ) )
			j += 1
		# They do align!
		# Add it to an old one
		else:
			# Another annotation at time t already aligned to an active precursor.
			if pci.endscan == t:
				oldix = pci.indices[-1]			
				pci.indices.append( -order[j] ) # If it's negative, it means that it corresponds to the same scan as the previous one
			else: # It didn't
#				for temp_t in xrange(pci.endscan+1,t):
#					pci.indices.append(None) # Placeholder if it's nothing
				pci.indices.append( order[j] )
			pci.endscan = t
			pci.minmz = min(pci.minmz,anns['min m/z'][j])   
			pci.maxmz = min(pci.maxmz,anns['max m/z'][j]) # alex
			ixes_to_continue.add( i )

			# See if anns[*][j] connects to subsequent Is. If the subsequent ones do, combine them with this
			# original i.
			k = i+1
			while k < I:
				pck = precursors[k]	
				if pck.minmz <= anns['max m/z'][j] and pck.charge == anns['charge'][j]:
					# Combine pci and precursors[k]!
					#print i, k
					pci.incorporate( pck )	
				k += 1
			j += 1

	# Add the remaining ones!!!
	while j < J: # anns[x][j] does not yet exist		ixes_to_continue.add( len(precursors) )
		ixes_to_continue.add( len(precursors) )
		precursors.append( precursor( anns['min m/z'][j], anns['max m/z'][j], anns['charge'][j], t, j  ) )
		j += 1

	precursors = [ precursors[k] for k in ixes_to_continue ]

print 'Unaligned'

# Output the last precursors
if True:
	numscansperblock = 5
	precursors.sort( key = lambda p: (p.charge,p.minmz) )
	for pci in precursors: # If precursor[i] is finished!
		#finished.append( pci )
		flen = pci.endscan - pci.startscan
		#if numscansperblock >= minlen or flen >= minlen-1:
		if flen >= minlen-1:
			output_precursor( fout, pci ) #
			lens.append(flen)	
			numfinished += 1

fout.close()


