import os
import sys

execfile('read_rowmajor_matrix.py')

# Monoisotopic masses
masses = { 'C':12.0000000, 'H':1.0078246, 'N':14.0030732,'O':15.9949141, 'S':27.976927 } # Masses of elements
averagine = { 'C':4.9384, 'H':7.7583, 'N':1.3577, 'O':1.4773, 'S':0.0417 } # Elemental composition of the average amino acid resideu
amass = sum( [ averagine[e]*masses[e] for e in averagine ] ) # mass of monoisotopic averagine residue

formulafile = 'test_formulas.txt' # Where the averagine formulas will be outputted for each model
specfile = 'test_spectra.txt'

# Avgmz includes the charge-bearing hydrogens
def get_elemental_composition( avgmz, charge ):
	residuemass = (avgmz-masses['H'])*charge
	aunits = float(residuemass)/amass # Compute the number of averagine residues are in this precursor mass
	acomp = { a:int(round(averagine[a]*aunits)) for a in averagine } # Compute the elemental composition
	averagine_mass = sum( [acomp[a]*masses[a] for a in acomp ] ) # Compute the uncorrected masses
	numH = int((residuemass-averagine_mass)/masses['H'])
	acomp['H'] += numH#+charge
	averagine_mass += (numH)*masses['H'] # Add back the hydrogens
	return (acomp,averagine_mass)

# Peaks is a 2d list
# each element is a list of peaks
def normalize_peaks( peaks ):
	for n in range(0,len(peaks)):
		peaklist = peaks[n]
		scalar = 1.0/np.sqrt(sum([ p[1]*p[1] for p in peaklist]))		
		for p in range(0,len(peaklist)):
			peaklist[p] = (peaklist[p][0],peaklist[p][1]*scalar)

# Generate the elemental formulas for the precursor ions for which
# we want to generate isotope distributions. Then call the computems1
# binary to generate these isotope distributions
def run_computems1_on_features( featureinfo, d = '' ):

	to_print = False

	ffile = d + formulafile
	sfile = d + specfile

	command = './ComputeMS1 ' + ffile + ' > ' + sfile 
	if to_print:
		print 'Len featureinfo: ' + str(len(featureinfo))
	fout = open(ffile,'w')
	submitted = [] # A list of all the mz/charge pairs for which to compute an isotope distribution
	formulanum = -1
	duplicates = 0

	formulas = set([])
	#for (mz1,mz2,charge) in featureinfo:
	for (n,charge,mz1,mz2) in featureinfo: # This is if you give it features, not featureinfo
		smz = (mz1+mz2)/2.0
		#print smz
		(comp,mass) = get_elemental_composition( smz, charge )
		compstring = ''
		for c in ['C','H','N','O','S']:
			compstring += c+str(comp[c])
		compstring += '\t'+str(charge)
		if not compstring in formulas: # Sometimes multiple mz/charge pairs will be rounded to the same elemental formula. Make sure you only output each formula once.				
			fout.write( compstring + '\n' )
			formulanum += 1
		else:
			duplicates += 1
		submitted.append( (smz,charge,formulanum) ) 
	fout.close()

	if to_print:
		print command
	formulanum += 1
	if to_print:
		print 'Submitted', formulanum
		print 'Duplicates', duplicates
		print 'Running computems1 command to ' + str(sfile)
		print command
	s = os.system( command )

	# Read the output of compute ms1 and store them into the list of lists "peaks"
	# Where each element in peaks is a list of tuples
	fin = open(sfile,'r')
	charges = []
	peaks = [[] for i in xrange(0,formulanum)]
	n = -1
	isvals = False
	for line in fin:
		if 'successful' in line:
			isvals = False
			continue
	
		if 'Sequence' in line:
			n += 1

		if 'Average Integer' in line:
			tokens = line.rstrip().split(',')
			mass = float(tokens[0].split()[-1])
			mz = float(tokens[-1].split()[-1])
			charge = round(mass/mz)
			charges.append(charge)

		if 'Calculation' in line:
			isvals = True
			continue

		if isvals:
			tokens = line.rstrip().split('\t')
			mz = float(tokens[0])
			#if mzbin >= M:
			#	continue 
			intensity = float(tokens[1])
			if intensity < 0.1:
				continue
			#peps[-1].append( mz )
			peaks[n].append( (mz,intensity))
	fin.close()
	normalize_peaks(peaks) # Normalize each list of peaks such that the magnitude is 1
	if to_print:
		print 'Done: ' + str(s)
		print 'Specfile: ' + sfile + ' lenpeaks: ' + str(len(peaks))
	return(peaks)

