import os
import sys
import numpy as np
import scipy.sparse as sp

def write_tab_delimited( d, outfile, keys=[] ):

	if keys == []:
		keys = sorted( d.keys() )
	fout = open(outfile,'w')
	fout.write( keys[0] )
	for k in keys[1:]:
		fout.write( '\t' + k )
	fout.write('\n')

	N = len(d[k])
	for n in range(N):
		fout.write( str( d[keys[0]][n] ) )
		for k in keys[1:]:
			fout.write( '\t' + str(d[k][n]) )
		fout.write('\n')
	fout.close()


def read_tab_delimited( f, sep="\t", header=True ): # returns a dictionary
	fin = open(f,'r')
	if header:	
		header = fin.readline()
		if header == '':
			return {}
		if header[-1] == '\n':
			header = header[:-1]
		header = header.split(sep)
		data = { h:[] for h in header }
		lh = len(header)
		for line in fin:
			if line[-1] == '\n':
				line = line[:-1]
			tokens = line.split(sep)
			if lh != len(tokens):
				print 'Mismatch!', lh, len(tokens)
				print line
			for i in xrange(0,len(tokens)):
				try:
					t = float(tokens[i])
					if int(t) == t:
						t = int(t)
					data[header[i]].append(t)
				except ValueError:
					data[header[i]].append(tokens[i])
				except OverflowError:
					data[header[i]].append(tokens[i])

	else:
		data = [ ]	
		for line in fin:
			if line == '':
				return {}
			tokens = line.rstrip().split()
			for i in range(len(tokens)):
				if len(data) <= i:
					data.append([])
				try: 
					t = float(tokens[i])
					if int(t) == t:
						t = int(t)
		                        data[i].append(t)
				except ValueError:
		                        data[i].append(tokens[i])				
				except OverflowError:
		                        data[i].append(tokens[i])				
		fin.close()

	fin.close()
	return data
		
def print_row( d, r, keys = [], column = True ):

	if column:
		if keys == []:
			keys = sorted(d.keys())
		sys.stdout.write( keys[0] + ': ' + str(d[keys[0]][r]) + '\n')
		for k in keys[1:]:
			sys.stdout.write( k + ': ' + str(d[k][r]) + '\n' )
	
	else:
		if keys == []:
			keys = sorted(d.keys())
		sys.stdout.write( keys[0] )
		for k in keys[1:]:
			sys.stdout.write( '\t' + k )
		sys.stdout.write( '\n'+str(d[keys[0]][r]) )
	
		for k in keys[1:]:
			sys.stdout.write( '\t'+str(d[k][r]) )
		sys.stdout.write('\n')

# This is intended to just read the last column of a
# matrix outputted from LASSO Lars of the regularization path
def read_rowmajor_matrix_lastcol( filename, mtype='' ):
	fin = open(filename,'r')
	dim = fin.readline().rstrip().split()
	dim = (int(dim[0]),int(dim[1]))

	nrow = dim[0]
	ncol = dim[1]
	lastcol = ncol-1
	mat = np.zeros( nrow )

	row = -1
	for line in fin:
		if '#' == line[0]:
			continue
		tokens = line.rstrip().split()
		if '>' == tokens[0]:
			row = int(tokens[1])
			continue
		col = int(tokens[0])
		if col == lastcol:
			mat[row,col] += float(tokens[1])
	fin.close()

	if 'csr' in str(mtype):
		mat = sp.csr_matrix(mat)
	elif 'csc' in str(mtype):
		mat = sp.csc_matrix(mat)
	return mat

# First line is nrows ncols
# Row index is given in line that starts with >
# Column index and element value follows
def read_rowmajor_matrix( filename, mtype='' ):
	fin = open(filename,'r')
	dim = fin.readline().rstrip().split()
	colvector = False
	if len(dim) > 2 and dim[2] == 'S':
		 colvector = True
	dim = (int(dim[0]),int(dim[1]))

	nrow = dim[0]
	ncol = dim[1]
	if mtype != '':
		mat = sp.lil_matrix( dim )
	else:
		mat = np.zeros( dim )
	if colvector:
		for line in fin:
			tokens = line.rstrip().split()
			col = 0
			row = int(tokens[0])
			mat[row,col] += float(tokens[1])

	else:
		row = -1
		col = -1
		for line in fin:
			if '#' == line[0]:
				continue
			tokens = line.rstrip().split()
			if '>' == tokens[0]:
				row = int(tokens[1])
				continue
			col = int(tokens[0])
			if col >= mat.shape[1]:
				continue
			if row >= nrow or col >= ncol:
				continue
			mat[row,col] += float(tokens[1])
		fin.close()



	if 'csr' in str(mtype):
		mat = sp.csr_matrix(mat)
	elif 'csc' in str(mtype):
		mat = sp.csc_matrix(mat)
	return mat


def get_binbounds( filename ):
	fin = open(filename,'r')
	line = fin.readline()
	dim = line.rstrip().split()
	dim = (int(dim[0]),int(dim[1]))
	nrow = dim[0]
	binbounds = [ (0,0,-1) ]*nrow
	row = -1
	col = -1
	for line in fin:
		if '>' == line[0]:
			tokens = line.rstrip().split()
			row = int(tokens[1])
			if len(tokens) == 5:
				binbounds[row]= (float(tokens[2]),float(tokens[3]), int(tokens[4])) 
			else:
				binbounds[row]= (float(tokens[2]),float(tokens[3])) 
	fin.close()
	return binbounds



def read_rowmajor_matrix_flexbins( filename, mtype='' ):
	fin = open(filename,'r')

	line = fin.readline()

	dim = line.rstrip().split()
	dim = (int(dim[0]),int(dim[1]))
	nrow = dim[0]
	ncol = dim[1]
	binbounds = [ (0,0,-1) ]*nrow
	if mtype != '':
		mat = sp.dok_matrix( dim )
	else:
		mat = np.zeros( dim )
	row = -1
	col = -1
	for line in fin:
		if '#' == line[0]:
			continue
		tokens = line.rstrip().split()
		if '>' == tokens[0]:
			row = int(tokens[1])
			if len(tokens) == 4:
				binbounds[row]= (float(tokens[2]),float(tokens[3])) 
			elif len(tokens) > 4:
				binbounds[row]= (float(tokens[2]),float(tokens[3]),int(tokens[4])) 
			continue
		col = int(tokens[0])
		if col >= mat.shape[1]:
			continue
		if row >= nrow or col >= ncol:
			continue
		mat[row,col] += float(tokens[1])
	fin.close()
	if 'csr' in str(mtype):
		mat = sp.csr_matrix(mat)
	elif 'csc' in str(mtype):
		mat = sp.csc_matrix(mat)
	return (mat,binbounds)

def discretize(mz,bw,offset):
	return int((mz/bw) + 1.0 - offset)

# Read peaklist
# Into a sparse, discretized matrix
def read_peaklist( filename, res=1.0005079 ):
	fin = open(filename,'r')
	lines = fin.readlines()
	fin.close()

	M = discretize( float(lines[-1].rstrip().split()[0]), res, 0.0 )+1

	dim = (M,1)
	mat = np.zeros( dim )
	row = -1
	col = -1
	for line in lines:
		tokens = line.rstrip().split()
		r = discretize( float(tokens[0]), res, 0.0 )
		mat[r,0] = float(tokens[1])
	return mat

def read_rowmajor_matrix_coo( filename ):
	fin = open(filename,'r')
	dim = fin.readline().rstrip().split()
	dim = (int(dim[0]),int(dim[1]))
	row = -1
	col = -1

	rows = []
	cols = []
	vals = []

	for line in fin:
		if '#' == line[0]:
			continue
		tokens = line.rstrip().split()
		if '>' in line:
			row = int(tokens[1])
			continue
		col = int(float(tokens[0]))
		if row >= dim[0] or col >= dim[1]:
			print line
		else:
			rows.append(row)
			cols.append(col)
			vals.append( float(tokens[1]) )

	fin.close()
	mat = sp.coo_matrix((vals, (rows, cols)), shape=dim)
	return mat

def read_rowmajor_matrix_sparse( filename ):
	fin = open(filename,'r')
	dim = fin.readline().rstrip().split()
	dim = (int(dim[0]),int(dim[1]))
	mat = sp.lil_matrix( dim )
	row = -1
	col = -1

	if dim[1] == 1:
		for line in fin:
			if '#' == line[0]:
				continue
			tokens = line.rstrip().split()
			row = int(tokens[0])
			if row >= dim[0]:
				print 'Does not fit: ' + line
			else:
				mat[row,0] += float(tokens[1])
	else:
		for line in fin:
			if '#' == line[0]:
				continue
			tokens = line.rstrip().split()
			if '>' in line:
				row = int(tokens[1])
				continue
			col = int(float(tokens[0]))
			if row >= dim[0] or col >= dim[1]:
				print line
			else:
				mat[row,col] += float(tokens[1])
	fin.close()
	return mat

def read_and_coarsen_rowmajor_matrix_sparse( filename, factor ):
	fin = open(filename,'r')
	dim = fin.readline().rstrip().split()
	dim = (int(dim[0])/factor,int(dim[1]))
	mat = sp.lil_matrix( dim )
	for line in fin:
		if '#' == line[0]:
			continue
		tokens = line.rstrip().split()
		if '>' in line:
			row = int(tokens[1])
			row = min(row/factor,dim[0]-1)
			continue
		col = int(tokens[0])
		mat[row,col] += float(tokens[1])
	fin.close()
	return mat

def output_lilmatrix( m, filename ):
        fout = open( filename, 'w' )
        fout.write( str(m.shape[0]) + '\t' + str(m.shape[1]) + '\n' )
        for r in xrange(0,m.shape[0]):
                nonzerocolumns = m[r,:].nonzero()[1]
                if nonzerocolumns.shape[0] > 0:
                        fout.write( '> ' + str(r) + '\n' )
                for c in nonzerocolumns.tolist():
                        fout.write( str(c) + '\t' + str(m[r,c]) + '\n' )
        fout.close()

def output_matrix( M, outfile ):
        nrow = M.shape[0]
        ncol = M.shape[1]
        (I,J) = M.nonzero()
        fout = open( outfile, 'w')
        fout.write( str(nrow) + '\t' + str(ncol) + '\n' )
        for n in xrange(0,I.shape[0]):
                i = I[n]
                j = J[n]
                fout.write( str(i) + '\t' + str(j) + '\t' + str( M[i,j] ) + '\n' )
        fout.close()

def output_rowmajor_matrix( M, outfile, binbounds = [] ):
        nrow = M.shape[0]
        ncol = M.shape[1]
        M = sp.csr_matrix(M)
        (I,J) = M.nonzero()
        fout = open( outfile, 'w')
        previ = -1
        fout.write( str(nrow) + '\t' + str(ncol) + '\n' )
        for n in xrange(0,I.shape[0]):
                i = I[n]
                if previ != i:
			if len(binbounds) == 0:
	                        fout.write( '> ' + str(i) + '\n' )
			else:
	                        fout.write( '> ' + str(i) + '\t' + str(binbounds[i][0]) + '\t' + str(binbounds[i][1]) + '\n' )
                previ = i
                j = J[n]
                fout.write( str(j) + '\t' + str( M[i,j] ) + '\n' )
        fout.close()

def combine_dict_of_lists( a, b ):
	c = {}
	lena = 0
	lenb = 0	
	for k in a.keys():
		if not k in b:
			continue
		c[k] = a[k] + b[k]	
		lena = len(a[k])
		lenb = len(b[k])
	return c



