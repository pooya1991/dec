import os
import sys
from sklearn import linear_model
#from matplotlib import pyplot as plt
import scipy
import time

execfile('read_rowmajor_matrix.py')

# The objective functioned learned here is this:
# (1 / (2 * n_samples)) * ||y - Xw||^2_2 + alpha * ||w||_1

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


def siren_lasso(y1file, x1file, outdir):
	# y1file = sys.argv[1]
	# x1file = sys.argv[2]
	# outdir = sys.argv[3]


	if True:

		print 'Learning path for scan ' + str(y1file)
		y = read_rowmajor_matrix( y1file )
		(x,binbounds) = read_rowmajor_matrix_flexbins( x1file )
		#b = read_rowmajor_matrix( bfile )


	binblocks = find_bin_blocks( binbounds )
	for (t,f,l) in binblocks:

		if t == -1:
			continue

		# Subdivide the x and y matrices to correspond to a single scan
		suby = y[f:l,:]
		subx = x[f:l,:]
		nzs = subx.nonzero()[1]
		x1 = nzs.min()
		x2 = nzs.max()+1
		scannum = t

		outfile = outdir + 'path_' + str(scannum) + '_' + str(x1) + 'to' + str(x2) + '.txt'

		if os.path.exists( outfile ):
			continue
		print 'Start time:', time.time()
		print 'Learning path for scan ' + str(scannum)
		sys.stdout.flush()
		subx = subx[:,x1:x2]
		lars = linear_model.LassoLars(fit_intercept=False, positive=True,max_iter=np.inf,alpha=0)
		lars.fit(subx,suby)

		alphaout = open( outdir + 'alphas_' + str(scannum) + '.txt', 'w' )
		for alpha in lars.alphas_:
			alphaout.write( str(alpha) + '\n' )
		alphaout.close()

		output_rowmajor_matrix( lars.coef_path_, outdir + 'path_' + str(scannum) + '_' + str(x1) + 'to' + str(x2) + '.txt')
		print 'End time:', time.time()

