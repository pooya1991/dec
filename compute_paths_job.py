import os
import sys
from compute_paths import siren_lasso
base = 'python compute_paths.py ' 

d = sys.argv[1]
ydir = d
xdir = d
outdir = sys.argv[2]

if not os.path.exists(outdir):
	os.system( 'mkdir ' + outdir )

n = 0
for yfile in os.listdir( ydir ):
	if not '_y.txt' in yfile:
		continue
	#if not '126_' in yfile:
	#	continue
	t = yfile.split('_')[0]
	xfile = xdir + yfile.replace('y.txt','xdecoy.txt')
	yfile = ydir + yfile
	#print t
	#print yfile
	#print xfile
	od = outdir + str(t) + '/'
	if not os.path.exists(od):
		os.system( 'mkdir ' + od )

	# command = base.replace('OUTPUTFILE',od+'log.txt').replace('MEM','8')
	# command += yfile + ' ' + xfile + ' ' + od
	siren_lasso(yfile, xfile, od)
	# print command
	# os.system(command)
	n += 1
	
	#y1file=$1
	#x1file=$2
	#start_t=$3
	#output_dir=$4  
#print str(n) + ' jobs submitted.'
