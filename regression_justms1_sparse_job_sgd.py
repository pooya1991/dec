import os
import sys

ydir = sys.argv[1]
xdir = ydir
outdir = sys.argv[2]

xfiles = os.listdir(xdir)
x1files = [ xdir+f for f in xfiles if 'xdecoy.txt' in f] # yes decoy 
x1files.sort( key = lambda f: float( f.split('/')[-1].split('_')[0] ) )
y1files = [ f.replace('xdecoy.txt','y.txt') for f in x1files ]

if not os.path.exists(outdir):
	os.system( 'mkdir ' + outdir )

# base = '../scripts/ms1_ols.sh '

numjobs = 0

od = outdir 
if not os.path.exists(od):
	os.system( 'mkdir ' + od )
for (y1f,x1f) in zip(y1files,x1files):
	if not os.path.exists(y1f):
		print y1f + ' does not exist'
		continue 
	t = x1f.split('/')[-1].split('_')[0]
	outputdir = od + t + '/'
	print t
	if not os.path.exists( outputdir ):
		os.system('mkdir ' + outputdir )

	# command = base + x1f + ' ' + y1f + ' ' + outputdir
	command = "./ratest_batch -p " + x1f + " -q " + y1f + " -a both -o " + outputdir + " -i 500 -s 0 -u -z 0.01"
	# command = command.replace('OUTPUTFILE', outputdir+'output.txt' )
	#command += ' > ' + outputdir+'output.txt'
	print(command)
	os.system( command )
	numjobs += 1

print numjobs


