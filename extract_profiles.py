import os
import sys

execfile('read_rowmajor_matrix.py')
execfile('sg_filter.py')
execfile('smooth.py')

bfile = sys.argv[1]
afile = sys.argv[2]
ms1file = sys.argv[3]
snumfile = ms1file.replace('.ms1', '_scannums.txt')
# Just get the S lines to get the scan numbers for the ms1 file
os.system('grep ^S ' + ms1file + ' > ' + snumfile)

scannums = []
fin = open(snumfile)
for line in fin:
    scannums.append(line.rstrip().split()[-1])
fin.close()


def get_time_bounds(b, t):
    T = b.shape[1]
    start = t - 1
    end = min(T, t + 1)
    while start >= 1:
        if b[0, start - 1] == 0:
            break
        if b[0, start - 1] > b[0, start]:
            break
        start -= 1

    while end < (T - 1):
        if b[0, end + 1] == 0:
            break
        if b[0, end + 1] > b[0, end]:
            break
        end += 1
    # Start is inclusive, end is inclusive
    # returns a tuple where it is (inclusive,exclusive)
    return (start, end + 1)


if True:
    print 'Reading ' + bfile
    b = read_rowmajor_matrix_sparse(bfile)
    b = sp.csr_matrix(b)
    anns = read_tab_delimited(afile)

# ID | precursor number | monoiostopic m/z | charge | peak intensity | elution profile

outfile = ms1file.replace('.ms1', '_precursorprofiles.txt')

fout = open(outfile, 'w')
fout.write('ID\tprecursor number\tmonoisotopic m/z\tcharge\tpeak alpha\telution profile\n')
i = 0
P = len(anns['charge'])
T = b.shape[0]
for p in xrange(P):
    if p % 5000 == 0:
        print p
    z = anns['charge'][p]
    mz = 0.5 * (anns['start m/z'][p] + anns['end m/z'][p])

    alphas = anns['target alphas'][p]
    peak_times = []
    peak_alphas = []
    if alphas == '':
        continue
    elif ',' in alphas:
        alphas = alphas.split(',')
    else:
        alphas = [alphas]

    for a in alphas:
        tokens = a.split(':')
        pt = int(tokens[1])
        alpha = float(tokens[0])
        peak_times.append(pt)
        peak_alphas.append(alpha)
    # subb = b[p,:].todense()
    subb = b[:, p].todense().transpose()
    for (pt, a) in zip(peak_times, peak_alphas):
        (start, end) = get_time_bounds(subb, pt)
        if end - start < 3:
            continue
        fout.write(str(i) + '\t' + str(p) + '\t' + str(mz) + '\t' + str(z) + '\t' + str(a) + '\t')
        i += 1
        fout.write(str(scannums[start]) + ':' + str(subb[0, start]))
        for t in range(start + 1, min(T, end)):
            fout.write(',' + str(scannums[t]) + ':' + str(subb[0, t]))
        fout.write('\n')
fout.close()



