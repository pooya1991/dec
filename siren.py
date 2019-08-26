import os
import sys 

ms1file = sys.argv[1]
mzxmlfile = ''
if len(sys.argv) > 2:
	mzxmlfile = sys.argv[2]

sys.stdout.write('Generating features from ' + ms1file )
os.system('python generate_theoretical_features_from_ms1_combine_individual.py ' + ms1file )
ms1dir = ms1file.replace('.ms1','_ms1_sparsebinned_singlescans_tannotated/')

sys.stdout.write('Generating decoy features (reversing sequences of top isotope intensities)' + ms1file )
os.system('python make_siren_decoys.py ' + ms1dir )

sys.stdout.write('Linking features over elution time')
os.system('python align_ms1_annotations_individual.py ' + ms1dir )
 
sys.stdout.write('Generating un-regularized solutions')
bdir = ms1dir.replace('_ms1_sparsebinned_singlescans_tannotated/','_peptide_abundances/')
os.system('python regression_justms1_sparse_job_sgd.py ' + ms1dir + ' ' + bdir )

sys.stdout.write('Combining solutions into elution profiles')
alignmentfile = ms1dir + 'joint_annfile_minlen4.txt'
os.system('python combine_ms1_matrices_from_alignment_individual.py ' + ms1dir + ' ' + alignmentfile + ' ' + bdir)

sys.stdout.write('Finding elution peaks')
bfile = bdir+'combined_B_sp0.txt'
print 'python smooth_extract_withintensities_ms1matrices.py ' + bfile + ' 4 ' + alignmentfile
os.system('python smooth_extract_withintensities_ms1matrices.py ' + bfile + ' 4 ' + alignmentfile ) 
intensitiesfile =  bfile.replace('.txt','_intensities_annfile.txt')

sys.stdout.write('Generating regularization paths')
os.system('python compute_paths_job.py ' + ms1dir + ' ' + bdir ) 

sys.stdout.write('Outputting regularization alphas for elution peaks')
os.system('python add_alphas_withintensities_to_annotations.py ' + ms1dir + ' ' + intensitiesfile + ' ' + bdir  ) 
alphaannfile = intensitiesfile.replace('.txt','_alphas_withintensities.txt')

sys.stdout.write('Extracting elution profiles for DIA-Umpire')
os.system('python extract_profiles_for_diaumpire.py ' + bfile + ' ' + alphaannfile + ' ' + ms1file )
precursorprofilefile = ms1file.replace('.ms1','_precursorprofiles.txt')

# Run DIA-Umpire on the learned elution profiles
# if an mzxml file is provided
if mzxmlfile != '':

    # Generate diaumpire parameter file with the correct alphaannfile
    fin = open('../java/diaumpire_nobackground_template.se','r')
    newparamfile = intensitiesfile.replace('.txt','_diaumpire_params.txt')
    fout = open( newparamfile, 'w')
    for line in fin:
        if 'OwnPeakURL' == line[:10]:
            fout.write(line.replace('OWNPEAKFILE', precursorprofilefile) )
        else:
            fout.write(line)
    fout.close()

    # Run DIA-Umpire
    sys.stdout.write('Running DIA-Umpire')
    os.system('java -jar -Xmx8G ../java/DIA_Umpire_SE2.java ' + mzxmlfile + ' ' + newparamfile)





