## Concise Siren pipeline

This repo is a concise version of Siren. Specifically, all the files in the C++ folder, the java folder and the scripts folder are thrown away. Also some files in the python folder that were not used in the pipeline are thrown away as well. Instead some files are added to this repo that are not present in the original Siren repo. Specifically there are some R scripts, a sample data and the results that siren generates for that data and some compiled C++ codes. Below is the list of all files and their description:
* **15b.ms1:** the sample data
* **ComputeMS1:** a Linux executable file which is obtained by compiling the contents of Mercury8 folder.
* **ISOTOPE.DAT:** this file is used by ComputeMS1 internally
* **ratest_batch:** this is a Linux executable file which is obtained by compiling the contents of the C++ folder which were related to regression. This file is used for the OLS regression step in Siren
* The python scripts are the same as python scripts in the original Siren repo. There are two main differences: 1) in comupute_paths_job.py which does the Lasso step of Siren, compute_paths.py is evoked as a system call, but I modified compute_paths_job.py to import the functions in compute_paths.py. 2) In regression_justms1_sparse_job_sgd.py, instead of using a shell command to do the regression, I directly called ratest_batch with the arguments that were stored in the shell script
* There are R scripts in this repo which have python counterparts. They do the same thing that the corresponding python script does. Since these files are not well tested, it's better to work with the python scripts instead. There's only one important feature that one has to be aware of. **analysis.R** script does the step of Siren concerning generating X and Y matrices. At first I intended to write my R version of Siren to make it more functional. But later I abandoned the idea due to lack of time.
* **matrix_operations.R:** include functions for working with sparse matrices defined in Siren
* **regression_analysis.R:** this scripts intends to substitute the OLS step of siren with Lasso regression per scan
* **15b_ms1_sparsebinned_singlescans_tannotated:** this folder contains the results concerning making X and Y matrices, decoy X matrices, annotation files and joint annotation file
* **15b_peptide_abundances:** this folder contains the results regarding combining OLS regression coefficients, building a combined B matrix, smoothing and Lasso regression
* **15b_precursorprofiles.txt:** this is the final elution profiles extracted by Siren
* **15b_scannums.txt:** real scan numbers

### Detailed description of pipeline execution

1. Install python2, numpy, scipy, sklearn
    * "sudo apt update"
    * "sudo apt upgrade"
    * "sudo apt install python2.7 python-pip"
    * "pip install numpy"
    * "pip install scipy"
    * "pip install sklearn"

2. clone the git repo

3. In the python session execute the following lines one by one (these are contents of siren.py)
    * import os
    * import sys
    * ms1file = "name of .ms1 file"
    * os.system('python generate_theoretical_features_from_ms1_combine_individual.py ' + ms1file )
    * ms1dir = ms1file.replace('.ms1','_ms1_sparsebinned_singlescans_tannotated/')
    * os.system('python make_siren_decoys.py ' + ms1dir)
    * os.system('python align_ms1_annotations_individual.py ' + ms1dir)
    * bdir = ms1dir.replace('_ms1_sparsebinned_singlescans_tannotated/','_peptide_abundances/')
    * In order to use original Siren OLS: os.system('python regression_justms1_sparse_job_sgd.py ' + ms1dir + ' ' + bdir )
    * In order to use R and glmnet OLS: os.system('./regression_justms1_sparse_job_sgd.R ' + ms1dir + ' ' + bdir)
    * alignmentfile = ms1dir + 'joint_annfile_minlen4.txt'
    * os.system('python combine_ms1_matrices_from_alignment_individual.py ' + ms1dir + ' ' + alignmentfile + ' ' + bdir)
    * bfile = bdir+'combined_B_sp0.txt'
    * os.system('python smooth_extract_withintensities_ms1matrices.py ' + bfile + ' 4 ' + alignmentfile) 
    * intensitiesfile =  bfile.replace('.txt','_intensities_annfile.txt')
    * os.system('python compute_paths_job.py ' + ms1dir + ' ' + bdir )
    * os.system('python add_alphas_withintensities_to_annotations.py ' + ms1dir + ' ' + intensitiesfile + ' ' + bdir) 
    * alphaannfile = intensitiesfile.replace('.txt','_alphas_withintensities.txt')
    * os.system('python extract_profiles_for_diaumpire.py ' + bfile + ' ' + alphaannfile + ' ' + ms1file)