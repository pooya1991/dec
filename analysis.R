library(stringr)
source("generate_theoretical_features_from_ms1_combine_individual.R")
source("peakclusterforbinning.R")
source("make_averagine_models.R")

# setting some parameters -------------------------------------------------

ms1filename <- "15b.ms1"
outdir <- stringr::str_replace(ms1filename, "\\.ms1$", '_ms1_sparsebinned_singlescans_tannotated/')
formulafile <- "test_formulas.txt"
specfile <- "test_spectra.txt"
maxfeatures <- 15000
maxbins <- maxfeatures
# consider charges 4, 3, 2, and 1
maxcharge <- 4L

# monoisotopic masses -----------------------------------------------------

# Masses of elements
masses <- c(C = 12, H = 1.0078246, N = 14.0030732, O = 15.9949141, S = 31.972072)
# Elemental composition of the average amino acid resideu
averagine <- c(C = 4.9384, H = 7.7583, N = 1.3577, O = 1.4773, S = 0.0417)
# Mass of monoisotopic averagine residue
amass <- sum(masses * averagine)

# start of the analysis ---------------------------------------------------

# Read the ms1 file into a list of peaks
# Each element of the list contains the retention time and a tuple in the form of (mz, intensity)
peaklist <- read_ms1(ms1filename)

cat("Deleting unduplicated peaks\n")
tol <- 0.005
justpeaks <- purrr::map(peaklist, 2)
justpeaks <- delete_unduplicated_peaklist(justpeaks, tol)

if (TRUE) {
	if (!dir.exists(outdir)) {
		dir.create(outdir)
	}
}

# Next we want to see how many clusters we get if we combine consecutive scans into blocks
# and cluster the combined peaks rather than cluster the peaks scan by scan
all_binbounds <- list()
ratios <- integer()
mspans <- integer()
t <- 1
toprint <- TRUE

while (t <= length(justpeaks)) {
	# browser()
	numfeats <- 0
	ov <- output_vector$new(paste0(outdir, t))
	while ((numfeats < maxfeatures) && (ov$M < maxbins)) {

		if (t > length(justpeaks)) break()

		cat("Numfeats:", numfeats, "\n")

		blockpeaks <- purrr::map(seq_len(length(justpeaks[[t]])),
								 ~c(justpeaks[[t]][[.x]][1], justpeaks[[t]][[.x]][2], t))

		if (length(blockpeaks) == 0) {
			t <- t + 1
			next()
		}

		tol <- 0.01
		clusts_cixes <- peaklist_to_cluster(blockpeaks, tol)
		clusts <- clusts_cixes[[1]]
		cixes <- clusts_cixes[[2]]

		if (toprint) {
			cat("num clusters:", length(clusts), length(clusts) / length(blockpeaks), "\n")
			cat("Hypothesizing isotope features \n")
		}

		tol <- .01
		# Each element of featues is a tuple
		# The tuple is of the form (index of peakcluster correspond to the monoisotope, charge)
		features <- construct_features_from_clusters(clusts, maxcharge, tol)

		if (toprint) cat(length(features), "features!", "\n")
		isopeaks <- run_computems1_on_features(features, outdir, ffile = formulafile, sfile = specfile,
											   masses, amass, averagine)
		justobsclusts <- length(clusts)
		if (toprint) cat("Before:", length(clusts), "\n")

		clusts <- generate_theoretical_clusts_from_observed(features, clusts, isopeaks)

		if (toprint) cat("After:", length(clusts), "\n")

		ordered_idx <- order(purrr::map_dbl(clusts, "minmz"))
		clusts <- clusts[ordered_idx]
		binbounds_cids <- assign_bins(clusts)
		binbounds <- binbounds_cids[[1]]
		cids <- binbounds_cids[[2]]
		spans <- purrr::map_dbl(binbounds, ~ .x[2] - .x[1])

		if (toprint) {
			cat("Lenbins:", length(binbounds), "# observed peak clusters:", justobsclusts, "\n")
			ratio <- length(binbounds) / justobsclusts
			ratios <- c(ratios, ratio)
			cat("Ratio of bins to observed clusters:", ratio, "\n")
			cat("Max span:", max(spans), "\n")
			mspans <- c(mspans, max(spans))
		}

		# Here, construct x into rowmajor matrices
		if (toprint) cat("Filling the matrix \n")
		N <- length(features)
		numfeats <- numfeats + N
		M <- length(binbounds)

		x <- construct_x_from_features(clusts, cids, .nrow = M, .ncol = N)

		footprints <- purrr::map(1:ncol(x), ~get_footprint(x, .x))
		ns <- do.call(order, purrr::flatten(apply(do.call(cbind, footprints), 1, list)))

		if (length(ns) == 0) next()

		to_preserve <- ns[1]
		for (i in 2: length(ns)) {
			if (!identical(footprints[[ns[i]]], footprints[[ns[i - 1]]])) {
				to_preserve <- c(to_preserve, ns[i])
			}
		}
		x <- x[, to_preserve]
		ov$output_xy(x, features[to_preserve], cids, binbounds, t)
		t <- t + 1
	}
	rm(ov)
}

