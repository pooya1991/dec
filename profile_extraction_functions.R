alignment_status <- function(feat_var, feat_fix, mass_accuracy) {
	if (feat_fix$scan %in% feat_var$scan) return(2L)
	meanmz <- mean(feat_var$meanmz)
	if (
		between(feat_fix$meanmz,
				meanmz * (1 - mass_accuracy),
				meanmz * (1 + mass_accuracy))
	) {
		if (feat_var$charge[1] == feat_fix$charge) return(0L)
		return(1L)
	}

	return(2L)
}
