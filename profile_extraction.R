sub_features_list <- map2(sub_features_list_reserve, coefs_list,
						  ~add_column(.x, coef = .y)) %>%
	map2(coefs2_list, ~add_column(.x, coef2 = .y)) %>%
	map(~ mutate(.x, hook =  near(coef, coef2, 1e-5)))

sub_features_list <- map(sub_features_list, filter, coef2 > 0) %>%
	map(~ mutate(.x, peak = map2(scan, coef, ~list(c(.x, .y))))) %>%
	map(~ mutate(.x, peak2 = map2(scan, coef2, ~list(c(.x, .y))))) %>%
	map(select, meanmz, charge, peak = peak2, hook)

features <- sub_features_list %>%
	bind_rows(.id = "scan") %>%
	mutate(scan = as.integer(scan)) %>%
	arrange(scan, meanmz) %>%
	pmap(list)

update_meanmz <- function(x) mean(x)

alignment_status <- function(feat_var, feat_fix, mass_accuracy) {
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

features_aligned <- list(c(features[[2]], reach = features[[2]][["scan"]] + 12))
for (feature_curr in features[3:length(features)]) {
	add_new <- TRUE

	for (i in seq_along(features_aligned)) {
		feature <- features_aligned[[i]]
		alignment_status_code <- alignment_status(feature, feature_curr, mass_accuracy)

		if ((alignment_status_code <= 1L) && (feature_curr$scan <= feature$reach)) {
			feature$scan <- c(feature$scan, feature_curr$scan)
			feature$meanmz <- c(feature$meanmz, feature_curr$meanmz)
			feature$charge <- c(feature$charge, feature_curr$charge)
			feature$peak <- c(feature$peak, feature_curr$peak)
			feature$hook <- c(feature$hook, feature_curr$hook)

			if (alignment_status_code == 0L && feature_curr$hook) {
				feature$reach <- feature_curr$scan + 12L
			}
			add_new <- FALSE
			features_aligned[[i]] <- feature
		}
	}

	if (add_new && feature_curr$hook) {
		features_aligned <- c(features_aligned, list(c(feature_curr, reach = feature_curr$scan + 12)))
	}
}

# profile extraction ------------------------------------------------------

idx <-  map_int(features_aligned, ~length(.x$peak)) %>%
	(function(x) x > 5)

features_aligned2 <- features_aligned[idx]
idx_ord <- map(features_aligned2, "meanmz") %>%
	map_dbl(mean) %>%
	order()

features_aligned2 <- features_aligned2[idx_ord]

profiles_mat <- matrix(NA_real_, ncol = length(features_aligned2), nrow = length(subX_list))
for (i in seq_along(features_aligned2)) {
	peaks <- features_aligned2[[i]][["peak"]]
	for (peak in peaks) {
		profiles_mat[as.integer(peak[1]) + 1, i] <- peak[2]
	}
}
