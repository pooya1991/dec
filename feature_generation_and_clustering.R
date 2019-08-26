read_ms1 <- function(ms1_file) {
	ms1_file <- readLines("15c.ms1")
	ms1_file <- stringr::str_subset(ms1_file, "^(D|H|I|Z)", negate = TRUE)
	split_scan <- stringr::str_detect(ms1_file, "^S") %>% cumsum()
	peaks_dfs <- split(ms1_file, split_scan) %>%
		map(readr::read_delim, " ", skip = 1, col_names = c("mz", "intensity"), col_types = "dd")
	peaks_dfs
}

delete_unduplicated_peaklist <- function(peaks_dfs, noisetol) {
	idx_unduplicated <- map(peaks_dfs, ~rep(FALSE, nrow(.x)))
	x <- peaks_dfs[[1]][["mz"]]
	for (t in 2:length(peaks_dfs)) {
		y <- peaks_dfs[[t]][["mz"]]
		i <- j <- 1
		I <- length(x)
		J <- length(y)

		while (i <= I && j <= J) {
			diff <- x[i] - y[j]

			if (abs(diff) <= noisetol) {
				idx_unduplicated[[t - 1]][i] <- TRUE
				idx_unduplicated[[t]][j] <- TRUE
			}

			if (diff > 0) {
				j <- j + 1
			} else {
				i <- i + 1
			}
		}
		x <- y
	}
	map2(peaks_dfs, idx_unduplicated, ~filter(.x, .y))
}

peaklist_to_clusters <- function(peaks_df, mass_accuracy, .sort = FALSE) {
	if (.sort == TRUE) {
		peaks_df <- arrange(peaks_df, mz)
	}

	peaks_list <- transpose(peaks_df)
	clusters <- vector("list", length(peaks_list))
	clusters[[1]] <- list(
		minmz = peaks_list[[1]]$mz,
		maxmz = peaks_list[[1]]$mz,
		intensity = peaks_list[[1]]$intensity
	)

	i <- 1
	for (peak in peaks_list[-1]) {
		mz <- peak$mz; intensity <- peak$intensity
		clust_curr <- clusters[[i]]
		if (clust_curr$maxmz > mz * (1 - mass_accuracy)) {
			clust_curr$maxmz <- mz
			clust_curr$intensity <- c(clust_curr$intensity, intensity)
			clusters[[i]] <- clust_curr
		} else {
			i <- i + 1
			clusters[[i]] <- list(minmz = mz, maxmz = mz, intensity = intensity)
		}
	}

	compact(clusters) %>% map(map_at, .at = "inetnsity", sum) %>% bind_rows()
}

clusters_to_features <- function(clusters, maxcharge, mass_accuracy) {
	clusters <- transpose(clusters)
	features <- vector("list", maxcharge * length(clusters))
	for (i in seq_along(clusters)) {
		clust_curr <- clusters[[i]]
		isops <- map(maxcharge:1L, ~ list(
			minmz = clust_curr$minmz + 1.003355 / .x,
			maxmz = clust_curr$maxmz + 1.003355 / .x,
			charge = .x
		))

		for (isop in isops) {

			for (clust_next in clusters[-1:-i]) {
				if (clust_next$maxmz * (1 + mass_accuracy) < isop$minmz) next()
				if (isop$maxmz * (1 + mass_accuracy) < clust_next$minmz) break()
				# There is a match. Make this feature and add it to a list
				features[[(i - 1) * maxcharge + isop$charge]] <- list(
					charge = isop$charge,
					minmz  = clust_curr$minmz * (1 - mass_accuracy),
					maxmz  = clust_curr$maxmz * (1 + mass_accuracy)
				)

				break()
			}
		}
	}

	compact(features) %>% bind_rows()
}

generate_theoretical_clusts <- function(features, isopeaks) {
	features <- transpose(features)
	clusts <- vector("list", length(features))
	iso_intensities_list <- map(isopeaks, "intensity")
	i <- 1
	for (feature in features) {
		minmzs <- (0:6) * 1.003355 / feature$charge + feature$minmz
		maxmzs <- (0:6) * 1.003355 / feature$charge + feature$maxmz
		iso_intensities <- iso_intensities_list[[i]]
		idx <- 1:min(length(minmzs), length(iso_intensities))
		clusts[[i]] <- tibble(
			minmz = minmzs[idx],
			maxmz = maxmzs[idx],
			intensity = iso_intensities[idx]
		)
		i <- i + 1
	}

	bind_rows(clusts, .id = "id") %>% mutate(id = as.integer(id))
}

generate_outputs <- function(clusts_theo, clusts_obs, features) {
	df <- clusts_obs %>%
		tibble::add_column(id = 0L, .before = TRUE) %>%
		bind_rows(clusts_theo) %>%
		arrange(minmz) %>%
		mutate(
			maxmz2 = lag(maxmz, default = first(maxmz)) %>% cummax(),
			row_num = cumsum(maxmz2 < minmz) + 1L
		) %>%
		group_by(row_num) %>%
		mutate(minmz = min(minmz), maxmz = max(maxmz)) %>%
		ungroup() %>%
		select(id, row_num, minmz, maxmz, intensity)

	binbounds <- distinct(df[c("minmz", "maxmz")])
	n <- nrow(binbounds)
	X <- filter(df, id > 0) %>%
		with(Matrix::sparseMatrix(i = row_num, j = id, x = intensity, dims = c(n, max(id))))

	Y <- filter(df, id == 0L) %>%
		with(Matrix::sparseVector(x = intensity, i = row_num, length = n))

	footprints <- df %>% filter(id > 0) %>%
		select(id, row_num) %>%
		unstack(row_num ~ id) %>%
		map(`[`, 1:4)

	idx_ord <- transpose(footprints) %>% map(unlist) %>% do.call(what = order)
	to_preserve <- rep(TRUE, length(idx_ord))
	for (i in 2:length(idx_ord)) {
		if (identical(footprints[[idx_ord[i]]], footprints[[idx_ord[i - 1]]])) {
			to_preserve[i] <- FALSE
		}
	}
	idx_ord <- idx_ord[to_preserve]

	list(X = X[, idx_ord], Y = Y, binbounds = binbounds, features_idx = idx_ord)
}
