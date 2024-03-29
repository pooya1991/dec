# install.packages("tidyverse", "glmnet", "igraph", "Metrics", "egg")
library(tidyverse)
source("feature_generation_and_clustering.R")
source("make_averagine_models.R")
source("regressors.R")
source("profile_extraction_functions.R")

# set parameters ----------------------------------------------------------

noisetol <- 5e-3
mass_accuracy <- 6e-6
maxcharge <- 6L
ms1_file <- "15c.ms1"

# generate features and matrices-------------------------------------------

peaks_dfs <- read_ms1(ms1_file)
peaks_dfs <- delete_unduplicated_peaklist(peaks_dfs, noisetol)
clusters_list <- map(peaks_dfs, peaklist_to_clusters_fast, mass_accuracy, .fun = sum)
features_list <- map(clusters_list, clusters_to_features, maxcharge, mass_accuracy)
isopeaks_list <- map(features_list, run_computems1_on_features)
clusters_theo_list <- map2(features_list, isopeaks_list, generate_theoretical_clusts_fast, maxcharge)
outputs_list <- map2(clusters_theo_list, clusters_list, generate_outputs, .fun = sum)
features_idx_list <- map(outputs_list, "features_idx")
features_list <- map2(features_list, features_idx_list, ~.x[.y, ])
X_list <- map(outputs_list, "X")
Y_list <- map(outputs_list, "Y")
binbounds_list <- map(outputs_list, "binbounds")

features_list <- map(features_list, ~mutate(.x, meanmz = ((minmz + maxmz) / 2)))
features_list <- map2(
	features_list, X_list,
	~add_column(.x, isop_pattern = apply(.y, 2,
										 function(x) {idx <- which(x > 0); x[idx[1]:idx[length(idx)]]})
	)
)

# clustering --------------------------------------------------------------

clustering_info_list <- vector("list", length(X_list))
for (i in seq_along(X_list)) {
	X <- X_list[[i]] %>% as.matrix()
	binbounds <- binbounds_list[[i]]
	mzbins <- rowMeans(binbounds)
	# Y <- Y_list[[i]] %>% as.vector() %>% as.matrix(ncol = 1)
	# idx_nonzero_y <- Y[, 1] > 0
	# Y <- Y[idx_nonzero_y, , drop = FALSE]
	# X <- X[idx_nonzero_y, ]
	# mzbins <- mzbins[idx_nonzero_y]

	X_df <- as_tibble(X) %>%
		mutate(mzbin = mzbins) %>%
		tidyr::gather("feature", "value", -mzbin) %>%
		filter(value > 0) %>%
		mutate(feature = str_remove(feature, "^V") %>% as.integer())

	graph_df <- X_df %>%
		select(mzbin, feature) %>%
		group_by(mzbin) %>%
		nest(.key = "features") %>%
		mutate(features = map(features, ~.x[[1]])) %>%
		mutate(len = map_int(features, length)) %>%
		filter(len > 1) %>%
		mutate(features = map(features, combn, m = 2, simplify = FALSE)) %>%
		unnest(features) %>%
		mutate(
			from = map_int(features, 1),
			to = map_int(features, 2)
		) %>%
		select(mzbin, from, to)

	g <- igraph::graph_from_data_frame(graph_df[, c(2, 3)],
									   vertices = unique(X_df$feature), directed = FALSE)

	ceb <- igraph::cluster_edge_betweenness(g)

	clustering_info <- tibble(feature = ceb$names, cluster = ceb$membership) %>%
		mutate_all(as.integer) %>%
		right_join(X_df[c("mzbin", "feature")], "feature") %>%
		select(mzbin, feature, cluster)

	clustering_info_list[[i]] <- clustering_info
}

# regression --------------------------------------------------------------

coefs_list <- vector("list", length(X_list))
coefs2_list <- vector("list", length(X_list))
for (i in seq_along(X_list)) {
	X <- X_list[[i]] %>% as.matrix()
	Y <- Y_list[[i]] %>% as.vector() %>% as.matrix(ncol = 1)
	clustering_info <- clustering_info_list[[i]]
	coefs <- regressor(X, Y, clustering_info)
	coefs2_list[[i]] <- coefs2
	coefs_list[[i]] <- coefs
}

features_list <- map2(features_list, coefs_list,
						  ~add_column(.x, coef = .y)) %>%
	map2(coefs2_list, ~add_column(.x, coef2 = .y)) %>%
	map(~ mutate(.x, anchor =  near(coef, coef2, 1e-5)))

# regression results presentation -----------------------------------------

if (!dir.exists("./regression_plots")) dir.create("./regression_plots")
drop_zero_rows <- TRUE
for (i in seq_along(X_list)) {
	X <- X_list[[i]] %>% as.matrix()
	Y <- Y_list[[i]] %>% as.vector() %>% as.matrix(ncol = 1)
	binbounds <- binbounds_list[[i]]
	mzbins <- rowMeans(binbounds)
	idx_nonzero_y <- Y[, 1] > 0
	Y <- Y[idx_nonzero_y, , drop = FALSE]
	X <- X[idx_nonzero_y, ]
	mzbins <- mzbins[idx_nonzero_y]
	features_info <- features_list[[i]] %>%
		select(meanmz, charge) %>%
		mutate(feature = row_number())

	clustering_info <- clustering_info_list[[i]] %>%
		filter(mzbin %in% mzbins)

	feature_cluster <- distinct(clustering_info[c("feature", "cluster")])

	valid_clusters <- feature_cluster %>%
		count(cluster) %>%
		filter(n > 1) %>% pull(cluster) %>% as.integer()

	coefs <- coefs_list[[i]]
	Y_hat <- X %*% matrix(coefs, ncol = 1)
	act_pred <- cbind(mzbins, Y, Y_hat)
	deconv_mat <- apply(X, 1, "*", coefs) %>% t()

	if (drop_zero_rows) {
		idx_nonzero_y_hat <- Y_hat[, 1] > 0
		act_pred <- act_pred[idx_nonzero_y_hat, ]
		deconv_mat <- deconv_mat[idx_nonzero_y_hat, ]
	}

	if (nrow(deconv_mat) == 0) next()

	deconv_long <- as_tibble(deconv_mat) %>%
		mutate(mzbin = act_pred[, 1]) %>%
		gather(feature, intensity, -mzbin) %>%
		filter(intensity > 0) %>%
		mutate(feature = str_remove(feature, "^V") %>% as.integer()) %>%
		left_join(feature_cluster, "feature")

	for (clust in valid_clusters) {
		mzbins_sub <- clustering_info %>%
			filter(cluster == clust) %>% pull(mzbin) %>% intersect(act_pred[, 1])

		deconv_long_sub <- deconv_long %>%
			filter(cluster == clust)

		if (length(mzbins_sub) == 0) next()
		act_pred_sub <- act_pred[act_pred[, 1] %in% mzbins_sub, , drop = FALSE]
		rmse_sub <- Metrics::rmse(act_pred_sub[ ,2], act_pred_sub[, 3])
		mape_sub <- Metrics::mape(act_pred_sub[ ,2], act_pred_sub[, 3])
		perc_err = ((act_pred_sub[, 2] - act_pred_sub[, 3]) / act_pred_sub[, 2]) * 100

		legend_data <- filter(features_info, feature %in% unique(deconv_long_sub$feature)) %>%
			select(feature, meanmz, charge)

		fake_data <-tibble(mzbin = mzbins_sub,
						   feature = deconv_long_sub$feature[1],
						   intensity = 0)

		ylim = c(0,  max(act_pred_sub[, -1]))
		label <- tibble(
			mzbin = Inf,
			intensity = Inf,
			label = paste0("scan = ", i, " cluster = ", clust)
		)

		p1 <- tibble(mzbin = act_pred_sub[, 1], intensity = act_pred_sub[, 2], perc_err = perc_err) %>%
			ggplot(aes(mzbin, intensity)) +
			geom_col() +
			geom_text(aes(label = paste0(round(perc_err, 1), "%")), vjust = -0.3, size = 3) +
			geom_text(aes(label = label), data = label, vjust = "top", hjust = "right") +
			coord_cartesian(ylim = ylim, expand = TRUE) +
			labs(title = "Actual")

		label <- tibble(
			mzbin = Inf,
			intensity = Inf,
			label = paste0("RMSE = ", round(rmse_sub), "\nMAPE = ", round(mape_sub, 3))
		)

		p2 <- deconv_long_sub %>%
			select(mzbin, feature, intensity) %>%
			bind_rows(fake_data) %>%
			group_by(mzbin, feature) %>%
			summarise(intensity = sum(intensity)) %>%
			filter(!is.na(feature)) %>%
			ggplot(aes(mzbin, intensity)) +
			geom_col(aes(fill = factor(feature))) +
			labs(fill = "Feature", title = "Predicted", y = NULL) +
			geom_text(aes(label = label), data = label, vjust = "top", hjust = "right") +
			scale_fill_discrete(
				breaks = legend_data$feature,
				labels = pmap(legend_data, paste, sep = ":")
			) +
			coord_cartesian(ylim = ylim, expand = TRUE) +
			theme(legend.title = element_text(size = 8),
				  legend.text = element_text(size = 8),
				  axis.ticks.y = element_blank(),
				  axis.title.y = element_blank(),
				  axis.text.y = element_blank())

		p <- egg::ggarrange(p1, p2, nrow = 1)
		plot_name <- paste0(i, "_", clust, ".png")
		ggsave(plot_name, plot = p, path = "./regression_plots", width = 16, height = 8)
	}
}

# feature construction ----------------------------------------------------

features <- map(features_list, filter, coef2 > 0) %>%
	map(select, charge, meanmz, coef2, anchor) %>%
	bind_rows(.id = "scan") %>%
	mutate(scan = as.integer(scan)) %>%
	arrange(scan, meanmz) %>%
	mutate(peak = map2(scan, coef2, ~ list(c(.x, .y)))) %>%
	select(-coef2)

# feature alignment -------------------------------------------------------

idx_first_anchor <- detect_index(features$anchor, ~.x)
features <- transpose(features)
features_aligned <- list(c(features[[idx_first_anchor]], reach = features[[idx_first_anchor]][["scan"]] + 12))
for (feature_curr in features[(idx_first_anchor + 1):length(features)]) {
	add_new <- TRUE

	for (i in seq_along(features_aligned)) {
		feature <- features_aligned[[i]]
		if (feature_curr$scan <= max(feature$scan)) next()
		alignment_status_code <- alignment_status(feature, feature_curr, mass_accuracy)

		if ((alignment_status_code <= 1L) && (feature_curr$scan <= feature$reach)) {
			feature$scan <- c(feature$scan, feature_curr$scan)
			feature$meanmz <- c(feature$meanmz, feature_curr$meanmz)
			feature$charge <- c(feature$charge, feature_curr$charge)
			feature$peak <- c(feature$peak, feature_curr$peak)
			feature$anchor <- c(feature$anchor, feature_curr$anchor)

			if (alignment_status_code == 0L && feature_curr$anchor) {
				feature$reach <- feature_curr$scan + 12L
			}
			add_new <- FALSE
			features_aligned[[i]] <- feature
			break()
		}
	}

	if (add_new && feature_curr$anchor) {
		features_aligned <- c(features_aligned, list(c(feature_curr, reach = feature_curr$scan + 12)))
	}
}

# profile extraction ------------------------------------------------------

profiles <- map(features_aligned, ~as_tibble(.x) %>% arrange(desc(scan)) %>%
						  	mutate(trailing_non_anchor = as.logical(cumsum(anchor))) %>%
						  	filter(trailing_non_anchor) %>%
						  	select(-reach, -trailing_non_anchor) %>%
						  	arrange(scan) %>%
						  	as.list()
						 )

idx <-  map_int(profiles, ~length(.x$peak)) %>%
	(function(x) x > 4)

profiles <- profiles[idx]

idx_ord <- map(profiles, "meanmz") %>%
	map_dbl(mean) %>%
	order()

profiles <- profiles[idx_ord]

profiles_mat <- matrix(NA_real_, ncol = length(profiles), nrow = length(X_list))
for (i in seq_along(profiles)) {
	peaks <- profiles[[i]][["peak"]]
	for (peak in peaks) {
		profiles_mat[as.integer(peak[1]), i] <- peak[2]
	}
}

# profiles representation -------------------------------------------------

anchors_mat <- matrix(NA, ncol = length(profiles), nrow = length(X_list))
for (i in seq_along(profiles)) {
  scans <- profiles[[i]][["scan"]]
  anchors_mat[scans, i] <- profiles[[i]][["anchor"]]
}

input <- 0.5
while (input > 0) {
	input <- readline(prompt="Enter profile number (Enter 0 to exit): " ) %>% as.integer()
	if (input > 0) {
		plot(profiles_mat[ , input],
			 type = "l", lty = 2,
			 xlab = "Scan", ylab = "Intensity")
		points(profiles_mat[ , input],
			col = c("blue", "red")[anchors_mat[, input] + 1], pch = 19, type = "b",
		)
	}
}
