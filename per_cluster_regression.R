for(i in 1:length(subX_list)){
	m <- subX_list[[i]]
	mzbins <- sub_binbounds_list[[i]] %>%
		rowMeans()

	y1 <- subY_list[[i]]
	rownames(y1) <- mzbins
	valid_y_rows <- as_tibble(y1) %>%
		mutate(mzbin = mzbins) %>%
		filter(V1 > 0) %>%
		pull(mzbin) %>% as.character()
	y1 <- y1[valid_y_rows,]

	rownames(m) <- mzbins
	m2 <- m[, 1:(ncol(m) / 2)]
	m2 <- m2[valid_y_rows,]

	df <- as_tibble(m2) %>%
		mutate(mzbin = valid_y_rows) %>%
		tidyr::gather("feature", "value", -mzbin) %>%
		filter(value > 0) %>%
		mutate(feature = str_remove(feature, "^V") %>% as.integer())

	graph_df <- df %>%
		select(mzbin, feature) %>%
		group_by(mzbin) %>%
		nest(.key = "features") %>%
		mutate(features = map(features, ~.x[[1]])) %>%
		mutate(len = map_int(features, length)) %>%
		filter(len > 1) %>%
		mutate(features = map(features, combn, m = 2, simplify = FALSE)) %>%
		unnest(features) %>%
		mutate(
			source = map_int(features, 1),
			target = map_int(features, 2)
		) %>%
		select(mzbin = mzbin, from = source, to = target)

	graph_df_2 <- df %>%
		select(mzbin, feature) %>%
		group_by(mzbin) %>%
		nest(.key = "features") %>%
		mutate(features = map(features, ~.x[[1]])) %>%
		mutate(len = map_int(features, length))

	g <- igraph::graph_from_data_frame(graph_df[,c(2,3)], vertices = unique(df$feature), directed = FALSE)

	ceb <- igraph::cluster_edge_betweenness(g)

	features_cluster <- stack(ceb)

	unique_clusters <- features_cluster %>%
		group_by(ind) %>%
		summarise(feature_numbers = n()) %>%
		filter(feature_numbers > 1) %>%
		pull(ind) %>% as.numeric()

	mzbin_cluster <- df %>%
		mutate(cluster = features_cluster[feature, 2]) %>%
		select(mzbin, feature, cluster)

	for (clust in unique_clusters) {
		idx_col <- filter(mzbin_cluster, cluster == clust) %>% pull(feature) %>% unique()
		idx_row <- rowSums(m2[, idx_col, drop = FALSE]) > 0
		m3 <- m2[idx_row, idx_col]
		y3 <- y1[names(y1) %in% rownames(m3), drop = FALSE]
		if(class(m3) == "matrix" && dim(m3) >= 3) {
			CVGLM <- glmnet::cv.glmnet(m3, y3, nfolds = nrow(m3), lower.limits=0)
		} else {
			next()
		}

	}
}
