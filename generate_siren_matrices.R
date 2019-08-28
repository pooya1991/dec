library(dplyr)
library(purrr)
source("feature_generation_and_clustering.R")
source("make_averagine_models.R")

noisetol <- 5e-3
mass_accuracy <- 6e-4
maxcharge <- 6L
ms1_file <- "15c.ms1"

peaks_dfs <- read_ms1(ms1_file)
peaks_dfs <- delete_unduplicated_peaklist(peaks_dfs, noisetol)
clusters_list <- map(peaks_dfs, peaklist_to_clusters, mass_accuracy, .fun = max)
features_list <- map(clusters_list, clusters_to_features, maxcharge, mass_accuracy)
isopeaks_list <- map(features_list, run_computems1_on_features)
clusters_theo_list <- map2(features_list, isopeaks_list, generate_theoretical_clusts)
outputs_list <- map2(clusters_theo_list, clusters_list, generate_outputs)
features_idx_list <- map(outputs_list, "features_idx")
features_list <- map2(features_list, features_idx_list, ~.x[.y, ])
