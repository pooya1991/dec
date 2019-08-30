for (i in seq_along(subX_list)) {
	mx1 <- outputs_list[[i]]$X
	mx1 <- as.matrix(mx1)
	idx_row <- rowSums(mx1) > 0
	mx1 <- mx1[idx_row, ]
	my1 <- outputs_list[[i]]$Y
	my1 <- as.vector(my1)
	my1 <- my1[idx_row]

	mx2 <- subX_list[[i]]
	my2 <- subY_list[[i]]

	print(all(near(drop(my2), my1)))
	print(all(near(mx1, mx2)))
}






