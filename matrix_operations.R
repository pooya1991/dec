# This function converts regular matrices to siren sparse format. When S is TRUE, a "S" character
# would be written at the first line. To understand the meaning of this, one should consult siren's
# documentation
regular_to_sparse <- function(mat, S = FALSE) {
	mat <- unname(mat)
	m <- nrow(mat)
	n <- ncol(mat)

	if (S) {
		idx_row <- which(mat[, 1] != 0)
		sp <- paste(idx_row - 1, mat[idx_row, 1], sep = "\t")
		return(c(paste(m, n, "S", sep = "\t"), sp))
	}

	sp <- purrr::map(apply(mat, 1, list), 1)
	idx_nonzero <- purrr::map(sp, ~ which(.x != 0))
	# idx_valid determines which rows have at least one nonzero element
	idx_valid <- purrr::map_lgl(idx_nonzero, ~length(.x) > 0)
	# This is intended to convert R's indices to python's indices
	idx_valid_int <- which(idx_valid) - 1
	sp <- purrr::map2(idx_nonzero[idx_valid], sp[idx_valid], ~paste0(.x - 1, "\t", .y[.x])) %>%
		purrr::map2(idx_valid_int, ~ list(paste(">", .y, sep = "\t"), .x))

	sp <- unlist(sp[idx_valid])
	c(paste0(m, "\t", n), sp)
}

# This function converts siren sparse format to regular matrices.
sparse_to_regular <- function(sp) {
	# The first line is information about the size and type of the matrix
	m_n_S <- stringr::str_split(sp[1], "\t")[[1]]
	m_n <- as.integer(m_n_S[1:2])
	sp <- sp[-1]
	mat <- matrix(0, nrow = m_n[1], ncol = m_n[2])

	if (length(m_n_S) == 3 && m_n_S[3] == "S") {
		mat_raw <- stringr::str_split(sp, "\t")
		idx_row <- purrr::map_chr(mat_raw, 1) %>%
			as.integer() %>%
			`+`(1)

		vals <- purrr::map_chr(mat_raw, 2) %>%
			as.numeric() %>%
			split(idx_row) %>%
			purrr::map_dbl(sum)

		mat[as.integer(names(vals)), 1] <- vals
		return(mat)
	}

	mat_raw <- split(sp, cumsum(stringr::str_detect(sp, ">")))
	idx_row <- purrr::map(mat_raw, 1) %>%
		stringr::str_extract("(?<=>\t)\\d+") %>%
		as.integer() %>%
		`+`(1)

	mat_raw <- purrr::map(mat_raw, `[`, -1) %>%
		purrr::map(stringr::str_split, "\t")

	idx_cols <- mat_raw %>%
		purrr::map_depth(2, 1) %>%
		purrr::map(purrr::reduce, c) %>%
		purrr::map(~as.integer(.x) + 1)

	vals <- mat_raw %>%
		purrr::map_depth(2, 2) %>%
		purrr::map(purrr::reduce, c) %>%
		purrr::map(as.numeric)

	purrr::pmap(list(idx_row, idx_cols, vals), function(x, y, z) {
		mat[x, y] <<- z
	})
	return(mat)
}

# This function returns a vector which determines the scan number for each row of a sparse matrix
get_scannum <- function(sp) {
	new_row_raw <- stringr::str_subset(sp, "^>") %>%
		stringr::str_split("\t")

	scan_num <- purrr::map_chr(new_row_raw, 5) %>% as.integer()
	row_num <- purrr::map_chr(new_row_raw, 2) %>%
		as.integer() %>% `+`(1)

	names(scan_num) <- row_num
	return(scan_num)
}

# This function returns the binbounds for each row of a sparse matrix
get_binbounds <- function(sp) {
	new_row_raw <- stringr::str_subset(sp, "^>") %>%
		stringr::str_split("\t")

	row_num <- purrr::map_chr(new_row_raw, 2) %>%
		as.integer() %>% `+`(1)

	binbounds <- new_row_raw %>%
		purrr::map(`[`, c(3, 4)) %>%
		do.call(what = rbind) %>%
		apply(2, as.numeric)

	rownames(binbounds) <- row_num
	return(binbounds)
}

# This function gives the range of nonzoro columns for each row.
get_row_col_ranges <- function(sp) {
	idx_new_row <- which(str_detect(sp, "^>"))
	idx_first_col <- idx_new_row + 1
	idx_last_col <- c((idx_new_row - 1)[-1], length(sp))


	first_cols <- sp[idx_first_col] %>%
		stringr::str_split("\t") %>%
		purrr::map_chr(1) %>%
		as.integer()

	last_cols <- sp[idx_last_col] %>%
		stringr::str_split("\t") %>%
		purrr::map_chr(1) %>%
		as.integer()

	col_ranges <- cbind(first_cols, last_cols) + 1L
	rownames(col_ranges) <- stringr::str_extract(sp[idx_new_row], "(?<=>\t)\\d+") %>%
		as.integer() %>% `+`(1)

	return(col_ranges)
}

# This function gives the range of nonzero columns for each scan
get_scan_nonzero_col_ranges <- function(sp) {
	row_col_ranges <- get_row_col_ranges(sp)
	first_cols <- row_col_ranges[, 1]
	last_cols <- row_col_ranges[, 2]
	scan_num <- get_scannum(sp)

	first_non_zero <- split(first_cols, scan_num) %>%
		purrr::map_int(min)

	last_non_zero <- split(last_cols, scan_num) %>%
		purrr::map_int(max)

	cbind(first_non_zero, last_non_zero)
}
