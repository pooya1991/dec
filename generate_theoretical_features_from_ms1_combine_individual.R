read_ms1 <- function(ms1filename) {
    cat("Reading the ms1file\n")
    scannum <- -1
    peaklist <- list()
    rtime <- 0
    numpeaks <- 0
    ms1file <- readLines(ms1filename)
    for (line in ms1file) {
        chr1 <- str_sub(line, 1, 1)
        if (chr1 %in% c("H", "Z", "D")) next()
        if (chr1 == "I") {
            if (str_detect(line, "RTime")) {
                rtime <- str_extract(line, "(?<=RTime\t).*$") %>% as.numeric()
                peaklist <- c(peaklist, list(list(rtime, list())))
            }
            next()
        }
        if (chr1 == "S") {
            if (scannum %% 200 == 0) cat("Reading scan", scannum, "\n")
            scannum <- scannum + 1
            next()
        }
        if (str_length(line) > 1) {
            tokens <- str_split(line, " ")[[1]]
            mz <- as.numeric(tokens[1])
            intensity <- as.numeric(tokens[2])
            numpeaks <- numpeaks + 1
            peaklist[[length(peaklist)]][[2]] <- c(peaklist[[length(peaklist)]][[2]], list(c(mz, intensity)))
        }
    }
    cat("Numpeaks:", numpeaks, "\n")
    return(peaklist)
}

# For each one of the input clusters, this function generate a monoisotopic precursor ion of each charge
# It only create the feature if both the monoisotopic and the 2nd isotopic peak is present
construct_features_from_clusters <- function(clusts, maxcharge, tol) {
    N <- length(clusts)
    features <- list()

    # For each one of these clusters, generate a monoisotopic precursor ion of each charge
    # Only create the feature if both the monoisotopic and the 2nd isotopic peak is present
    for (n in seq_len(N)) {
        p <- clusts[[n]]
        n2 <- n + 1L

        # There must be a peakcluster that matches isop
        # where isop is the second isotopic peak
        for (charge in maxcharge:1L) {
            isop <- PeakCluster$new(p$minmz + 1 / charge, p$maxmz + 1 / charge, 0L, 0L)
            while (n2 <= N) {
                compare <- compare_clusters(clusts[[n2]], isop, tol)
                if (compare > 0) n2 <- n2 + 1L
                if (compare < 0) break()
                # There is a match. Make this feature and add it to a list
                if (compare == 0) {
                    features <- c(features, list(c(n, charge, p$minmz - tol, p$maxmz + tol)))
                    break()
                }
            }

        }
    }
    return(features)
}

# Here, construct x into rowmajor matrices
construct_x_from_features <- function(features, cids, .nrow, .ncol){
    x <- matrix(0, nrow = .nrow, ncol = .ncol)
    n <- 1

    for (cid in cids) {
        # This has both observed peaks and theoretical peaks
        pc <- features[[n]]
        # It's a theoretical peak
        if (pc$idnum < 0) {
            x[cid, get_n(pc$idnum)] <- pc$peaks[[1]][2]
        }
        n <- n + 1
    }
    return(x)
}

# This outputs a series of x and y matrices in sparse format
# The Y matrix is one long vector
# The X matrix is one big matrix
output_vector <- R6::R6Class(
    "output_vector",
    public = list(
        M = 0,
        N = 0,
        suffix = NULL,
        xout = NULL,
        yout = NULL,
        annout = NULL,
        initialize = function(outfilesuffix) {
            cat("\nOPENING", outfilesuffix, "\n\n")
            self$suffix <- outfilesuffix
            self$xout <- file(paste0(outfilesuffix, "_x.txt"), open = "w", raw = TRUE)
            self$yout <- file(paste0(outfilesuffix, "_y.txt"), open = "w", raw = TRUE)
            self$annout <- file(paste0(outfilesuffix, "_annfile.txt"), open = "w", raw = TRUE)
            write("t\tmin m/z\tmax m/z\tcharge" ,self$annout)
        },

        # X starts out as
        # here, the Y doesn't even contain the binbounds or the row number
        # because it's 1 dimensional!
        output_xy = function(x, features, cids, binbounds, t) {
            I_J <- apply(which(t(x) != 0, arr.ind = TRUE), 1, rev)
            I <- I_J[1,]
            J <- I_J[2,]
            previ <- 0
            for (n in seq_along(I)) {
                i <- I[n]
                if (previ != i) {
                    write(
                        paste0(">\t",(i + self$M), "\t", binbounds[[i]][1], "\t", binbounds[[i]][1],
                               "\t", t),
                        self$xout, append = TRUE
                    )
                }

                previ <- i
                j <- J[n]
                write(paste0((j + self$N), "\t", x[i, j]), self$xout, append = TRUE)
            }

            self$N <- self$N + ncol(x)

            for (vec in features) {
                n <- vec[1]; charge <- vec[2]; mz1 <- vec[3]; mz2 <- vec[4]
                write(paste0(t, "\t", mz1, "\t", mz2, "\t", charge), self$annout, append = TRUE)
            }

            n <- 1

            for (cid in cids) {
                # this has both observed peaks and theoretical peaks
                pc <- clusts[[n]]
                # it's an observed peak
                if (pc$idnum >= 0) {
                    npeaks <- length(pc$peaks)
                    peaks2 <- pc$peaks[do.call(order, purrr::flatten(apply(do.call(cbind, pc$peaks), 1, list)))]
                    for (i in 1:npeaks) {
                        if ((i == npeaks) || (peaks2[[i]][1] != peaks2[[i + 1]][1])) {
                            write(paste0((self$M + cid), "\t", peaks2[[i]][2]), self$yout, append = TRUE)
                            if (peaks2[[i]][1] < t) {
                                pc$p()
                                cat(pc$peaks)
                            }
                        }
                    }
                }
                n <- n + 1
            }
            self$M <- self$M + nrow(x)
        },

        finalize = function() {
            cat("\n\nFINISHING AND CLOSING!!!\n")
            cat("M:", self$M, "N:", self$N, "\n\n")

            close(self$xout)
            temp_file <- paste0(self$suffix, "_tempx.txt")
            xname <- paste0(self$suffix, "_x.txt")
            fin <- readLines(xname)
            write(paste0(self$M, "\t", self$N), temp_file)
            write(fin, temp_file, append = TRUE)
            file.remove(xname)
            file.rename(temp_file, xname)

            close(self$yout)
            temp_file <- paste0(self$suffix, "_tempy.txt")
            yname <- paste0(self$suffix, "_y.txt")
            fin <- readLines(yname)
            # The S means that the y matrix is in sparse format
            write(paste0(self$M, "\t1\tS"), temp_file)
            write(fin, temp_file, append = TRUE)
            file.remove(yname)
            file.rename(temp_file, yname)

            close(self$annout)
        }
    )
)

get_ix <- function(n) -n - 1

get_n <- function(ix) -ix - 1

get_footprint <- function(x1, n) {
    nz <- c(which(x[, n] != 0), 0, 0, 0, 0)
    return(nz[1:4])
}
