delete_unduplicated_peaklist <- function(peaklists, tol = 0.005) {
    svs <- vector("list", length(peaklists)) %>%
        purrr::map2(1:length(peaklists), ~vector("list", 0))

    for (t in 1:(length(svs) - 1)) {
        if (t == 1) {
            a <- peaklists[[t]]
        } else {
            a <- b
        }

        b <- peaklists[[t + 1]]

        if (t %% 100 == 1) cat("t:", t, length(svs[[t]]), length(a), length(b), "\n")

        i <- 1
        j <- 1
        I <- length(a)
        J <- length(b)

        while (i <= I && j <= J) {
            diff <- a[[i]][1] - b[[j]][1]
            if (abs(diff) <= tol) {
                svs[[t]] <- c(svs[[t]], list(a[[i]]))
                svs[[t + 1]] <- c(svs[[t + 1]],list(b[[j]]))
            }
            if (diff > 0) {
                j <- j + 1
            } else {
                i <- i + 1
            }
        }
        svs[[t]] <- unique(svs[[t]])
        ordered_idx <- order(purrr::map_dbl(svs[[t]], 1))
        svs[[t]] <- svs[[t]][ordered_idx]
    }
    return(svs)
}

peaklist_to_cluster <- function(peaklist, tol = .005) {
    # the sorting at the first line is omited because peaklist is already sorted
    clusters <- list()
    mz <- peaklist[[1]][1]
    intensity <- peaklist[[1]][2]
    realt <- peaklist[[1]][3]
    clusters <- c(clusters, PeakCluster$new(mz, mz, realt, intensity))

    cids <- rep(1, length(peaklist))
    ix <- 2
    for (vec in peaklist[-1]) {
        mz <- vec[1]; intensity <- vec[2]; realt <- vec[3]
        if (clusters[[length(clusters)]]$maxmz > (mz - tol)) {
            clusters[[length(clusters)]]$maxmz <- mz
            cids[ix] <- cids[ix - 1]
            clusters[[length(clusters)]]$peaks <- c(clusters[[length(clusters)]]$peaks, list(c(realt, intensity)))
        } else {
            clusters <- c(clusters, PeakCluster$new(mz, mz, realt, intensity))
            cids[ix] <- cids[ix - 1] + 1
        }
        ix <- ix + 1
    }
    return(list(clusters, cids))
}

PeakCluster <- R6::R6Class(
    "PeakCluster", lock_objects = FALSE,
    private = list(
        shared_env = new.env(),
        get_numpcs = function() {
            .numpcs <- private$shared_env$.numpcs
            if (is.null(.numpcs)) .numpcs <- 0
            return(.numpcs)
        }
    ),
    active = list(
        numpcs = function(values) {
            if (missing(values)) {
                private$shared_env$.numpcs
            } else {
                stop("`numpcs` is read-only", call. = FALSE)
            }
        }
    ),
    public = list(
        minmz = NULL,
        maxmz = NULL,
        peaks = NULL,
        idnum = NULL,
        initialize = function(mz1, mz2, t, intensity) {
            self$minmz <- mz1
            self$maxmz <- mz2
            self$peaks <- list(c(t, intensity))
            private$shared_env$.numpcs <- private$get_numpcs() + 1
            self$idnum <- private$get_numpcs()
        },

        add = function(mz, tol) {
            if (mz + tol >= self$minmz) {
                if (mz - tol <= self$maxmz) {
                    self$maxmz <- max(self$maxmz, mz)
                    self$minmz <- min(self$minmz, mz)
                    return(TRUE)
                }
            }
            return(FALSE)
        },

        p = function() {
            cat(private$shared_env$.numpcs, " (", self$minmz, ", ", self$maxmz, ")", sep = "")
        },

        pfile = function(fout) {
            write(paste(self$minmz, "\t", self$maxmz, "\n", sep = ""), fout)
        }
    )
)

compare_clusters <- function(a, b, tol) {
    # b is greater than a
    if (a$maxmz + tol < b$minmz) return(1)
    #b is less than a
    if (b$maxmz + tol < a$minmz) return(-1)
    return(0) # they overlap
}

assign_bins <- function(clist) {
    binbounds <- list()
    N <- length(clist)
    cids <- rep(1, N)
    startmz <- clist[[1]]$minmz
    endmz <- clist[[1]]$maxmz

    for (n in 2:N) {
        if (clist[[n]]$minmz < endmz) {
            cids[n] <- cids[n - 1]
            endmz <- max(endmz, clist[[n]]$maxmz)
        } else {
            binbounds <- c(binbounds, list(c(startmz, endmz)))
            cids[n] <- cids[n-1] + 1
            startmz <- clist[[n]]$minmz
            endmz <- clist[[n]]$maxmz
        }
    }
    binbounds <- c(binbounds, list(c(startmz, endmz)))
    return(list(binbounds, cids))
}

segregate_annotations_by_time <- function(anns) {
    anns <- as.data.frame(anns)
    ts <- as.character(anns[,"t"])
    split(anns, ts)
}

precursor <- R6::R6Class(
    "precursor", lock_objects = FALSE,
    private = list(
        shared_env = new.env(),
        get_numpcs = function() {
            .numpcs <- private$shared_env$.numpcs
            if (is.null(.numpcs)) .numpcs <- 0
            return(.numpcs)
        }
    ),

    active = list(
        numpcs = function(values) {
            if (missing(values)) {
                private$shared_env$.numpcs
            } else {
                stop("`numpcs` is read-only", call. = FALSE)
            }
        }
    ),

    public = list(
        minmz = NULL,
        maxmz = NULL,
        charge = NULL,
        idnum = NULL,
        startscan = NULL,
        endscan = NULL,
        indices = NULL,
        minmass = NULL,
        to_output = TRUE,

        initialize = function(mz1, mz2, charge, t, index) {
            self$minmz <- mz1
            self$maxmz <- mz2
            self$charge <- charge
            private$shared_env$.numpcs <- private$get_numpcs() + 1
            self$idnum <- private$get_numpcs()
            self$startscan <- t
            self$endscan <- t
            self$indices <- index
            self$minmass <- (mz1 - massHMono) * charge
            self$maxmass <- (mz2 - massHMono) * charge
        },
        p = function() {
            text <- paste0("ID: ", self$idnum, ", charge: ", self$charge, ", m/z: ", self$minmz
                           , " - ", self$maxmz, ", scannums: ", self$startscan, ", ", self$endscan)
            cat(text, "\n")
        },

        incorporate = function(p) {

            if (self$idnum == p$idnum) return()

            if (self$startscan < p$startscan) {
                p$indices = c(rep.int(NA, p$startscan - self$startscan), p$indices)
            } else if (self.startscan > p$startscan){
                self$indices = c(rep.int(NA, self$startscan - p$startscan), self$indices)
            }

            if (self$endscan > p$endscan){
                p$indices <- c(p$indices, rep.int(NA, self$endscan - p$endscan))
            } else if (self$endscan < p$endscan){
                self$indices = c(self$indices, rep.int(NA, p$endscan - self$endscan))
            }

            p$idnum = self$idnum
            self$minmz = min(self$minmz, p$minmz)
            self$maxmz = max(self$maxmz, p$maxmz)
            self$startscan = min( self$startscan, p$startscan)
            self$endscan = max(self$endscan, p$endscan)

            if(FALSE) {
                #if self.idnum == 878:
                #if None in self.indices or None in p.indices:
                cat('Padded indices\n')
                cat(self$indices, "\n")
                cat(p$indices, "\n")
                p(self)
            }

            # Combine the padded indices
            newEnd <- self$endscan - self$startscan + 1
            newindices = rep.int(NA, newEnd)
            i = 1
            j = 1
            I = length(self$indices)
            J = length(p$indices)

            for (t in 1:newEnd) {
                z = TRUE
                if(!is.na(self$indices[i])) {
                    newindices[t] = self$indices[i]
                    z = FALSE
                    i = i + 1
                    while(i <= I & !is.na(self$indices[i]) & self$indices[i] < 0){
                        newindices[t] = self$indices[i]
                        i = i + 1
                    }
                } else {
                    i = i + 1
                }
                if(is.na(p$indices[j])) {
                    if(z) {
                        newindices[t] = p$indices[j]
                    } else {
                        newindices[t] = -p$indices[j]
                    }
                    j = j + 1
                    while(j <= J & !is.na(p$indices[i]) & p$indices[i] < 0) {
                        newindices[t] = p$indices[j]
                        j = j + 1
                    }
                } else {
                    j = j + 1
                }
            }
            self$indices <- NULL
            for (t in 1:newEnd) {
                if(is.na(newindices)){
                    self$indices <- NA
                } else {
                    self$indices <- newindices[t]
                }
            }
        }
    )
)

generate_theoretical_clusts_from_observed <- function(features, clusts, isopeaks) {
    fn <- 1
    for (vec in features) {
        n <- vec[1]; charge <- vec[2]; mz1 <- vec[3]; mz2 <- vec[4]
        mz1s <- purrr::map_dbl(0:6, ~ .x / charge + mz1)
        mz2s <- purrr::map_dbl(0:6, ~ .x / charge + mz2)
        ipeaks <- isopeaks[[fn]]
        for (i in 1:min(length(ipeaks), length(mz1s))) {
            clusts <- c(clusts, PeakCluster$new(mz1s[i], mz2s[i], t, ipeaks[[i]][2]))
            clusts[[length(clusts)]]$idnum <- get_ix(fn)
        }
        fn <- fn + 1
    }
    return(clusts)
}
