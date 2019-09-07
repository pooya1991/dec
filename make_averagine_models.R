get_elemental_composition <- function(features) {
    masses <- c(C = 12, H = 1.0078246, N = 14.0030732, O = 15.9949141, S = 31.972072)
    # Elemental composition of the average amino acid resideu
    averagine <- c(C = 4.9384, H = 7.7583, N = 1.3577, O = 1.4773, S = 0.0417)
    # Mass of monoisotopic averagine residue
    amass <- sum(masses * averagine)
    df <- features %>%
        mutate(meanmz = (minmz + maxmz) / 2,
               residumass = (meanmz - masses["H"]) * charge,
               aunits = residumass / amass,
               acomp = map(aunits, ~ averagine * .x)
        )

    acomp <- df$acomp %>%
        unlist() %>% round() %>% as.integer() %>%
        matrix(ncol = 5, byrow = TRUE)

    colnames(acomp) <- names(masses)

    averagine_mass <- apply(acomp, 1, function(x) sum(x * masses))
    numH <- as.integer((df$residumass - averagine_mass) / masses["H"])
    acomp[, "H"] <- acomp[, "H"] + numH
    acomp
}

run_computems1_on_features <- function(features) {
    comp <- get_elemental_composition(features)
    elements <- colnames(comp)
    comp_str <- apply(comp, 1, function(x) paste0(elements, x, collapse = "")) %>%
        paste(features$charge, sep = "\t")

    ffile <- tempfile()
    writeLines(comp_str, ffile)
    computems1_file <- ifelse(.Platform$OS.type == "windows", "./ComputeMS1.exe", "./ComputeMS1")
    sfile <- system(paste(computems1_file, ffile), intern = TRUE)

    tibble(mz_intensity = stringr::str_subset(sfile, "^(Sequence|\\d+)")) %>%
        mutate(
            idx_seq_start = stringr::str_detect(mz_intensity, "^Sequence"),
            id = cumsum(idx_seq_start)
            ) %>%
        filter(!idx_seq_start) %>%
        tidyr::separate(mz_intensity, into = c("mz", "intensity"), sep = "\t", convert = TRUE) %>%
        select(id, mz, intensity) %>%
        filter(intensity >= 0.1) %>%
        group_by(id) %>%
        mutate(intensity = intensity / sqrt(sum(intensity ^ 2))) %>%
        ungroup()
}
