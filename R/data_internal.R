# Internal: coefficient loading and cache

.tanaka2018_cache          <- new.env(parent = emptyenv())
.lehallier2019_cache       <- new.env(parent = emptyenv())
.sathyan2020_cache         <- new.env(parent = emptyenv())
.oh2023_conventional_cache <- new.env(parent = emptyenv())
.wang2024_aric_cache      <- new.env(parent = emptyenv())
.kuo2024_pac_cache        <- new.env(parent = emptyenv())
.goeminne2025_cache       <- new.env(parent = emptyenv())
.olink_map_cache          <- new.env(parent = emptyenv())

.read_clock_csv <- function(path, ...) {
  coefs <- utils::read.csv(
    path, stringsAsFactors = FALSE, colClasses = "character", ...
  )
  if ("Weight" %in% names(coefs)) {
    coefs$Weight <- as.numeric(coefs$Weight)
  }
  if ("weight" %in% names(coefs)) {
    coefs$weight <- as.numeric(coefs$weight)
  }
  coefs
}

.seqid_dot_lookup <- function(lk) {
  bare_names <- sub("^seq[.]", "", names(lk), ignore.case = TRUE)
  bare <- lk
  names(bare) <- bare_names
  prefixed <- lk
  names(prefixed) <- paste0("seq.", bare_names)
  c(bare, prefixed)
}

load_tanaka2018_coefs <- function() {
  if (is.null(.tanaka2018_cache$coefs)) {
    path <- system.file("extdata", "tanaka2018_coefs.csv",
                        package = "proteomicAge", mustWork = TRUE)
    coefs <- .read_clock_csv(path)
    .tanaka2018_cache$intercept <- coefs$Weight[coefs$SOMAID == "(Intercept)"]
    .tanaka2018_cache$proteins  <- coefs[coefs$SOMAID != "(Intercept)", ]
    col_map <- c(seqid_sl = "seqid_sl", gene = "Gene", uniprot = "UniProt", seqid_dot = "seqid_dot")
    for (mb in names(col_map)) {
      col <- col_map[mb]
      lk <- stats::setNames(.tanaka2018_cache$proteins$SOMAID,
                             .tanaka2018_cache$proteins[[col]])
      lk <- lk[!is.na(names(lk)) & names(lk) != "" & lk != ""]
      .tanaka2018_cache[[paste0("lookup_", mb)]] <- lk
      # For seqid_dot: also add "seq." prefixed keys (matching common SOMAscan column names)
      if (mb == "seqid_dot") {
        .tanaka2018_cache[["lookup_seqid_dot"]] <- .seqid_dot_lookup(lk)
      }
    }
    .tanaka2018_cache$lookup_Weight <- stats::setNames(
      .tanaka2018_cache$proteins$Weight, .tanaka2018_cache$proteins$SOMAID)
  }
  invisible()
}

load_lehallier2019_coefs <- function() {
  if (is.null(.lehallier2019_cache$coefs)) {
    path <- system.file("extdata", "lehallier2019_coefs.csv",
                        package = "proteomicAge", mustWork = TRUE)
    coefs <- .read_clock_csv(path)
    coefs <- coefs[!grepl("^#", coefs$SOMAID), ]
    .lehallier2019_cache$intercept <- coefs$Weight[coefs$SOMAID == "(Intercept)"]
    .lehallier2019_cache$proteins  <- coefs[coefs$SOMAID != "(Intercept)", ]
    col_map <- c(seqid_sl = "seqid_sl", gene = "Gene", uniprot = "UniProt", seqid_dot = "seqid_dot")
    for (mb in names(col_map)) {
      col <- col_map[mb]
      lk <- stats::setNames(.lehallier2019_cache$proteins$SOMAID,
                             .lehallier2019_cache$proteins[[col]])
      lk <- lk[!is.na(names(lk)) & names(lk) != "" & lk != ""]
      .lehallier2019_cache[[paste0("lookup_", mb)]] <- lk
      if (mb == "seqid_dot") {
        .lehallier2019_cache[["lookup_seqid_dot"]] <- .seqid_dot_lookup(lk)
      }
    }
    .lehallier2019_cache$lookup_Weight <- stats::setNames(
      .lehallier2019_cache$proteins$Weight, .lehallier2019_cache$proteins$SOMAID)
  }
  invisible()
}

load_sathyan2020_coefs <- function() {
  if (is.null(.sathyan2020_cache$coefs)) {
    path <- system.file("extdata", "sathyan2020_coefs.csv",
                        package = "proteomicAge", mustWork = FALSE)
    if (!file.exists(path)) {
      .sathyan2020_cache$proteins  <- data.frame()
      .sathyan2020_cache$intercept <- 0
      return(invisible())
    }
    coefs <- .read_clock_csv(path)
    coefs <- coefs[!grepl("^#", coefs$SOMAID), ]
    .sathyan2020_cache$intercept <- coefs$Weight[coefs$SOMAID == "(Intercept)"]
    .sathyan2020_cache$proteins  <- coefs[coefs$SOMAID != "(Intercept)", ]
    col_map <- c(seqid_sl = "seqid_sl", gene = "Gene", uniprot = "UniProt", seqid_dot = "seqid_dot")
    for (mb in names(col_map)) {
      col <- col_map[mb]
      lk <- stats::setNames(.sathyan2020_cache$proteins$SOMAID,
                             .sathyan2020_cache$proteins[[col]])
      lk <- lk[!is.na(names(lk)) & names(lk) != "" & lk != ""]
      .sathyan2020_cache[[paste0("lookup_", mb)]] <- lk
      if (mb == "seqid_dot") {
        .sathyan2020_cache[["lookup_seqid_dot"]] <- .seqid_dot_lookup(lk)
      }
    }
    .sathyan2020_cache$lookup_Weight <- stats::setNames(
      .sathyan2020_cache$proteins$Weight, .sathyan2020_cache$proteins$SOMAID)
  }
  invisible()
}

load_oh2023_conventional_coefs <- function() {
  if (is.null(.oh2023_conventional_cache$coefs)) {
    model_path <- system.file("extdata", "oh2023_conventional_models.csv",
                              package = "proteomicAge", mustWork = TRUE)
    protein_path <- system.file("extdata", "oh2023_conventional_proteins.csv",
                                package = "proteomicAge", mustWork = TRUE)
    coefs <- utils::read.csv(model_path, stringsAsFactors = FALSE, check.names = FALSE)
    proteins <- .read_clock_csv(protein_path, check.names = FALSE)
    protein_cols <- setdiff(names(coefs), c("organ", "bootstrap_seed", "y_intercept", "Sex_F"))
    .oh2023_conventional_cache$coefs <- coefs
    .oh2023_conventional_cache$proteins <- proteins
    .oh2023_conventional_cache$protein_cols <- protein_cols
    .oh2023_conventional_cache$intercept <- coefs$y_intercept
    .oh2023_conventional_cache$sex_weight <- coefs$Sex_F
    .oh2023_conventional_cache$weight_matrix <- as.matrix(coefs[, protein_cols, drop = FALSE])
    col_map <- c(seqid_sl = "seqid_sl", gene = "Gene", uniprot = "UniProt", seqid_dot = "seqid_dot")
    for (mb in names(col_map)) {
      col <- col_map[mb]
      lk <- stats::setNames(.oh2023_conventional_cache$proteins$seqid_dot,
                             .oh2023_conventional_cache$proteins[[col]])
      lk <- lk[!is.na(names(lk)) & names(lk) != "" & lk != ""]
      .oh2023_conventional_cache[[paste0("lookup_", mb)]] <- lk
      if (mb == "seqid_dot") {
        prefixed <- lk
        names(prefixed) <- sub("^seq[.]", "", names(lk), ignore.case = TRUE)
        full <- stats::setNames(.oh2023_conventional_cache$proteins$seqid_dot,
                                 .oh2023_conventional_cache$proteins$seqid_dot_full)
        full <- full[!is.na(names(full)) & names(full) != "" & full != ""]
        .oh2023_conventional_cache[["lookup_seqid_dot"]] <- c(lk, prefixed, full)
      }
    }
    .oh2023_conventional_cache$lookup_Weight <- stats::setNames(
      .oh2023_conventional_cache$proteins$Weight, .oh2023_conventional_cache$proteins$seqid_dot)
  }
  invisible()
}

load_wang2024_aric_coefs <- function() {
  if (is.null(.wang2024_aric_cache$coefs)) {
    path <- system.file("extdata", "wang2024_aric_midlife_coefs.csv",
                        package = "proteomicAge", mustWork = TRUE)
    coefs <- .read_clock_csv(path)
    missing_dot <- coefs$SOMAID != "(Intercept)" &
      (is.na(coefs$seqid_dot) | coefs$seqid_dot == "")
    coefs$seqid_dot[missing_dot] <- paste0("seq.", gsub("-", ".", coefs$SOMAID[missing_dot]))
    .wang2024_aric_cache$intercept <- coefs$Weight[coefs$SOMAID == "(Intercept)"]
    .wang2024_aric_cache$proteins  <- coefs[coefs$SOMAID != "(Intercept)", ]
    col_map <- c(seqid_sl = "seqid_sl", gene = "Gene", uniprot = "UniProt", seqid_dot = "seqid_dot")
    for (mb in names(col_map)) {
      col <- col_map[mb]
      lk <- stats::setNames(.wang2024_aric_cache$proteins$SOMAID,
                             .wang2024_aric_cache$proteins[[col]])
      lk <- lk[!is.na(names(lk)) & names(lk) != "" & lk != ""]
      .wang2024_aric_cache[[paste0("lookup_", mb)]] <- lk
      if (mb == "seqid_dot") {
        .wang2024_aric_cache[["lookup_seqid_dot"]] <- .seqid_dot_lookup(lk)
      }
    }
    .wang2024_aric_cache$lookup_Weight <- stats::setNames(
      .wang2024_aric_cache$proteins$Weight, .wang2024_aric_cache$proteins$SOMAID)
  }
  invisible()
}

load_kuo2024_pac_coefs <- function() {
  if (is.null(.kuo2024_pac_cache$coefs)) {
    coef_path <- system.file("extdata", "kuo2024_pac_coefs.csv",
                             package = "proteomicAge", mustWork = TRUE)
    param_path <- system.file("extdata", "kuo2024_pac_params.csv",
                              package = "proteomicAge", mustWork = TRUE)
    coefs <- utils::read.csv(coef_path, stringsAsFactors = FALSE)
    params <- utils::read.csv(param_path, stringsAsFactors = FALSE)
    .kuo2024_pac_cache$coefs <- coefs
    .kuo2024_pac_cache$params <- stats::setNames(params$value, params$parameter)
    .kuo2024_pac_cache$proteins <- coefs$predictor[coefs$predictor != "age"]
    .kuo2024_pac_cache$weights <- stats::setNames(coefs$weight, coefs$predictor)
  }
  invisible()
}

load_goeminne2025_coefs <- function() {
  if (is.null(.goeminne2025_cache$coefs)) {
    path <- system.file("extdata", "goeminne2025_organaging_coefs.csv",
                        package = "proteomicAge", mustWork = TRUE)
    coefs <- .read_clock_csv(path)
    .goeminne2025_cache$coefs <- coefs
    .goeminne2025_cache$organs <- sort(unique(coefs$organ))
    .goeminne2025_cache$proteins <- sort(unique(coefs$protein[coefs$protein != "Intercept"]))
  }
  invisible()
}

load_olink_protein_map <- function() {
  if (is.null(.olink_map_cache$map)) {
    path <- system.file("extdata", "olink_protein_map.csv",
                        package = "proteomicAge", mustWork = TRUE)
    .olink_map_cache$map <- utils::read.csv(path, stringsAsFactors = FALSE)
  }
  invisible()
}
