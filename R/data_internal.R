# Internal: coefficient loading and cache

.tanaka2018_cache     <- new.env(parent = emptyenv())
.lehallier2019_cache  <- new.env(parent = emptyenv())

load_tanaka2018_coefs <- function() {
  if (is.null(.tanaka2018_cache$coefs)) {
    path <- system.file("extdata", "tanaka2018_coefs.csv",
                        package = "proteomicAge", mustWork = TRUE)
    coefs <- utils::read.csv(path, stringsAsFactors = FALSE)
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
        prefixed <- lk
        names(prefixed) <- paste0("seq.", names(lk))
        .tanaka2018_cache[["lookup_seqid_dot"]] <- c(lk, prefixed)
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
    coefs <- utils::read.csv(path, stringsAsFactors = FALSE)
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
        prefixed <- lk
        names(prefixed) <- paste0("seq.", names(lk))
        .lehallier2019_cache[["lookup_seqid_dot"]] <- c(lk, prefixed)
      }
    }
    .lehallier2019_cache$lookup_Weight <- stats::setNames(
      .lehallier2019_cache$proteins$Weight, .lehallier2019_cache$proteins$SOMAID)
  }
  invisible()
}
