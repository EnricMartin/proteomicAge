# Internal: coefficient loading and cache
# This file is loaded before other R files (alphabetical order)

.tanaka2018_cache     <- new.env(parent = emptyenv())
.lehallier2019_cache  <- new.env(parent = emptyenv())

load_tanaka2018_coefs <- function() {
  if (is.null(.tanaka2018_cache$coefs)) {
    path <- system.file("extdata", "tanaka2018_coefs.csv",
                        package = "proteomicAge", mustWork = TRUE)
    coefs <- utils::read.csv(path, stringsAsFactors = FALSE)
    intercept_row <- coefs[coefs$SOMAID == "(Intercept)", ]
    protein_rows  <- coefs[coefs$SOMAID != "(Intercept)", ]
    .tanaka2018_cache$intercept <- intercept_row$Weight
    .tanaka2018_cache$proteins  <- protein_rows
  }
  invisible()
}

load_lehallier2019_coefs <- function() {
  if (is.null(.lehallier2019_cache$coefs)) {
    path <- system.file("extdata", "lehallier2019_coefs.csv",
                        package = "proteomicAge", mustWork = FALSE)
    if (!file.exists(path)) {
      .lehallier2019_cache$proteins  <- data.frame()
      .lehallier2019_cache$intercept <- 0
      return(invisible())
    }
    coefs <- utils::read.csv(path, stringsAsFactors = FALSE,
                             comment.char = "#")
    intercept_row <- coefs[coefs$SOMAID == "(Intercept)", ]
    protein_rows  <- coefs[coefs$SOMAID != "(Intercept)" & !grepl("^#", coefs$SOMAID), ]
    .lehallier2019_cache$intercept <- if (nrow(intercept_row) > 0) intercept_row$Weight[1] else 0
    .lehallier2019_cache$proteins  <- protein_rows
  }
  invisible()
}
