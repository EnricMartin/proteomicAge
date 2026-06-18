# Internal: coefficient loading and cache
# This file is loaded before tanaka2018.R (alphabetical order)

.tanaka2018_cache <- new.env(parent = emptyenv())

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
