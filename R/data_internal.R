# Internal: coefficient loading and cache

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
    up_lookup <- stats::setNames(protein_rows$SOMAID, protein_rows$UniProt)
    .tanaka2018_cache$uniprot_lookup <- up_lookup[!is.na(names(up_lookup)) & names(up_lookup) != ""]
    .tanaka2018_cache$weight_lookup  <- stats::setNames(protein_rows$Weight, protein_rows$SOMAID)
  }
  invisible()
}

load_lehallier2019_coefs <- function() {
  if (is.null(.lehallier2019_cache$coefs)) {
    path <- system.file("extdata", "lehallier2019_coefs.csv",
                        package = "proteomicAge", mustWork = TRUE)
    coefs <- utils::read.csv(path, stringsAsFactors = FALSE)
    intercept_row <- coefs[coefs$SOMAID == "(Intercept)", ]
    protein_rows  <- coefs[coefs$SOMAID != "(Intercept)", ]
    protein_rows  <- protein_rows[!grepl("^#", protein_rows$SOMAID), ]
    .lehallier2019_cache$intercept <- intercept_row$Weight[1]
    .lehallier2019_cache$proteins  <- protein_rows
    up_lookup <- stats::setNames(protein_rows$SOMAID, protein_rows$UniProt)
    .lehallier2019_cache$uniprot_lookup <- up_lookup[!is.na(names(up_lookup)) & names(up_lookup) != ""]
    .lehallier2019_cache$weight_lookup  <- stats::setNames(protein_rows$Weight, protein_rows$SOMAID)
  }
  invisible()
}
