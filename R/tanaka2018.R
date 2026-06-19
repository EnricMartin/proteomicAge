#' Tanaka et al. (2018) 76-protein proteomic aging clock
#'
#' Implements the proteomic aging clock from Tanaka T, et al. (2018).
#' Aging Cell 17(5):e12799. Coefficients from eLife (2020) Supplementary File 1D.
#'
#' Elastic net regression (alpha=0.5, lambda=0.8767859), log2 transform.
#' Input columns must be named by UniProt accession (e.g., "P36222").
#'
#' @name tanaka2018
#' @keywords internal
"_PACKAGE"

# ---------- tanaka2018_proteins ----------

#' List Tanaka 2018 clock proteins with coefficients and UniProt IDs
#' @return data.frame with SOMAID, Gene, UniProt, Weight
#' @export
tanaka2018_proteins <- function() {
  load_tanaka2018_coefs()
  .tanaka2018_cache$proteins
}

# ---------- compute_tanaka2018_age ----------

#' Compute proteomic age via Tanaka 2018 clock
#'
#' @param data data.frame with UniProt-named protein columns + Age column
#' @param id_col Sample ID column name
#' @param age_col Chronological age column name
#' @param log_transform Apply log2 transform (default TRUE)
#' @return data.frame with proteomic_age, age_acceleration, etc.
#' @export
compute_tanaka2018_age <- function(data,
                                    id_col = "SampleID",
                                    age_col = "Age",
                                    log_transform = TRUE) {

  if (!is.data.frame(data)) stop("'data' must be a data.frame")
  if (!id_col %in% names(data)) stop("id_col '", id_col, "' not found")
  if (!age_col %in% names(data)) stop("age_col '", age_col, "' not found")

  load_tanaka2018_coefs()
  up_lookup <- .tanaka2018_cache$uniprot_lookup
  weight_lk <- .tanaka2018_cache$weight_lookup
  intercept <- .tanaka2018_cache$intercept

  # Match user's UniProt columns to clock proteins
  user_cols <- setdiff(names(data), c(id_col, age_col, "Sex"))
  matched_up <- intersect(user_cols, names(up_lookup))
  unmatched   <- setdiff(names(up_lookup), user_cols)

  if (length(unmatched) > 0) {
    warning(length(unmatched), " clock proteins not found in data. Missing UP: ",
            paste(head(unmatched, 10), collapse = ", "),
            if (length(unmatched) > 10) "..." else "")
  }

  ids <- data[[id_col]]
  chron_age <- as.numeric(data[[age_col]])
  n <- nrow(data)
  prot_age <- numeric(n)

  for (i in seq_len(n)) {
    pred <- intercept
    row <- data[i, ]
    for (up in matched_up) {
      val <- row[[up]]
      if (!is.na(val) && is.numeric(val)) {
        if (log_transform && val > 0) {
          val <- log2(val)
        } else if (log_transform && val <= 0) {
          val <- NA_real_
        }
        if (!is.na(val)) {
          sid <- up_lookup[up]
          pred <- pred + weight_lk[sid] * val
        }
      }
    }
    prot_age[i] <- pred
  }

  fit <- stats::lm(prot_age ~ chron_age)
  age_accel <- as.numeric(stats::residuals(fit))

  data.frame(
    id                = ids,
    chronological_age = chron_age,
    proteomic_age     = prot_age,
    age_acceleration  = age_accel,
    n_proteins_matched = length(matched_up),
    n_proteins_missing = length(unmatched),
    stringsAsFactors   = FALSE
  )
}
