#' Lehallier et al. (2020) 491-protein ultra-predictive aging clock
#'
#' Implements the ultra-predictive LASSO aging clock from:
#' Lehallier B, Shokhirev MN, Wyss-Coray T, Johnson AA. "Data mining of human
#' plasma proteins generates a multitude of highly predictive aging clocks
#' that reflect different aspects of aging." Aging Cell. 2020;19(11):e13256.
#' doi: 10.1111/acel.13256
#'
#' The clock uses LASSO regression on 2,978 SOMAscan proteins, selecting
#' 491 SOMAmer entries. Training: n=2,178, Test: n=1,123.
#' Performance: r=0.98 (training), r=0.96 (test), MAE=1.84 years.
#'
#' @references
#' Lehallier B, et al. (2020) Aging Cell 19(11):e13256.
#' Table S7 (protein list) and Table S8 (coefficients + intercept).
#'
#' @name lehallier2020
#' @keywords internal
NULL

# ---------------------------------------------------------------------------
# lehallier2020_proteins()
# ---------------------------------------------------------------------------

#' List proteins in the Lehallier 2020 clock
#'
#' Returns the 491 SOMAmer entries and their LASSO coefficients
#' from the ultra-predictive Lehallier et al. (2020) aging clock.
#'
#' @return A data.frame with columns: \code{SOMAID}, \code{UniProt},
#'   \code{Gene}, \code{Weight}
#' @export
lehallier2020_proteins <- function() {
  load_lehallier2020_coefs()
  .lehallier2020_cache$proteins
}

# ---------------------------------------------------------------------------
# compute_lehallier2020_age()
# ---------------------------------------------------------------------------

#' Compute proteomic age using Lehallier et al. (2020) ultra-predictive clock
#'
#' Applies the 491-SOMAmer LASSO model to SOMAscan proteomic data.
#' Input values are log-transformed (log2 by default, matching common
#' SOMAscan preprocessing). This clock achieved r=0.98 in the training set
#' and r=0.96 in the test set.
#'
#' @param data A data.frame of SOMAscan proteomic data.
#' @param id_col Name of the sample identifier column.
#' @param age_col Name of the chronological age column.
#' @param log_transform Logical. If TRUE (default), log2-transform values.
#'
#' @return A data.frame with predicted proteomic age and age acceleration.
#' @export
compute_lehallier2020_age <- function(data,
                                       id_col = "SampleID",
                                       age_col = "Age",
                                       log_transform = TRUE) {

  if (!is.data.frame(data)) stop("'data' must be a data.frame")
  if (!id_col %in% names(data)) stop("id_col not found")
  if (!age_col %in% names(data)) stop("age_col not found")

  load_lehallier2020_coefs()
  proteins <- .lehallier2020_cache$proteins
  intercept <- .lehallier2020_cache$intercept

  if (is.null(proteins) || nrow(proteins) == 0) {
    stop(
      "Lehallier 2020 coefficients not available.\n",
      "Download Tables S7 and S8 from:\n",
      "https://pmc.ncbi.nlm.nih.gov/articles/PMC7681068/\n",
      "and save as inst/extdata/lehallier2020_coefs.csv"
    )
  }

  required_somaids <- proteins$SOMAID
  coef_vec <- stats::setNames(proteins$Weight, proteins$SOMAID)

  present  <- required_somaids %in% names(data)
  missing  <- required_somaids[!present]

  if (any(!present)) {
    warning(sum(!present), " of ", length(required_somaids),
            " required SOMAmer entries not found")
  }

  n_samples <- nrow(data)
  prot_age <- numeric(n_samples)
  ids <- data[[id_col]]
  chron_age <- as.numeric(data[[age_col]])

  for (i in seq_len(n_samples)) {
    row <- data[i, ]
    pred <- intercept
    for (sid in required_somaids) {
      if (sid %in% names(row)) {
        val <- row[[sid]]
        if (!is.na(val) && is.numeric(val)) {
          if (log_transform && val > 0) {
            val <- log2(val)
          } else if (log_transform && val <= 0) {
            val <- NA_real_
          }
          if (!is.na(val)) {
            pred <- pred + coef_vec[sid] * val
          }
        }
      }
    }
    prot_age[i] <- pred
  }

  data.frame(
    id = ids,
    chronological_age = chron_age,
    proteomic_age = prot_age,
    age_acceleration = prot_age - chron_age,
    n_proteins_matched = sum(present),
    n_proteins_missing = length(missing),
    stringsAsFactors = FALSE
  )
}
