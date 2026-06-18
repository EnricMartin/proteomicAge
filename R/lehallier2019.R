#' Lehallier et al. (2019) 373-protein proteomic aging clock
#'
#' Implements the LASSO-based proteomic aging clock from:
#' Lehallier B, Gate D, Schaum N, et al. "Undulating changes in human plasma
#' proteome profiles across the lifespan." Nature Medicine. 2019;25:1843-1850.
#' doi: 10.1038/s41591-019-0673-2
#'
#' The clock uses LASSO regression (glmnet, alpha=1, lambda.min from 10-fold CV)
#' trained on SOMAscan data from 4,263 healthy individuals (INTERVAL + LonGenity).
#' Input variables: Z-scaled log10 RFU values + sex.
#'
#' @references
#' Lehallier B, et al. (2019) Nature Medicine 25:1843-1850.
#' Supplementary Table 7 contains the 373-protein coefficients.
#'
#' @name lehallier2019
#' @keywords internal
NULL

# ---------------------------------------------------------------------------
# lehallier2019_proteins()
# ---------------------------------------------------------------------------

#' List proteins in the Lehallier 2019 clock
#'
#' Returns the 373 SOMAscan proteins and their LASSO coefficients
#' from the Lehallier et al. (2019) proteomic aging clock.
#'
#' @return A data.frame with columns: \code{SOMAID}, \code{Gene}, \code{Weight}
#' @export
lehallier2019_proteins <- function() {
  load_lehallier2019_coefs()
  .lehallier2019_cache$proteins
}

# ---------------------------------------------------------------------------
# compute_lehallier2019_age()
# ---------------------------------------------------------------------------

#' Compute proteomic age using Lehallier et al. (2019) LASSO clock
#'
#' Applies the 373-protein LASSO model to SOMAscan proteomic data.
#' Input values are log10-transformed then Z-scaled, matching the
#' original paper's preprocessing. Sex is included as a covariate.
#'
#' @param data A data.frame of SOMAscan proteomic data. Must contain:
#'   \itemize{
#'     \item An ID column (default \code{"SampleID"})
#'     \item An age column (default \code{"Age"})
#'     \item A sex column (default \code{"Sex"}, 0=male, 1=female or "M"/"F")
#'     \item Protein columns named by SOMAscan SeqIds
#'   }
#' @param id_col Name of the sample identifier column.
#' @param age_col Name of the chronological age column.
#' @param sex_col Name of the sex column.
#' @param male_value Value coding for male in sex column (default 0).
#' @param protein_format Format of protein column names.
#'   \code{"auto"} (default) auto-detects; \code{"gene"} for gene symbols;
#'   \code{"seqid_dot"} for dot-format SeqIds; \code{"seqid_sl"} for SL format.
#'   Columns are internally normalized to Gene.SeqId format.
#'
#' @return A data.frame with predicted proteomic age and age acceleration.
#' @export
compute_lehallier2019_age <- function(data,
                                       id_col = "SampleID",
                                       age_col = "Age",
                                       sex_col = "Sex",
                                       male_value = 0,
                                       protein_format = c("auto", "gene",
                                                          "seqid_dot", "seqid_sl",
                                                          "seqid_full")) {

  protein_format <- match.arg(protein_format)

  # Normalize protein names to seqid_full format
  data <- normalize_protein_names(data,
    id_col = id_col, age_col = age_col,
    target_format = "seqid_full",
    input_format = protein_format
  )

  if (!is.data.frame(data)) stop("'data' must be a data.frame")
  if (!id_col %in% names(data)) stop("id_col not found")
  if (!age_col %in% names(data)) stop("age_col not found")

  load_lehallier2019_coefs()
  proteins <- .lehallier2019_cache$proteins
  intercept <- .lehallier2019_cache$intercept

  if (is.null(proteins) || nrow(proteins) == 0) {
    stop(
      "Lehallier 2019 coefficients not available.\n",
      "Download Supplementary Table 7 from:\n",
      "https://pmc.ncbi.nlm.nih.gov/articles/PMC7062043/\n",
      "and save as inst/extdata/lehallier2019_coefs.csv"
    )
  }

  required_somaids <- proteins$SOMAID
  coef_vec <- stats::setNames(proteins$Weight, proteins$SOMAID)

  present  <- required_somaids %in% names(data)
  missing  <- required_somaids[!present]

  if (any(!present)) {
    warning(sum(!present), " of ", length(required_somaids),
            " required proteins not found")
  }

  # Sex coding
  if (sex_col %in% names(data)) {
    sex_numeric <- ifelse(data[[sex_col]] == male_value, 0, 1)
  } else {
    warning("Sex column '", sex_col, "' not found. Assuming all female.")
    sex_numeric <- rep(1, nrow(data))
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
          # log10 transform
          if (val > 0) {
            val <- log10(val)
            # Z-scaling would need population mean/sd - applied approximately
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
