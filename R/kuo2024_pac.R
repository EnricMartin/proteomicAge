#' Kuo et al. (2024) PAC Olink proteomic aging clock
#'
#' Kuo CL, Chen Z, Liu P, et al. Proteomic aging clock (PAC) predicts
#' age-related outcomes in middle-aged and older adults. Aging Cell. 2024.
#'
#' PAC uses chronological age plus 128 Olink proteins to estimate a
#' mortality-calibrated proteomic age.
#'
#' @name kuo2024_pac
NULL

#' List Kuo 2024 PAC predictors
#'
#' @return data.frame with predictor names and coefficients.
#' @export
kuo2024_pac_proteins <- function() {
  load_kuo2024_pac_coefs()
  .kuo2024_pac_cache$coefs[.kuo2024_pac_cache$coefs$predictor != "age", ]
}

#' Compute PAC proteomic age via Kuo et al. (2024)
#'
#' @param data data.frame with Olink NPX protein columns plus age.
#' @param id_col Sample ID column name.
#' @param age_col Chronological age column name.
#' @param match_by Protein column naming convention: "gene" or "uniprot".
#' @param protein_map Optional data.frame mapping gene/assay names to UniProt
#'   IDs. If NULL, the built-in Olink map is used.
#' @return data.frame with PAC proteomic age and age acceleration.
#' @export
compute_kuo2024_pac_age <- function(data,
                                    id_col = "SampleID",
                                    age_col = "Age",
                                    match_by = c("gene", "uniprot"),
                                    protein_map = NULL) {
  match_by <- match.arg(match_by)
  if (!is.data.frame(data)) stop("'data' must be a data.frame")
  if (!id_col %in% names(data)) stop("id_col not found")
  if (!age_col %in% names(data)) stop("age_col not found")

  load_kuo2024_pac_coefs()
  predictors <- .kuo2024_pac_cache$coefs$predictor
  proteins <- .kuo2024_pac_cache$proteins
  input_lookup <- .olink_match_lookup(names(data), proteins, match_by, protein_map)
  missing <- setdiff(toupper(proteins), names(input_lookup))
  if (length(missing) > 0) {
    stop("Missing PAC predictors: ", paste(missing, collapse = ", "))
  }

  chron_age <- as.numeric(data[[age_col]])
  design <- matrix(NA_real_, nrow = nrow(data), ncol = length(predictors))
  colnames(design) <- predictors
  design[, "age"] <- chron_age
  for (protein in proteins) {
    design[, protein] <- as.numeric(data[[input_lookup[[toupper(protein)]]]])
  }

  weights <- .kuo2024_pac_cache$weights[predictors]
  params <- .kuo2024_pac_cache$params
  b_x <- as.numeric(design %*% weights)
  rate_new <- params[["rate"]] * exp(b_x)
  cdf_10_year <- 1 - exp(-(rate_new / params[["shape"]]) *
                           (exp(params[["shape"]] * 10) - 1))
  prot_age <- (1 / params[["beta_age"]]) *
    log(params[["shape0"]] * log(1 - cdf_10_year) /
          (params[["rate0"]] * (1 - exp(10 * params[["shape0"]]))))
  age_accel <- as.numeric(stats::residuals(stats::lm(prot_age ~ chron_age)))

  data.frame(
    id = data[[id_col]],
    chronological_age = chron_age,
    proteomic_age = prot_age,
    age_acceleration = age_accel,
    n_proteins_matched = length(proteins),
    n_proteins_missing = 0L,
    match_by = match_by,
    clock = "Kuo2024_PAC",
    stringsAsFactors = FALSE
  )
}
