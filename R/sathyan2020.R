#' Sathyan et al. (2020) proteomic aging clock
#'
#' Sathyan S, Ayers E, Gao T, et al. Plasma proteomic profile of age,
#' health span, and all-cause mortality in older adults.
#' Aging Cell. 2020;19(11):e13250. doi: 10.1111/acel.13250
#'
#' Elastic net regression (alpha=0.5, 10-fold CV) on SOMAscan v4
#' (4,265 SOMAmers), 1,025 older adults from LonGenity cohort.
#' Model: log(SOMAmer) ~ age + gender + cohort.
#'
#' @name sathyan2020
NULL

#' List Sathyan 2020 clock proteins
#' @export
sathyan2020_proteins <- function() {
  load_sathyan2020_coefs()
  .sathyan2020_cache$proteins
}

#' Compute proteomic age via Sathyan 2020 clock
#'
#' @param data data.frame with protein columns + Age column.
#' @param id_col Sample ID column name.
#' @param age_col Chronological age column name.
#' @param sex_col Sex column name.
#' @param match_by How to match: "uniprot" (default), "gene", "seqid_sl", "seqid_dot".
#' @return data.frame with proteomic_age, age_acceleration
#' @export
compute_sathyan2020_age <- function(data,
                                     id_col = "SampleID",
                                     age_col = "Age",
                                     sex_col = "Sex",
                                     match_by = c("uniprot", "gene", "seqid_sl", "seqid_dot")) {

  match_by <- match.arg(match_by)
  if (!is.data.frame(data)) stop("'data' must be a data.frame")
  if (!id_col %in% names(data)) stop("id_col not found")
  if (!age_col %in% names(data)) stop("age_col not found")

  load_sathyan2020_coefs()
  lk_name <- paste0("lookup_", match_by)
  lookup <- .sathyan2020_cache[[lk_name]]
  weight_lk <- .sathyan2020_cache$lookup_Weight
  intercept <- .sathyan2020_cache$intercept

  if (is.null(lookup) || length(lookup) == 0) {
    stop("Sathyan 2020 coefficients not loaded. Download supplementary tables from\n",
         "https://onlinelibrary.wiley.com/doi/10.1111/acel.13250\n",
         "and save to inst/extdata/sathyan2020_coefs.csv")
  }

  user_cols <- setdiff(names(data), c(id_col, age_col, sex_col))
  matched <- intersect(user_cols, names(lookup))
  matched_somaids <- unique(unname(lookup[matched]))
  matched <- matched[!duplicated(unname(lookup[matched]))]
  unmatched <- setdiff(names(weight_lk), matched_somaids)

  if (length(unmatched) > 0) {
    warning(length(unmatched), " clock proteins not found. Match by: ", match_by)
  }

  ids <- data[[id_col]]
  chron_age <- as.numeric(data[[age_col]])
  n <- nrow(data)
  prot_age <- numeric(n)

  for (i in seq_len(n)) {
    pred <- intercept
    row <- data[i, ]
    for (col in matched) {
      val <- row[[col]]
      if (!is.na(val) && is.numeric(val) && val > 0) {
        val <- log(val)
        sid <- lookup[col]
        pred <- pred + weight_lk[sid] * val
      }
    }
    prot_age[i] <- pred
  }

  fit <- stats::lm(prot_age ~ chron_age)
  age_accel <- as.numeric(stats::residuals(fit))

  data.frame(
    id = ids,
    chronological_age = chron_age,
    proteomic_age = prot_age,
    age_acceleration = age_accel,
    n_proteins_matched = length(matched),
    n_proteins_missing = length(unmatched),
    match_by = match_by,
    stringsAsFactors = FALSE
  )
}
