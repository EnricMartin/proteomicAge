#' Oh et al. (2023) conventional proteomic aging clock
#'
#' Oh HS, Rutledge J, et al. Organ aging signatures in the plasma proteome
#' track health and disease. Nature. 2023;624:164-172.
#' doi: 10.1038/s41586-023-06802-1
#'
#' Bagged LASSO ensemble (500 models) on SOMAscan v4 (4,778 proteins),
#' trained on 1,398 healthy participants from Knight-ADRC.
#'
#' @name oh2023
NULL

#' List Oh 2023 conventional clock proteins
#' @export
oh2023_conventional_proteins <- function() {
  load_oh2023_conventional_coefs()
  .oh2023_conventional_cache$proteins
}

#' Compute proteomic age via Oh 2023 conventional clock
#'
#' @param data data.frame with protein columns + Age column.
#' @param id_col Sample ID column name.
#' @param age_col Chronological age column name.
#' @param match_by How to match: "seqid_dot" (default), "uniprot", "gene", "seqid_sl".
#' @return data.frame with proteomic_age, age_acceleration
#' @export
compute_oh2023_conventional_age <- function(data,
                                             id_col = "SampleID",
                                             age_col = "Age",
                                             match_by = c("seqid_dot", "uniprot", "gene", "seqid_sl")) {

  match_by <- match.arg(match_by)
  if (!is.data.frame(data)) stop("'data' must be a data.frame")
  if (!id_col %in% names(data)) stop("id_col not found")
  if (!age_col %in% names(data)) stop("age_col not found")

  load_oh2023_conventional_coefs()
  lk_name <- paste0("lookup_", match_by)
  lookup <- .oh2023_conventional_cache[[lk_name]]
  weight_lk <- .oh2023_conventional_cache$lookup_Weight
  intercept <- .oh2023_conventional_cache$intercept

  user_cols <- setdiff(names(data), c(id_col, age_col, "Sex"))
  matched <- intersect(user_cols, names(lookup))
  matched_somaids <- unique(unname(lookup[matched]))
  matched <- matched[!duplicated(unname(lookup[matched]))]
  unmatched <- setdiff(names(weight_lk), matched_somaids)

  if (length(unmatched) > 0) {
    message(length(unmatched), " clock proteins not found")
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
