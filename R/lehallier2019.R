#' Lehallier et al. (2019) 373-protein LASSO proteomic clock
#'
#' Lehallier B, et al. Nat Med 2019;25:1843-1850.
#'
#' @name lehallier2019
NULL

#' List Lehallier 2019 clock proteins
#' @export
lehallier2019_proteins <- function() {
  load_lehallier2019_coefs()
  .lehallier2019_cache$proteins
}

#' Compute proteomic age via Lehallier 2019 LASSO clock
#'
#' @param data data.frame with protein columns + Age column.
#' @param id_col Sample ID column name.
#' @param age_col Chronological age column name.
#' @param sex_col Sex column name (0=male, 1=female).
#' @param male_value Value coding for male.
#' @param match_by How to match: \code{"uniprot"} (default), \code{"gene"},
#'   \code{"seqid_sl"}, \code{"seqid_dot"}.
#' @return data.frame with proteomic_age, age_acceleration
#' @export
compute_lehallier2019_age <- function(data,
                                       id_col = "SampleID",
                                       age_col = "Age",
                                       sex_col = "Sex",
                                       male_value = 0,
                                       match_by = c("uniprot", "gene", "seqid_sl", "seqid_dot")) {

  match_by <- match.arg(match_by)
  if (!is.data.frame(data)) stop("'data' must be a data.frame")
  if (!id_col %in% names(data)) stop("id_col not found")
  if (!age_col %in% names(data)) stop("age_col not found")

  load_lehallier2019_coefs()
  lk_name <- paste0("lookup_", match_by)
  lookup <- .lehallier2019_cache[[lk_name]]
  weight_lk <- .lehallier2019_cache$lookup_Weight
  intercept <- .lehallier2019_cache$intercept

  user_cols <- setdiff(names(data), c(id_col, age_col, sex_col))
  matched <- intersect(user_cols, names(lookup))
  matched_somaids <- unique(unname(lookup[matched]))
  matched <- matched[!duplicated(unname(lookup[matched]))]
  unmatched <- setdiff(names(weight_lk), matched_somaids)

  if (length(unmatched) > 0) {
    warning(length(unmatched), " clock proteins not found. Match by: ", match_by)
  }

  sex_num <- if (sex_col %in% names(data)) {
    as.numeric(data[[sex_col]] != male_value)
  } else {
    rep(1, nrow(data))
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
        val <- log10(val)
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
