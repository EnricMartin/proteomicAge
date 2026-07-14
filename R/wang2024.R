#' Wang et al. (2024) ARIC midlife proteomic aging clock
#'
#' Wang S, Rao Z, Cao R, et al. Development, characterization, and
#' replication of proteomic aging clocks: Analysis of 2 population-based
#' cohorts. PLOS Medicine. 2024;21(9):e1004464.
#' doi: 10.1371/journal.pmed.1004464
#'
#' Elastic net (alpha=0.5, 10-fold CV), log2 transform.
#' 4,955 SOMAscan v4 aptamers → 788 selected.
#' Midlife: 11,761 participants aged 46-70 (ARIC Visit 2, 1990-92).
#'
#' @name wang2024
NULL

#' List Wang 2024 ARIC midlife clock proteins
#' @export
wang2024_aric_proteins <- function() {
  load_wang2024_aric_coefs()
  .wang2024_aric_cache$proteins
}

#' Compute proteomic age via Wang 2024 ARIC midlife clock
#'
#' @param data data.frame with protein columns + Age.
#' @param id_col Sample ID column name.
#' @param age_col Chronological age column name.
#' @param match_by How to match: "seqid_dot" (default), "uniprot", "gene", "seqid_sl".
#' @return data.frame with proteomic_age, age_acceleration
#' @export
compute_wang2024_aric_age <- function(data,
                                       id_col = "SampleID",
                                       age_col = "Age",
                                       match_by = c("seqid_dot", "uniprot", "gene", "seqid_sl")) {

  match_by <- match.arg(match_by)
  if (!is.data.frame(data)) stop("'data' must be a data.frame")
  if (!id_col %in% names(data)) stop("id_col not found")
  if (!age_col %in% names(data)) stop("age_col not found")

  load_wang2024_aric_coefs()
  lk_name <- paste0("lookup_", match_by)
  lookup <- .wang2024_aric_cache[[lk_name]]
  weight_lk <- .wang2024_aric_cache$lookup_Weight
  intercept <- .wang2024_aric_cache$intercept

  user_cols <- setdiff(names(data), c(id_col, age_col, "Sex"))
  matched <- intersect(user_cols, names(lookup))
  matched_somaids <- unique(unname(lookup[matched]))
  matched <- matched[!duplicated(unname(lookup[matched]))]
  unmatched <- setdiff(names(weight_lk), matched_somaids)

  if (length(unmatched) > 0) message(length(unmatched), " clock proteins not found")

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
        val <- log2(val)
        sid <- lookup[col]
        pred <- pred + weight_lk[sid] * val
      }
    }
    prot_age[i] <- pred
  }

  fit <- stats::lm(prot_age ~ chron_age)
  age_accel <- as.numeric(stats::residuals(fit))

  data.frame(
    id = ids, chronological_age = chron_age, proteomic_age = prot_age,
    age_acceleration = age_accel, n_proteins_matched = length(matched),
    n_proteins_missing = length(unmatched), match_by = match_by,
    stringsAsFactors = FALSE
  )
}
