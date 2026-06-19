#' Lehallier et al. (2019) 373-protein LASSO proteomic clock
#'
#' Implements the clock from Lehallier B, et al. (2019).
#' Nature Medicine 25:1843-1850. LASSO regression, log10 + Z-scaling.
#' Input columns must be named by UniProt accession.
#'
#' @name lehallier2019
#' @keywords internal
NULL

# ---------- lehallier2019_proteins ----------

#' List Lehallier 2019 clock proteins
#' @return data.frame with SOMAID, Gene, UniProt, Weight
#' @export
lehallier2019_proteins <- function() {
  load_lehallier2019_coefs()
  .lehallier2019_cache$proteins
}

# ---------- compute_lehallier2019_age ----------

#' Compute proteomic age via Lehallier 2019 LASSO clock
#'
#' @param data data.frame with UniProt-named columns + Age + Sex
#' @param id_col Sample ID column name
#' @param age_col Chronological age column name
#' @param sex_col Sex column name (0=male, 1=female)
#' @param male_value Value coding for male
#' @return data.frame with proteomic_age, age_acceleration
#' @export
compute_lehallier2019_age <- function(data,
                                       id_col = "SampleID",
                                       age_col = "Age",
                                       sex_col = "Sex",
                                       male_value = 0) {

  if (!is.data.frame(data)) stop("'data' must be a data.frame")
  if (!id_col %in% names(data)) stop("id_col not found")
  if (!age_col %in% names(data)) stop("age_col not found")

  load_lehallier2019_coefs()
  up_lookup <- .lehallier2019_cache$uniprot_lookup
  weight_lk <- .lehallier2019_cache$weight_lookup
  intercept <- .lehallier2019_cache$intercept

  user_cols <- setdiff(names(data), c(id_col, age_col, sex_col))
  matched_up <- intersect(user_cols, names(up_lookup))
  unmatched   <- setdiff(names(up_lookup), user_cols)

  if (length(unmatched) > 0) {
    warning(length(unmatched), " clock proteins not found. Missing: ",
            paste(head(unmatched, 10), collapse = ", "),
            if (length(unmatched) > 10) "..." else "")
  }

  ids <- data[[id_col]]
  chron_age <- as.numeric(data[[age_col]])
  sex_num <- if (sex_col %in% names(data)) {
    as.numeric(data[[sex_col]] != male_value)
  } else {
    warning("Sex column not found. Assuming all female.")
    rep(1, nrow(data))
  }

  n <- nrow(data)
  prot_age <- numeric(n)

  for (i in seq_len(n)) {
    pred <- intercept
    row <- data[i, ]
    for (up in matched_up) {
      val <- row[[up]]
      if (!is.na(val) && is.numeric(val) && val > 0) {
        val <- log10(val)
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
    id = ids,
    chronological_age = chron_age,
    proteomic_age = prot_age,
    age_acceleration = age_accel,
    n_proteins_matched = length(matched_up),
    n_proteins_missing = length(unmatched),
    stringsAsFactors = FALSE
  )
}
