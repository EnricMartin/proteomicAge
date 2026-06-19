#' Tanaka et al. (2018) 76-protein elastic net proteomic clock
#'
#' Tanaka T, et al. Aging Cell 2018;17(5):e12799.
#' Coefficients from eLife (2020) Supplementary File 1D.
#'
#' @name tanaka2018
NULL

#' List Tanaka 2018 clock proteins
#' @export
tanaka2018_proteins <- function() {
  load_tanaka2018_coefs()
  .tanaka2018_cache$proteins
}

#' Compute proteomic age via Tanaka 2018 clock
#'
#' @param data data.frame with protein columns + Age column.
#'   Protein columns can use gene symbols, seq.xxx.xx, SLxxxxxx, or UniProt.
#' @param id_col Sample ID column name.
#' @param age_col Chronological age column name.
#' @param log_transform Apply log2 transform (default TRUE).
#' @param match_by How to match user columns to clock proteins:
#'   \code{"uniprot"} (default), \code{"gene"}, \code{"seqid_sl"}, \code{"seqid_dot"}.
#' @return data.frame with proteomic_age, age_acceleration, etc.
#' @export
compute_tanaka2018_age <- function(data,
                                    id_col = "SampleID",
                                    age_col = "Age",
                                    log_transform = TRUE,
                                    match_by = c("uniprot", "gene", "seqid_sl", "seqid_dot")) {

  match_by <- match.arg(match_by)
  if (!is.data.frame(data)) stop("'data' must be a data.frame")
  if (!id_col %in% names(data)) stop("id_col not found")
  if (!age_col %in% names(data)) stop("age_col not found")

  load_tanaka2018_coefs()
  lk_name <- paste0("lookup_", match_by)
  lookup <- .tanaka2018_cache[[lk_name]]
  weight_lk <- .tanaka2018_cache$lookup_Weight
  intercept <- .tanaka2018_cache$intercept

  user_cols <- setdiff(names(data), c(id_col, age_col, "Sex"))
  matched <- intersect(user_cols, names(lookup))
  unmatched <- setdiff(names(lookup), user_cols)

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
      if (!is.na(val) && is.numeric(val)) {
        if (log_transform && val > 0) {
          val <- log2(val)
        } else if (log_transform && val <= 0) {
          next
        }
        sid <- lookup[col]
        pred <- pred + weight_lk[sid] * val
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
    n_proteins_matched = length(matched),
    n_proteins_missing = length(unmatched),
    match_by          = match_by,
    stringsAsFactors   = FALSE
  )
}
