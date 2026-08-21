#' Goeminne et al. organ-specific Olink aging models
#'
#' Plasma protein-based organ-specific aging and mortality models from
#' Goeminne et al. Full Olink Explore 3072 coefficients are included for
#' academic/non-commercial use under the source repository license.
#'
#' @name goeminne2025
NULL

#' List Goeminne 2025 organ aging model proteins
#'
#' @return data.frame with organ, model type, fold, protein, and weight.
#' @export
goeminne2025_organaging_proteins <- function() {
  load_goeminne2025_coefs()
  .goeminne2025_cache$coefs[.goeminne2025_cache$coefs$protein != "Intercept", ]
}

#' Compute Goeminne 2025 organ-specific Olink aging models
#'
#' @param data data.frame with Olink NPX protein columns.
#' @param id_col Sample ID column name.
#' @param age_col Chronological age column name.
#' @param organs Character vector of organs to compute, or "all".
#' @param model_type "chronological" or "mortality".
#' @param fold Model fold to use; defaults to 1, as suggested by authors.
#' @param match_by Protein column naming convention: "auto", "gene", or
#'   "uniprot".
#' @param protein_map Optional data.frame mapping gene/assay names to UniProt
#'   IDs. If NULL, the built-in Olink map is used.
#' @param group Optional grouping vector or column name for bundled QC plots.
#' @param sample_data Optional data.frame containing group metadata and IDs.
#' @param return_list If TRUE, return predictions, QC, plots, and group comparison
#'   as a list. This is automatically enabled when `group` is supplied.
#' @return Long data.frame with one row per sample and organ.
#' @export
compute_goeminne2025_organ_age <- function(data,
                                           id_col = "SampleID",
                                           age_col = "Age",
                                           organs = "all",
                                           model_type = c("chronological", "mortality"),
                                           fold = 1,
                                           match_by = c("auto", "gene", "uniprot"),
                                           protein_map = NULL,
                                           group = NULL,
                                           sample_data = NULL,
                                           return_list = !is.null(group)) {
  model_type <- match.arg(model_type)
  match_by <- match.arg(match_by)
  if (!is.data.frame(data)) stop("'data' must be a data.frame")
  if (!id_col %in% names(data)) stop("id_col not found")
  if (!age_col %in% names(data)) stop("age_col not found")
  if (!fold %in% 1:5) stop("fold must be an integer from 1 to 5")

  load_goeminne2025_coefs()
  all_organs <- .goeminne2025_cache$organs
  if (identical(organs, "all")) organs <- all_organs
  missing_organs <- setdiff(organs, all_organs)
  if (length(missing_organs) > 0) {
    stop("Unknown Goeminne organ model: ", paste(missing_organs, collapse = ", "))
  }

  coefs <- .goeminne2025_cache$coefs
  coefs <- coefs[coefs$model_type == model_type &
                   coefs$fold == fold &
                   coefs$organ %in% organs, ]
  chron_age <- as.numeric(data[[age_col]])
  out <- vector("list", length(organs))

  for (i in seq_along(organs)) {
    organ_coefs <- coefs[coefs$organ == organs[[i]], ]
    intercept <- sum(organ_coefs$weight[organ_coefs$protein == "Intercept"])
    protein_coefs <- organ_coefs[organ_coefs$protein != "Intercept", ]
    input_lookup <- .olink_match_lookup(
      names(data), protein_coefs$protein, match_by, protein_map
    )
    matched_by <- attr(input_lookup, "match_by")
    matched <- protein_coefs$protein[toupper(protein_coefs$protein) %in% names(input_lookup)]
    missing <- setdiff(protein_coefs$protein, matched)
    if (length(matched) == 0) {
      stop("No proteins matched for Goeminne organ model '", organs[[i]],
           "'. Check Olink protein column names or set match_by explicitly.")
    }
    score <- rep(intercept, nrow(data))
    for (protein in matched) {
      score <- score + as.numeric(data[[input_lookup[[toupper(protein)]]]]) *
        protein_coefs$weight[match(protein, protein_coefs$protein)]
    }
    if (model_type == "chronological") {
      out[[i]] <- data.frame(
        id = data[[id_col]],
        chronological_age = chron_age,
        organ = organs[[i]],
        proteomic_age = score,
        age_acceleration = as.numeric(stats::residuals(stats::lm(score ~ chron_age))),
        mortality_score = NA_real_,
        n_proteins_matched = length(matched),
        n_proteins_missing = length(missing),
        match_by = matched_by,
        clock = "Goeminne2025_organAging",
        stringsAsFactors = FALSE
      )
    } else {
      out[[i]] <- data.frame(
        id = data[[id_col]],
        chronological_age = chron_age,
        organ = organs[[i]],
        proteomic_age = NA_real_,
        age_acceleration = NA_real_,
        mortality_score = score,
        n_proteins_matched = length(matched),
        n_proteins_missing = length(missing),
        match_by = matched_by,
        clock = "Goeminne2025_organAging_mortality",
        stringsAsFactors = FALSE
      )
    }
  }

  result <- do.call(rbind, out)
  .clock_result_bundle(
    result, group, sample_data, data, id_col, return_list,
    "Goeminne2025_organAging"
  )
}
