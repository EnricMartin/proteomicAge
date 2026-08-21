#' Global proteomic age ensemble
#'
#' Computes Global proteomic age by averaging predicted ages and age gaps from
#' the five conventional proteomic clocks: Tanaka 2018, Lehallier 2019,
#' Sathyan 2020, Oh 2023 conventional, and Wang 2024 ARIC.
#'
#' @param data data.frame with protein columns + Age column.
#' @param id_col Sample ID column name.
#' @param age_col Chronological age column name.
#' @param sex_col Optional sex column name. If NULL, uses "Sex_F" when present,
#'   otherwise "Sex" when present.
#' @param male_value Value coding for male in the sex column.
#' @param match_by How to match proteins: "seqid_dot" (default), "uniprot",
#'   "gene", or "seqid_sl".
#' @param group Optional grouping vector or column name for bundled QC plots.
#' @param sample_data Optional data.frame containing group metadata and IDs.
#' @param return_list If TRUE, return predictions, QC, plots, and group comparison
#'   as a list. This is automatically enabled when `group` is supplied.
#' @return data.frame with global proteomic_age, age_acceleration, and component
#'   clock predictions.
#' @export
compute_global_age <- function(data,
                               id_col = "SampleID",
                               age_col = "Age",
                               sex_col = NULL,
                               male_value = 0,
                               match_by = c("seqid_dot", "uniprot", "gene", "seqid_sl"),
                               group = NULL,
                               sample_data = NULL,
                               return_list = !is.null(group)) {

  match_by <- match.arg(match_by)
  if (!is.data.frame(data)) stop("'data' must be a data.frame")
  if (!id_col %in% names(data)) stop("id_col not found")
  if (!age_col %in% names(data)) stop("age_col not found")

  if (is.null(sex_col)) {
    sex_col <- if ("Sex_F" %in% names(data)) {
      "Sex_F"
    } else if ("Sex" %in% names(data)) {
      "Sex"
    } else {
      "__missing_sex__"
    }
  }

  clock_results <- list(
    tanaka2018 = compute_tanaka2018_age(
      data, id_col = id_col, age_col = age_col, match_by = match_by
    ),
    lehallier2019 = compute_lehallier2019_age(
      data, id_col = id_col, age_col = age_col, sex_col = sex_col,
      male_value = male_value, match_by = match_by
    ),
    sathyan2020 = compute_sathyan2020_age(
      data, id_col = id_col, age_col = age_col, sex_col = sex_col,
      match_by = match_by
    ),
    oh2023_conventional = compute_oh2023_conventional_age(
      data, id_col = id_col, age_col = age_col, sex_col = sex_col,
      match_by = match_by
    ),
    wang2024_aric = compute_wang2024_aric_age(
      data, id_col = id_col, age_col = age_col, match_by = match_by
    )
  )

  age_mat <- do.call(cbind, lapply(clock_results, `[[`, "proteomic_age"))
  gap_mat <- do.call(cbind, lapply(clock_results, `[[`, "age_acceleration"))
  matched <- vapply(clock_results, function(x) x$n_proteins_matched[1], numeric(1))
  missing <- vapply(clock_results, function(x) x$n_proteins_missing[1], numeric(1))
  colnames(age_mat) <- paste0(names(clock_results), "_age")
  colnames(gap_mat) <- paste0(names(clock_results), "_age_acceleration")
  protein_qc <- as.data.frame(
    as.list(c(matched, missing)),
    stringsAsFactors = FALSE
  )
  names(protein_qc) <- c(
    paste0(names(clock_results), "_n_proteins_matched"),
    paste0(names(clock_results), "_n_proteins_missing")
  )

  base <- clock_results[[1]]
  out <- data.frame(
    id = base$id,
    chronological_age = base$chronological_age,
    proteomic_age = rowMeans(age_mat, na.rm = TRUE),
    age_acceleration = rowMeans(gap_mat, na.rm = TRUE),
    n_proteins_matched = sum(matched, na.rm = TRUE),
    n_proteins_missing = sum(missing, na.rm = TRUE),
    n_clocks = ncol(age_mat),
    clocks_used = paste(names(clock_results), collapse = ";"),
    match_by = match_by,
    stringsAsFactors = FALSE
  )

  result <- cbind(
    out,
    as.data.frame(age_mat),
    as.data.frame(gap_mat),
    protein_qc[rep(1, nrow(out)), , drop = FALSE]
  )
  .clock_result_bundle(
    result, group, sample_data, data, id_col, return_list, NULL
  )
}
