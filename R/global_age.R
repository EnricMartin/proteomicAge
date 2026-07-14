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
#' @return data.frame with global proteomic_age, age_acceleration, and component
#'   clock predictions.
#' @export
compute_global_age <- function(data,
                               id_col = "SampleID",
                               age_col = "Age",
                               sex_col = NULL,
                               male_value = 0,
                               match_by = c("seqid_dot", "uniprot", "gene", "seqid_sl")) {

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
  colnames(age_mat) <- paste0(names(clock_results), "_age")
  colnames(gap_mat) <- paste0(names(clock_results), "_age_acceleration")

  base <- clock_results[[1]]
  out <- data.frame(
    id = base$id,
    chronological_age = base$chronological_age,
    proteomic_age = rowMeans(age_mat, na.rm = TRUE),
    age_acceleration = rowMeans(gap_mat, na.rm = TRUE),
    n_clocks = ncol(age_mat),
    clocks_used = paste(names(clock_results), collapse = ";"),
    match_by = match_by,
    stringsAsFactors = FALSE
  )

  cbind(out, as.data.frame(age_mat), as.data.frame(gap_mat))
}
