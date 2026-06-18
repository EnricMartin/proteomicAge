# ---------------------------------------------------------------------------
# SOMAscan data preprocessing utilities
# ---------------------------------------------------------------------------

#' Preprocess SOMAscan proteomic data
#'
#' Performs standard preprocessing on SOMAscan RFU data:
#' log2-transformation, outlier handling, and missing value reporting.
#' Based on the preprocessing described in Tanaka et al. (2018):
#' "Each protein was log-transformed, and outliers +/-4 SD were removed."
#'
#' @param data A data.frame of SOMAscan RFU values (samples x proteins).
#'   Protein columns should be named with SOMAIDs (SeqIds, e.g. "SL003869").
#' @param protein_cols Character vector of protein column names. If NULL
#'   (default), all columns matching the SOMAID pattern (starting with "SL")
#'   are treated as protein columns.
#' @param log_transform Logical. If TRUE (default), log2-transform values.
#' @param handle_outliers Logical. If TRUE (default), winsorize values
#'   beyond 4 standard deviations from the mean.
#' @param report_missingness Logical. If TRUE (default), print a summary
#'   of missing values per protein.
#'
#' @return A data.frame with log-transformed (and optionally winsorized)
#'   protein values. Non-protein columns are preserved.
#' @export
#'
#' @examples
#' \dontrun{
#' raw <- read.csv("somascan_raw.csv")
#' processed <- preprocess_somascan(raw)
#' }
preprocess_somascan <- function(data,
                                 protein_cols = NULL,
                                 log_transform = TRUE,
                                 handle_outliers = TRUE,
                                 report_missingness = TRUE) {

  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame")
  }

  # Identify protein columns
  if (is.null(protein_cols)) {
    protein_cols <- grep("^SL[0-9]", names(data), value = TRUE)
    if (length(protein_cols) == 0) {
      stop("No protein columns found matching SOMAID pattern 'SL[0-9]...'.",
           " Specify protein_cols manually.")
    }
    message("Auto-detected ", length(protein_cols), " protein columns")
  }

  # Check that protein columns exist
  missing_cols <- setdiff(protein_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Protein columns not found: ", paste(missing_cols, collapse = ", "))
  }

  # Log-transform
  if (log_transform) {
    for (col in protein_cols) {
      vals <- data[[col]]
      if (is.numeric(vals)) {
        # Handle zeros/negatives: set to NA before log
        bad <- vals <= 0 & !is.na(vals)
        if (any(bad)) {
          warning("Found ", sum(bad), " non-positive values in '", col,
                  "'. Setting to NA before log2 transform.")
          vals[bad] <- NA_real_
        }
        data[[col]] <- log2(vals)
      }
    }
  }

  # Outlier handling (winsorize at +/- 4 SD)
  if (handle_outliers) {
    for (col in protein_cols) {
      vals <- data[[col]]
      if (is.numeric(vals)) {
        valid <- vals[!is.na(vals)]
        if (length(valid) > 3) {
          m <- mean(valid, na.rm = TRUE)
          s <- stats::sd(valid, na.rm = TRUE)
          lo <- m - 4 * s
          hi <- m + 4 * s
          n_lo <- sum(vals < lo, na.rm = TRUE)
          n_hi <- sum(vals > hi, na.rm = TRUE)
          if (n_lo + n_hi > 0) {
            vals[vals < lo & !is.na(vals)] <- lo
            vals[vals > hi & !is.na(vals)] <- hi
            data[[col]] <- vals
          }
        }
      }
    }
  }

  # Report missingness
  if (report_missingness) {
    na_counts <- vapply(data[protein_cols], function(x) sum(is.na(x)), integer(1))
    na_pct    <- round(na_counts / nrow(data) * 100, 2)
    summary_df <- data.frame(
      Protein  = protein_cols,
      N_missing = na_counts,
      Pct_missing = na_pct,
      stringsAsFactors = FALSE
    )
    summary_df <- summary_df[order(-summary_df$N_missing), ]
    cat("\n--- Missing Value Summary ---\n")
    cat("Total samples:", nrow(data), "\n")
    cat("Proteins with any missing values:", sum(na_counts > 0), "/", length(protein_cols), "\n")
    if (sum(na_counts > 0) > 0) {
      cat("\nTop proteins by missingness:\n")
      print(head(summary_df[summary_df$N_missing > 0, ], 10))
    }
    cat("-----------------------------\n\n")
  }

  return(data)
}


#' Validate SOMAscan input data for proteomic clock
#'
#' Checks whether a data.frame is suitable for proteomic clock computation.
#' Reports: required columns, protein coverage, missing values, data types.
#'
#' @param data A data.frame of SOMAscan proteomic data.
#' @param id_col Name of the identifier column.
#' @param age_col Name of the age column.
#' @param required_proteins Character vector of required SOMAIDs.
#'   Default: Tanaka 2018 76 proteins.
#'
#' @return Invisibly returns a list with validation results.
#' @export
validate_somascan_input <- function(data,
                                     id_col = "SampleID",
                                     age_col = "Age",
                                     required_proteins = NULL) {

  results <- list(
    valid = TRUE,
    warnings = character(0),
    errors = character(0),
    protein_coverage = NULL
  )

  # Check basic structure
  if (!is.data.frame(data)) {
    results$errors <- c(results$errors, "Input is not a data.frame")
    results$valid <- FALSE
    cat("ERROR: Input is not a data.frame\n")
    return(invisible(results))
  }

  # Check required columns
  if (!id_col %in% names(data)) {
    results$errors <- c(results$errors, paste("ID column '", id_col, "' not found"))
    results$valid <- FALSE
  }
  if (!age_col %in% names(data)) {
    results$errors <- c(results$errors, paste("Age column '", age_col, "' not found"))
    results$valid <- FALSE
  }

  # Check age values
  if (age_col %in% names(data)) {
    ages <- as.numeric(data[[age_col]])
    if (any(is.na(ages))) {
      results$warnings <- c(results$warnings,
        paste(sum(is.na(ages)), "missing age values"))
    }
    if (any(ages < 0 | ages > 120, na.rm = TRUE)) {
      results$warnings <- c(results$warnings,
        "Some age values are outside plausible range (0-120)")
    }
  }

  # Check protein coverage
  if (is.null(required_proteins)) {
    load_tanaka2018_coefs()
    required_proteins <- .tanaka2018_cache$proteins$SOMAID
  }

  present  <- required_proteins %in% names(data)
  coverage <- sum(present) / length(required_proteins) * 100
  results$protein_coverage <- coverage

  if (coverage < 100) {
    missing_prots <- required_proteins[!present]
    results$warnings <- c(results$warnings,
      paste0(sprintf("Protein coverage: %.1f%% (%d/%d). ", coverage,
                     sum(present), length(required_proteins)),
             "Missing: ", paste(missing_prots, collapse = ", ")))
  }

  # Print report
  cat("=== SOMAscan Input Validation ===\n")
  cat("Samples:", nrow(data), "\n")
  cat("Total columns:", ncol(data), "\n")

  if (length(results$errors) > 0) {
    cat("\nERRORS:\n")
    for (e in results$errors) cat("  -", e, "\n")
  }

  if (length(results$warnings) > 0) {
    cat("\nWARNINGS:\n")
    for (w in results$warnings) cat("  -", w, "\n")
  }

  if (results$valid) {
    cat(sprintf("\nProtein coverage: %.1f%% (%d/%d required proteins found)\n",
                coverage, sum(present), length(required_proteins)))
    cat("Overall: VALID\n")
  } else {
    cat("\nOverall: INVALID\n")
  }
  cat("==================================\n")

  invisible(results)
}
