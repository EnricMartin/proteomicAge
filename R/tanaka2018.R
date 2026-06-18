#' Tanaka et al. (2018) 76-protein proteomic aging clock
#'
#' Implements the proteomic aging clock from:
#' Tanaka T, et al. "Plasma proteomic signature of age in healthy humans."
#' Aging Cell. 2018 Oct;17(5):e12799. doi: 10.1111/acel.12799
#'
#' The clock uses elastic net regression (alpha=0.5, lambda=0.8767859)
#' trained on SOMAscan 1.3K data from 240 healthy adults (22-93 years).
#'
#' @references
#' Tanaka T, Basisty N, Fantoni G, et al. (2020). Plasma proteomic
#' biomarker signature of age predicts health and life span. eLife 9:e61073.
#' Supplementary File 1D contains the 76-protein signature weights.
#'
#' @name tanaka2018
#' @keywords internal
"_PACKAGE"

# ---------------------------------------------------------------------------
# tanaka2018_proteins()
# ---------------------------------------------------------------------------

#' List proteins in the Tanaka 2018 clock
#'
#' Returns the 76 SOMAscan proteins (SeqIds and gene symbols) and their
#' elastic net coefficients used in the Tanaka et al. (2018) proteomic
#' aging clock.
#'
#' @return A data.frame with columns: \code{SOMAID}, \code{Gene}, \code{Weight}
#' @export
#'
#' @examples
#' head(tanaka2018_proteins())
tanaka2018_proteins <- function() {
  load_tanaka2018_coefs()
  .tanaka2018_cache$proteins
}

# ---------------------------------------------------------------------------
# compute_tanaka2018_age()
# ---------------------------------------------------------------------------

#' Compute proteomic age using Tanaka et al. (2018) clock
#'
#' Applies the 76-protein elastic net model to SOMAscan proteomic data
#' to compute predicted biological age (PROage) and age acceleration
#' (PROaccel = residuals of predicted age on chronological age).
#'
#' The input data should be SOMAscan RFU values. The function will
#' log2-transform the values internally (matching the original paper's
#' pre-processing: "Each protein was log-transformed").
#'
#' @param data A data.frame of SOMAscan proteomic data. Must contain:
#'   \itemize{
#'     \item An ID column (user-specified, default \code{"SampleID"})
#'     \item An age column (user-specified, default \code{"Age"})
#'     \item Columns matching the 76 SOMAIDs in the Tanaka clock.
#'           Missing proteins will trigger a warning and be treated as NA.
#'   }
#' @param id_col Name of the sample identifier column.
#' @param age_col Name of the chronological age column.
#' @param log_transform Logical. If TRUE (default), log2-transform
#'   protein values before computing. Set to FALSE if data is already
#'   log-transformed.
#' @param protein_format Format of protein column names in the data.
#'   \code{"auto"} (default) auto-detects; \code{"gene"} for gene symbols;
#'   \code{"seqid_dot"} for dot-format SeqIds (e.g. 4374.45.2);
#'   \code{"seqid_sl"} for SL-format SeqIds (e.g. SL003869).
#'   Columns are internally normalized to SL format for the Tanaka clock.
#'
#' @return A data.frame with columns:
#'   \itemize{
#'     \item \code{id} — sample identifier
#'     \item \code{chronological_age} — input age
#'     \item \code{proteomic_age} — predicted biological age (PROage)
#'     \item \code{age_acceleration} — PROaccel (proteomic_age - chronological_age)
#'     \item \code{n_proteins_matched} — number of the 76 required proteins found
#'     \item \code{n_proteins_missing} — number of required proteins missing
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' # Read SOMAscan data (RFU values)
#' data <- read.csv("my_somascan_data.csv")
#'
#' # Compute proteomic ages
#' results <- compute_tanaka2018_age(data,
#'   id_col = "SampleID",
#'   age_col = "ChronologicalAge"
#' )
#' head(results)
#' }
compute_tanaka2018_age <- function(data,
                                    id_col = "SampleID",
                                    age_col = "Age",
                                    log_transform = TRUE,
                                    protein_format = c("auto", "gene",
                                                       "seqid_dot", "seqid_sl")) {

  protein_format <- match.arg(protein_format)

  # --- Normalize protein names to SL format ---
  data <- normalize_protein_names(data,
    id_col = id_col, age_col = age_col,
    target_format = "seqid_sl",
    input_format = protein_format
  )
  # --- Validate input ---
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame")
  }
  if (!id_col %in% names(data)) {
    stop("id_col '", id_col, "' not found in data")
  }
  if (!age_col %in% names(data)) {
    stop("age_col '", age_col, "' not found in data")
  }

  # --- Load coefficients ---
  load_tanaka2018_coefs()
  proteins <- .tanaka2018_cache$proteins
  intercept <- .tanaka2018_cache$intercept

  required_somaids <- proteins$SOMAID
  required_genes    <- proteins$Gene

  # --- Match SOMAIDs in data ---
  present  <- required_somaids %in% names(data)
  missing  <- required_somaids[!present]

  if (any(!present)) {
    warning(
      sum(!present), " of ", length(required_somaids),
      " required SOMAIDs not found in data. Missing proteins: ",
      paste(missing, collapse = ", ")
    )
  }

  # --- Build coefficient vector ---
  coef_vec <- stats::setNames(proteins$Weight, proteins$SOMAID)

  # --- Compute proteomic age for each sample ---
  ids <- data[[id_col]]
  chron_age <- as.numeric(data[[age_col]])

  n_samples <- nrow(data)
  prot_age <- numeric(n_samples)

  for (i in seq_len(n_samples)) {
    row <- data[i, ]
    pred <- intercept  # start with intercept
    n_found <- 0
    for (sid in required_somaids) {
      if (sid %in% names(row)) {
        val <- row[[sid]]
        if (!is.na(val) && is.numeric(val)) {
          if (log_transform && val > 0) {
            val <- log2(val)
          } else if (log_transform && val <= 0) {
            val <- NA_real_
          }
          if (!is.na(val)) {
            pred <- pred + coef_vec[sid] * val
            n_found <- n_found + 1
          }
        }
      }
    }
    prot_age[i] <- pred
    if (i == 1 && n_found < length(required_somaids)) {
      warning("Only ", n_found, " of ", length(required_somaids),
              " proteins found in the first sample. ",
              "Proteomic age estimates may be unreliable.")
    }
  }

  age_accel <- prot_age - chron_age

  result <- data.frame(
    id                = ids,
    chronological_age = chron_age,
    proteomic_age     = prot_age,
    age_acceleration  = age_accel,
    n_proteins_matched = length(required_somaids) - length(missing),
    n_proteins_missing = length(missing),
    stringsAsFactors   = FALSE
  )

  return(result)
}
