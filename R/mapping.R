# Protein name normalization for cross-format SOMAscan data
#
# Supports mapping between:
#   - Gene symbols (e.g., "GDF15")
#   - SeqId dot format (e.g., "4374.45.2")
#   - SeqId SL format (e.g., "SL003869")
#   - SeqId full format (e.g., "GDF15.4374.45.2")
#
# Uses built-in mapping table; optionally integrates with SomaScan.db
# for broader coverage.

.protein_mapping <- NULL

load_protein_mapping <- function() {
  if (!is.null(.protein_mapping)) return(invisible())
  path <- system.file("extdata", "protein_mapping.csv",
                      package = "proteomicAge", mustWork = FALSE)
  if (!file.exists(path)) {
    warning("protein_mapping.csv not found. Protein name normalization disabled.")
    .protein_mapping <<- data.frame()
    return(invisible())
  }
  .protein_mapping <<- utils::read.csv(path, stringsAsFactors = FALSE)
  invisible()
}

# Detect the format of protein column names
detect_column_format <- function(colnames_vec) {
  load_protein_mapping()
  mp <- .protein_mapping

  if (nrow(mp) == 0) return("unknown")

  n_sl    <- sum(colnames_vec %in% mp$seqid_sl, na.rm = TRUE)
  n_dot   <- sum(colnames_vec %in% mp$seqid_dot, na.rm = TRUE)
  n_full  <- sum(colnames_vec %in% mp$seqid_full, na.rm = TRUE)
  n_gene  <- sum(colnames_vec %in% mp$gene_tanaka, na.rm = TRUE) +
             sum(colnames_vec %in% mp$gene_lehallier, na.rm = TRUE)

  best <- which.max(c(n_sl, n_dot, n_full, n_gene))
  c("seqid_sl", "seqid_dot", "seqid_full", "gene")[best]
}

# Normalize protein column names to a target format
#
# @param data Data frame with protein columns
# @param id_col, age_col Non-protein columns to preserve
# @param target_format One of "seqid_sl", "seqid_full"
# @param input_format One of "auto", "gene", "seqid_sl", "seqid_dot", "seqid_full"
# @return Data frame with renamed protein columns
normalize_protein_names <- function(data,
                                     id_col = "SampleID",
                                     age_col = "Age",
                                     target_format = c("seqid_sl", "seqid_full"),
                                     input_format = c("auto", "gene", "seqid_sl",
                                                      "seqid_dot", "seqid_full")) {

  target_format <- match.arg(target_format)
  input_format  <- match.arg(input_format)

  load_protein_mapping()
  mp <- .protein_mapping

  if (nrow(mp) == 0) {
    warning("No protein mapping available. Returning data unchanged.")
    return(data)
  }

  # Identify protein columns (everything except id/age and known non-protein cols)
  skip_cols <- c(id_col, age_col)
  protein_cols <- setdiff(names(data), skip_cols)

  if (length(protein_cols) == 0) {
    stop("No protein columns found in data (excluding id/age columns)")
  }

  # Auto-detect input format
  if (input_format == "auto") {
    input_format <- detect_column_format(protein_cols)
    message("Auto-detected input format: ", input_format)
  }

  # Build lookup from input format to target format
  lookup <- switch(
    input_format,
    gene = {
      # Gene → target format. Use clock-specific gene column to avoid clashes.
      if (target_format == "seqid_sl") {
        stats::setNames(mp[["seqid_sl"]], mp[["gene_tanaka"]])
      } else {
        stats::setNames(mp[["seqid_full"]], mp[["gene_lehallier"]])
      }
    },
    seqid_sl = {
      if (target_format == "seqid_sl") return(data)  # already correct
      # SL → gene → seqid_full
      gene_lk <- stats::setNames(mp$gene_clean, mp$seqid_sl)
      full_lk <- stats::setNames(mp$seqid_full, mp$gene_clean)
      # Chain: SL → gene → full
      lk <- character()
      for (sl in names(gene_lk)) {
        gene <- gene_lk[sl]
        if (!is.na(gene) && gene %in% names(full_lk)) {
          lk[sl] <- full_lk[gene]
        }
      }
      lk
    },
    seqid_dot = {
      target_col <- if (target_format == "seqid_sl") "seqid_sl" else "seqid_full"
      stats::setNames(mp[[target_col]], mp$seqid_dot)
    },
    seqid_full = {
      if (target_format == "seqid_full") return(data)
      # Full → gene → seqid_sl
      gene_lk <- stats::setNames(mp$gene_clean, mp$seqid_full)
      sl_lk   <- stats::setNames(mp$seqid_sl, mp$gene_clean)
      lk <- character()
      for (full in names(gene_lk)) {
        gene <- gene_lk[full]
        if (!is.na(gene) && gene %in% names(sl_lk)) {
          lk[full] <- sl_lk[gene]
        }
      }
      lk
    }
  )

  # Remove NAs from lookup
  lookup <- lookup[!is.na(lookup) & lookup != ""]

  # Rename columns
  new_names <- names(data)
  renamed <- 0
  for (i in seq_along(new_names)) {
    col <- new_names[i]
    if (col %in% skip_cols) next
    if (col %in% names(lookup)) {
      new_names[i] <- lookup[col]
      renamed <- renamed + 1
    }
  }

  if (renamed > 0) {
    message("Renamed ", renamed, " protein columns from ", input_format,
            " to ", target_format)
  }

  names(data) <- new_names
  data
}
