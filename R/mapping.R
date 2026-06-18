# Protein name normalization using SomaScan.db (primary) + built-in fallback
#
# Architecture:
#   1. SomaScan.db is the PRIMARY mapper for gene, uniprot, PROBEID formats
#   2. Built-in protein_mapping.csv is a FALLBACK for:
#      - SL format (v3 SOMAscan, not in SomaScan.db)
#      - When SomaScan.db is not installed
#
# Supported input formats:
#   "gene"       - Gene symbols (e.g., "GDF15")
#   "uniprot"    - UniProt accessions (e.g., "Q99988")
#   "seqid_dot"  - Dot-format SeqIds (e.g., "seq.4374.45.2" or "4374.45.2")
#   "seqid_sl"   - SL-format SeqIds (e.g., "SL003869") - only via built-in
#   "seqid_full" - Gene.SeqId format (e.g., "GDF15.4374.45.2")

# ---------- Internal state (environment to avoid locked binding on reload) ----------

.mapping_state <- new.env(parent = emptyenv())

check_somascan_db <- function() {
  if (!is.null(.mapping_state$has_somascan_db)) return(.mapping_state$has_somascan_db)
  has <- requireNamespace("SomaScan.db", quietly = TRUE)
  .mapping_state$has_somascan_db <- has
  if (has) {
    message("Using SomaScan.db for protein name mapping")
  } else {
    message("SomaScan.db not available. Using built-in mapping (limited coverage).")
    message("Install with: BiocManager::install('SomaScan.db')")
  }
  has
}

load_protein_mapping <- function() {
  if (!is.null(.mapping_state$protein_mapping)) return(invisible())
  path <- system.file("extdata", "protein_mapping.csv",
                      package = "proteomicAge", mustWork = FALSE)
  if (!file.exists(path)) {
    .mapping_state$protein_mapping <- data.frame()
    return(invisible())
  }
  .mapping_state$protein_mapping <- utils::read.csv(path, stringsAsFactors = FALSE)
  invisible()
}

# ---------- Format converters ----------

dot_to_probeid <- function(dot_seqid) {
  # Strip "seq." prefix if present (e.g., "seq.4374.45.2" -> "4374.45.2")
  x <- as.character(dot_seqid)
  x <- sub("^seq[.]", "", x, ignore.case = TRUE)
  parts <- strsplit(x, "[.]")[[1]]
  if (length(parts) >= 3) {
    nums <- utils::tail(parts, 3)
    paste0(nums[1], "-", nums[2], "_", nums[3])
  } else {
    NA_character_
  }
}

probeid_to_dot <- function(probeid) {
  gsub("-", ".", gsub("_", ".", probeid, fixed = TRUE), fixed = TRUE)
}

# Extract the dot-format SeqId from a seqid_full string like "GDF15.4374.45.2"
extract_dot_from_full <- function(full) {
  sapply(strsplit(full, "[.]"), function(x) {
    if (length(x) >= 3) paste(utils::tail(x, 3), collapse = ".") else NA_character_
  })
}

# ---------- SomaScan.db based mapping ----------

map_via_somascan <- function(keys, from_type, to_type, multi_val = "first") {
  if (!check_somascan_db()) return(rep(NA_character_, length(keys)))

  keytype_map <- c(
    gene    = "ALIAS",
    uniprot = "UNIPROT",
    probeid = "PROBEID",
    symbol  = "SYMBOL"
  )
  column_map <- c(
    gene    = "SYMBOL",
    uniprot = "SYMBOL",
    probeid = "SYMBOL",
    symbol  = "PROBEID",
    probeid_id = "PROBEID"
  )

  from_keytype <- keytype_map[from_type]
  to_column    <- column_map[to_type]

  tryCatch({
    result <- AnnotationDbi::mapIds(
      SomaScan.db::SomaScan.db,
      keys = keys,
      column = to_column,
      keytype = from_keytype,
      multiVals = multi_val
    )
    as.character(result)
  }, error = function(e) {
    rep(NA_character_, length(keys))
  })
}

# ---------- Main normalization function ----------

#' Normalize protein column names to a target format
#'
#' Uses SomaScan.db (primary) for comprehensive mapping across all SOMAscan
#' proteins, falling back to a built-in mapping table for SL-format SeqIds
#' (which predate SomaScan.db's coverage).
#'
#' @param data Data frame with protein columns
#' @param id_col, age_col Non-protein columns to preserve
#' @param target_format Target format: "seqid_sl" (Tanaka) or "seqid_full" (Lehallier)
#' @param input_format Input format: "auto", "gene", "uniprot", "seqid_dot", "seqid_sl", "seqid_full"
#' @return Data frame with renamed protein columns
#' @export
normalize_protein_names <- function(data,
                                     id_col = "SampleID",
                                     age_col = "Age",
                                     target_format = c("seqid_sl", "seqid_full"),
                                     input_format = c("auto", "gene", "uniprot",
                                                      "seqid_dot", "seqid_sl",
                                                      "seqid_full")) {

  target_format <- match.arg(target_format)
  input_format  <- match.arg(input_format)

  load_protein_mapping()
  mp <- .mapping_state$protein_mapping

  skip_cols <- c(id_col, age_col)
  protein_cols <- setdiff(names(data), skip_cols)

  if (length(protein_cols) == 0) {
    stop("No protein columns found in data (excluding id/age columns)")
  }

  has_db <- check_somascan_db()

  # ---------- Step 1: Resolve each protein column to a (gene, probeid) pair ----------
  resolved <- data.frame(
    original = protein_cols,
    gene     = NA_character_,
    probeid  = NA_character_,
    stringsAsFactors = FALSE
  )

  if (input_format == "auto") {
    input_format <- detect_column_format(protein_cols)
    message("Auto-detected input format: ", input_format)
  }

  for (i in seq_along(protein_cols)) {
    col <- protein_cols[i]

    if (input_format == "gene") {
      resolved$gene[i] <- col
      if (has_db) {
        resolved$probeid[i] <- map_via_somascan(col, "gene", "probeid_id")
      }

    } else if (input_format == "uniprot") {
      if (has_db) {
        resolved$gene[i]    <- map_via_somascan(col, "uniprot", "gene")
        resolved$probeid[i] <- map_via_somascan(col, "uniprot", "probeid_id")
      }

    } else if (input_format == "seqid_dot") {
      resolved$probeid[i] <- dot_to_probeid(col)
      if (has_db) {
        resolved$gene[i] <- map_via_somascan(resolved$probeid[i], "probeid", "gene")
      }

    } else if (input_format == "seqid_full") {
      dot <- extract_dot_from_full(col)
      resolved$probeid[i] <- dot_to_probeid(dot)
      if (has_db) {
        resolved$gene[i] <- map_via_somascan(resolved$probeid[i], "probeid", "gene")
      } else {
        resolved$gene[i] <- strsplit(col, "[.]")[[1]][1]
      }

    } else if (input_format == "seqid_sl") {
      resolved$gene[i] <- col
      # SL format is NOT in SomaScan.db — use built-in mapping
      if (nrow(mp) > 0) {
        match_row <- mp[mp$seqid_sl == col & !is.na(mp$seqid_sl), ]
        if (nrow(match_row) > 0) {
          resolved$gene[i] <- match_row$gene_simple[1]
        }
        # Try to get PROBEID via gene → SomaScan.db
        if (has_db && !is.na(resolved$gene[i])) {
          resolved$probeid[i] <- map_via_somascan(resolved$gene[i], "gene", "probeid_id")
        }
      }
    }
  }

  # ---------- Step 2: Build target names ----------
  target_names <- character(length(protein_cols))

  for (i in seq_along(protein_cols)) {
    gene    <- resolved$gene[i]
    probeid <- resolved$probeid[i]

    if (target_format == "seqid_sl") {
      # Tanaka: need SL-format (only via built-in mapping)
      if (nrow(mp) > 0) {
        if (!is.na(gene)) {
          match_row <- mp[mp$gene_simple == gene & !is.na(mp$seqid_sl), ]
          if (nrow(match_row) == 0) {
            match_row <- mp[mp$gene_tanaka == gene & !is.na(mp$seqid_sl), ]
          }
          if (nrow(match_row) > 0) {
            target_names[i] <- match_row$seqid_sl[1]
            next
          }
        }
        # Also try direct SL lookup
        match_row <- mp[mp$seqid_sl == protein_cols[i] & !is.na(mp$seqid_sl), ]
        if (nrow(match_row) > 0) {
          target_names[i] <- match_row$seqid_sl[1]
          next
        }
      }
      # Fallback: keep original name
      target_names[i] <- protein_cols[i]

    } else {
      # Lehallier: need seqid_full format (Gene.ProbeIdAsDot)
      # First try built-in mapping for exact match
      target_set <- FALSE
      if (nrow(mp) > 0 && !is.na(gene)) {
        match_row <- mp[mp$gene_simple == gene & !is.na(mp$seqid_full), ]
        if (nrow(match_row) == 0) {
          match_row <- mp[mp$gene_lehallier == gene & !is.na(mp$seqid_full), ]
        }
        if (nrow(match_row) > 0) {
          target_names[i] <- match_row$seqid_full[1]
          target_set <- TRUE
        }
      }
      if (target_set) next

      probeid_dot <- NA_character_
      if (!is.na(probeid)) {
        probeid_dot <- probeid_to_dot(probeid)
      } else if (input_format %in% c("seqid_dot", "seqid_full")) {
        if (input_format == "seqid_full") {
          probeid_dot <- extract_dot_from_full(protein_cols[i])
        } else {
          # Strip seq. prefix: "seq.4374.45.2" -> "4374.45.2"
          probeid_dot <- sub("^seq[.]", "", protein_cols[i], ignore.case = TRUE)
        }
      }

      if (!is.na(gene) && !is.na(probeid_dot) && gene != "" && probeid_dot != "") {
        target_names[i] <- paste0(gene, ".", probeid_dot)
      } else if (!is.na(gene) && gene != "") {
        target_names[i] <- gene
      } else {
        target_names[i] <- protein_cols[i]
      }
    }
  }

  # ---------- Step 3: Rename columns ----------
  rename_map <- stats::setNames(target_names, protein_cols)
  rename_map <- rename_map[rename_map != protein_cols]
  rename_map <- rename_map[!is.na(rename_map) & rename_map != ""]

  if (length(rename_map) > 0) {
    for (old in names(rename_map)) {
      if (old %in% names(data)) {
        names(data)[names(data) == old] <- rename_map[old]
      }
    }
    message("Renamed ", length(rename_map), " protein columns from ",
            input_format, " to ", target_format, " format")
  }

  data
}

# ---------- Format detection ----------

#' Detect the format of protein column names
#' @export
detect_column_format <- function(colnames_vec) {
  load_protein_mapping()
  mp <- .mapping_state$protein_mapping

  n_sl    <- sum(colnames_vec %in% mp$seqid_sl, na.rm = TRUE)
  n_dot   <- sum(grepl("^(seq[.])?[0-9]+[.][0-9]+[.][0-9]+$", colnames_vec, ignore.case = TRUE))
  n_full  <- sum(grepl("^[A-Za-z0-9]+[.][0-9]+[.][0-9]+[.][0-9]+$", colnames_vec))
  n_gene  <- sum(colnames_vec %in% mp$gene_tanaka, na.rm = TRUE) +
             sum(colnames_vec %in% mp$gene_lehallier, na.rm = TRUE)
  n_uniprot <- sum(grepl("^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$", colnames_vec))

  scores <- c(seqid_sl = n_sl, seqid_dot = n_dot, seqid_full = n_full,
              gene = n_gene, uniprot = n_uniprot)
  best <- which.max(scores)
  if (scores[best] == 0) return("gene")
  names(scores)[best]
}
