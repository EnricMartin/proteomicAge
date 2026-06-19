# Protein mapping: UniProt-based normalization using prolist.csv
#
# Architecture:
#   prolist.csv (7,596 proteins) is the SOLE mapping source.
#   Input protein columns MUST be UniProt accessions.
#   Mapping flow: UniProt → prolist.csv → target SeqId format

.prolist <- NULL
.prolist_lookup <- NULL

load_prolist <- function() {
  if (!is.null(.prolist)) return(invisible())
  path <- system.file("extdata", "prolist.csv",
                      package = "proteomicAge", mustWork = TRUE)
  .prolist <<- utils::read.csv(path, stringsAsFactors = FALSE)
  .prolist_lookup <<- split(.prolist, .prolist$UniProt)
  invisible()
}

# Map a UniProt accession to the clock's target SeqId
map_uniprot_to_target <- function(uniprot, target_format) {
  load_prolist()
  col <- if (target_format == "seqid_sl") "seqid_sl" else "seqid_dot"
  pl <- .prolist
  idx <- which(pl$UniProt == uniprot)
  if (length(idx) == 0) return(NA_character_)
  sid <- pl[[col]][idx[1]]
  if (target_format == "seqid_full") {
    gene <- pl$EntrezGeneSymbol[idx[1]]
    if (!is.na(gene) && gene != "" && !is.na(sid) && sid != "") {
      return(paste0(gene, ".", sid))
    }
  }
  sid
}

#' Normalize UniProt protein columns to a clock-specific SeqId format
#'
#' Columns named by UniProt accessions (e.g., "P36222") are mapped to
#' the target clock's SeqId format using prolist.csv.
#'
#' @param data Data frame with UniProt-named protein columns
#' @param id_col, age_col Non-protein columns to preserve
#' @param target_format Target format: "seqid_sl" (Tanaka) or "seqid_full" (Lehallier)
#' @return Data frame with renamed protein columns
#' @export
normalize_protein_names <- function(data,
                                     id_col = "SampleID",
                                     age_col = "Age",
                                     target_format = c("seqid_sl", "seqid_full")) {

  target_format <- match.arg(target_format)
  load_prolist()

  skip_cols <- c(id_col, age_col, "Sex")
  all_cols <- names(data)
  protein_cols <- setdiff(all_cols, skip_cols)

  if (length(protein_cols) == 0) {
    stop("No protein columns found (excluding id, age, sex columns).")
  }

  col_to_use <- if (target_format == "seqid_sl") "seqid_sl" else "seqid_dot"
  gene_col <- "EntrezGeneSymbol"

  rename_map <- character()
  for (col in protein_cols) {
    idx <- which(.prolist$UniProt == col)
    if (length(idx) > 0) {
      sid <- .prolist[[col_to_use]][idx[1]]
      if (target_format == "seqid_full" && !is.na(sid) && sid != "") {
        gene <- .prolist[[gene_col]][idx[1]]
        if (!is.na(gene) && gene != "") {
          rename_map[col] <- paste0(gene, ".", sid)
        }
      } else if (!is.na(sid) && sid != "") {
        rename_map[col] <- sid
      }
    }
  }

  if (length(rename_map) > 0) {
    for (old in names(rename_map)) {
      names(data)[names(data) == old] <- rename_map[old]
    }
    message("Mapped ", length(rename_map), " UniProt columns to ", target_format)
  }

  data
}
