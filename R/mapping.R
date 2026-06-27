# Protein format detection and conversion using prolist.csv
#
# Supported naming conventions:
#   "gene"      - Entrez gene symbols (e.g., "GDF15", "CHI3L1")
#   "seqid_sl"  - SL-format SeqIds (e.g., "SL003869")
#   "seqid_dot" - Dot-format SeqIds (e.g., "seq.4374.45.2" or "4374.45.2")
#   "uniprot"   - UniProt accessions (e.g., "P36222", "Q99988")

.plist <- new.env(parent = emptyenv())

load_plist <- function() {
  if (!is.null(.plist$data)) return(invisible())
  path <- system.file("extdata", "prolist.csv",
                      package = "proteomicAge", mustWork = TRUE)
  .plist$data <- utils::read.csv(path, stringsAsFactors = FALSE)
  invisible()
}

#' Detect the naming convention of protein columns
#'
#' @param colnames_vec Character vector of column names
#' @return One of "gene", "seqid_sl", "seqid_dot", "uniprot"
#' @export
detect_format <- function(colnames_vec) {
  load_plist()
  pl <- .plist$data
  n_sl   <- sum(colnames_vec %in% pl$seqid_sl, na.rm = TRUE)
  n_dot  <- sum(grepl("^(seq[.])?[0-9]+[.][0-9]+([.][0-9]+)?$", colnames_vec, ignore.case = TRUE))
  n_full <- sum(grepl("^[A-Za-z0-9_]+[.][0-9]+[.][0-9]+[.][0-9]+$", colnames_vec))
  n_gene <- sum(colnames_vec %in% pl$EntrezGeneSymbol, na.rm = TRUE)
  n_up   <- sum(grepl("^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$", colnames_vec))
  scores <- c(gene = n_gene, seqid_sl = n_sl, seqid_dot = n_dot + n_full, uniprot = n_up)
  best <- which.max(scores)
  if (scores[best] == 0) "gene" else names(scores)[best]
}

#' Convert protein column names between formats
#'
#' Uses prolist.csv to map between gene symbols, SeqIds, and UniProt.
#'
#' @param data Data frame with protein columns
#' @param target_format Target naming convention
#' @param id_col, age_col Non-protein columns to preserve
#' @return Data frame with renamed protein columns
#' @export
convert_format <- function(data,
                            target_format = c("gene", "seqid_sl", "seqid_dot", "uniprot"),
                            id_col = "SampleID",
                            age_col = "Age") {

  target_format <- match.arg(target_format)
  load_plist()
  pl <- .plist$data

  skip <- c(id_col, age_col, "Sex")
  prot_cols <- setdiff(names(data), skip)

  target_col <- switch(target_format,
    gene      = "EntrezGeneSymbol",
    seqid_sl  = "seqid_sl",
    seqid_dot = "seqid_dot",
    uniprot   = "UniProt"
  )

  rename_map <- character()
  for (col in prot_cols) {
    new_name <- NA_character_
    for (src in c("seqid_sl", "seqid_dot", "EntrezGeneSymbol", "UniProt")) {
      idx <- which(pl[[src]] == col)
      if (length(idx) > 0) { new_name <- pl[[target_col]][idx[1]]; break }
    }
    if (is.na(new_name)) {
      stripped <- sub("^seq[.]", "", col, ignore.case = TRUE)
      idx <- which(pl$seqid_dot2 == stripped)
      if (length(idx) > 0) new_name <- pl[[target_col]][idx[1]]
    }
    if (!is.na(new_name) && new_name != "" && new_name != col) {
      rename_map[col] <- new_name
    }
  }

  if (length(rename_map) > 0) {
    for (old in names(rename_map)) {
      names(data)[names(data) == old] <- rename_map[old]
    }
    message("Converted ", length(rename_map), " columns to ", target_format)
  }
  data
}
