.olink_gene_lookup <- function(columns, genes) {
  lookup <- stats::setNames(columns, toupper(columns))
  syntactic <- stats::setNames(make.names(genes, unique = FALSE), toupper(genes))
  for (gene in names(syntactic)) {
    key <- toupper(syntactic[[gene]])
    if (key %in% names(lookup) && !(gene %in% names(lookup))) {
      lookup[[gene]] <- lookup[[key]]
    }
  }
  lookup
}

.olink_match_lookup <- function(columns, proteins, match_by, protein_map = NULL) {
  if (match_by == "gene") return(.olink_gene_lookup(columns, proteins))
  if (is.null(protein_map)) {
    load_olink_protein_map()
    protein_map <- .olink_map_cache$map
  }

  map <- .olink_standardize_map(protein_map)
  map <- map[toupper(map$gene) %in% toupper(proteins), ]
  input_lookup <- stats::setNames(columns, toupper(columns))
  matched_map <- map[toupper(map$uniprot) %in% names(input_lookup), ]
  stats::setNames(input_lookup[toupper(matched_map$uniprot)], toupper(matched_map$gene))
}

.olink_standardize_map <- function(protein_map) {
  if (!is.data.frame(protein_map)) stop("'protein_map' must be a data.frame")
  names_lower <- tolower(names(protein_map))
  gene_col <- match(TRUE, names_lower %in% c("gene", "assay", "olink.assay"))
  uniprot_col <- match(TRUE, names_lower %in% c("uniprot", "uniprot_id", "olink.uniprot_id"))
  if (is.na(gene_col) || is.na(uniprot_col)) {
    stop("protein_map must contain gene/assay and uniprot columns")
  }

  out <- data.frame(
    gene = as.character(protein_map[[gene_col]]),
    uniprot = as.character(protein_map[[uniprot_col]]),
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$gene) & out$gene != "" &
               !is.na(out$uniprot) & out$uniprot != "", ]
  out[!duplicated(paste(toupper(out$gene), toupper(out$uniprot))), ]
}

#' List built-in Olink gene-UniProt mappings
#'
#' @return data.frame with Olink assay/gene names, UniProt accessions, assay IDs,
#'   panels, and source.
#' @export
olink_protein_map <- function() {
  load_olink_protein_map()
  .olink_map_cache$map
}

#' Build an Olink gene-UniProt map from long-format Olink data
#'
#' @param data Long-format Olink data frame.
#' @param gene_col Column containing Olink assay or gene names.
#' @param uniprot_col Column containing UniProt accessions.
#' @return data.frame with gene and uniprot columns. This is optional; Olink
#'   clock functions use a built-in Olink map when `match_by = "uniprot"`.
#' @export
olink_protein_map_from_long <- function(data,
                                        gene_col = "olink.assay",
                                        uniprot_col = "olink.uniprot_id") {
  if (!is.data.frame(data)) stop("'data' must be a data.frame")
  if (!gene_col %in% names(data)) stop("gene_col not found")
  if (!uniprot_col %in% names(data)) stop("uniprot_col not found")

  .olink_standardize_map(data.frame(
    gene = data[[gene_col]],
    uniprot = data[[uniprot_col]],
    stringsAsFactors = FALSE
  ))
}
