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
