#' @export sort_gene_pairs

sort_gene_pairs <- function(genes1, genes2, sep = ";") {
  base::paste(base::pmin(genes1, genes2), base::pmax(genes1, genes2), sep = sep)
}
