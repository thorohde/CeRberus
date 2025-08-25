#' Alphabetically sorts gene pairs. 
#' 
#' @description
#' After pairwise sorting, `sort_gene_pairs` concatenates two vectors of gene names into one. 
#' 
#' @param g1 Character vector of gene names
#' @param g2 Character vector of gene names
#' @param sep Separator used to concatenate the gene pairs
#' @param pairs Vector of gene pair names, separated by the pair_sep separator
#' @param pair_sep Separator used when pairs is provided
#' @param invert inversts the sorting of gene pairs. If false, the output will be "GeneA;GeneB". If true: "GeneB;GeneA"
#' @examples
#' genes1 <- c("RB1", "NOTCH1", "TTN", "MSH2", "FANCD2")
#' genes2 <- c("RNF43", "NEK1", "HLTF", "REV3L", "PAPD7")
#' 
#' # Basic usage
#' sort_gene_pairs(genes1, genes2)
#' 
#' # Inverted sorting
#' sort_gene_pairs(genes1, genes2, invert = TRUE)

#' 
#' @export sort_gene_pairs

sort_gene_pairs <- function(g1, g2, sep = ";", pairs, pair_sep = ";", invert = F) {
  
  stopifnot("invert needs to be logical." = is.logical(invert), 
            "No input pairs given!" = !(missing(pairs) & missing(g1) & missing(g2)))
  
  if (!(missing(g1) | missing(g2))) {
    "g1 and g2 have to be of equal length!" = (!missing(g1) & !missing(g2) & length(g1) == length(g2))
  }
  
  stopifnot(
    "Either two gene vectors or a gene pair vector required" = !(missing(pairs) & (missing(g1) | missing(g2)))
  )
  
  if (!missing(pairs) & missing(g1) & missing(g2)) {
    g1 <- str_split_i(pairs, pattern = pair_sep, 1)
    g2 <- str_split_i(pairs, pattern = pair_sep, 2)
  }
  
  return(switch(as.character(invert),
                "TRUE" = paste(pmax(g1, g2), pmin(g1, g2), sep = sep), 
                "FALSE" = paste(pmin(g1, g2), pmax(g1, g2), sep = sep)))
  
}

