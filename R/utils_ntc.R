#' Label result pairs involving non-targeting controls
#'
#' @keywords internal
label_ntc_pairs <- function(results, non_targeting_controls = NULL) {
  stopifnot(
    "results must be a data frame or data table." = is.data.frame(results),
    "results must contain query_gene and library_gene columns." = all(
      c("query_gene", "library_gene") %in% colnames(results)
    ),
    "non_targeting_controls must be NULL or a character vector." = is.null(
      non_targeting_controls
    ) || is.character(non_targeting_controls)
  )

  output <- data.table::as.data.table(data.table::copy(results))
  output[, has_NTC :=
    query_gene %in% non_targeting_controls |
      library_gene %in% non_targeting_controls]

  output
}

#' Store non-targeting-control pair annotations
#'
#' @keywords internal
store_ntc_pair_annotations <- function(gi_obj) {
  results <- gi_df(gi_obj)
  annotations <- label_ntc_pairs(
    results,
    gi_obj@metadata$non_targeting_controls
  )[, .(gene_pair, has_NTC)]

  gi_obj@metadata$ntc_pair_annotations <- annotations
  gi_obj
}

#' Append stored non-targeting-control annotations
#'
#' @keywords internal
append_ntc_pair_annotations <- function(results, gi_obj) {
  annotations <- gi_obj@metadata$ntc_pair_annotations
  if (is.null(annotations)) {
    return(results)
  }

  output <- data.table::as.data.table(data.table::copy(results))
  output[annotations, has_NTC := i.has_NTC, on = "gene_pair"]
  output
}
