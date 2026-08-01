make_screen_design <- function(
  query_genes = character(),
  library_genes = character(),
  all_pairs = character(),
  observations_per_query = NULL,
  contrasts = NULL,
  unique_pairs = NULL,
  ...
) {
  query_genes <- unique(as.character(query_genes))
  library_genes <- unique(as.character(library_genes))
  all_genes <- union(query_genes, library_genes)

  if (is.null(observations_per_query)) {
    observations_per_query <- stats::setNames(
      integer(length(query_genes)),
      query_genes
    )
  } else {
    observations_per_query <- as.integer(observations_per_query)
    names(observations_per_query) <- query_genes
  }

  if (is.null(unique_pairs)) {
    if (length(all_pairs) == 0L) {
      unique_pairs <- character()
    } else {
      pair_parts <- strsplit(all_pairs, ";", fixed = TRUE)
      unique_pairs <- unique(sort_gene_pairs(
        purrr::map_chr(pair_parts, 1L),
        purrr::map_chr(pair_parts, 2L)
      ))
    }
  }

  methods::new(
    "ScreenDesign",
    contrasts = contrasts,
    query_genes = query_genes,
    library_genes = library_genes,
    all_genes = all_genes,
    query_genes_not_in_lib = setdiff(query_genes, library_genes),
    library_genes_not_in_query = setdiff(library_genes, query_genes),
    n_query_genes = length(query_genes),
    n_lib_genes = length(library_genes),
    n_all_genes = length(all_genes),
    observations_per_query = observations_per_query,
    all_pairs = as.character(all_pairs),
    unique_pairs = as.character(unique_pairs),
    ...
  )
}
