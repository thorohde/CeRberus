#' Remove selected principal components from a data matrix
#'
#' @description
#' Performs principal component analysis on a numeric matrix-like object, sets
#' selected principal component scores to zero, and reconstructs the data in the
#' original variable space. This can be used to remove dominant axes of variation
#' such as batch effects or other unwanted structure captured by specific PCs.
#'
#' If `to_remove` is `NA`, no PCs are removed and `.x` is returned unchanged.
#'
#' @param .x Numeric matrix-like object with observations in rows and variables
#'   in columns. Passed to [stats::prcomp()].
#' @param to_remove Integer vector of principal component indices to remove. Use
#'   `NA` to skip PC removal.
#' @param .center Logical scalar passed to the `center` argument of
#'   [stats::prcomp()]. If `TRUE`, variables are centered before PCA and the
#'   center is added back after reconstruction.
#' @param .scale Logical scalar passed to the `scale.` argument of
#'   [stats::prcomp()]. If `TRUE`, variables are scaled before PCA and the scale
#'   is applied back after reconstruction.
#'
#' @return A numeric matrix with the same dimensions as `.x`, reconstructed after
#'   setting the selected principal component scores to zero. If `to_remove` is
#'   `NA`, returns `.x` unchanged.
#'
#' @examples
#' x <- matrix(rnorm(50), nrow = 10)
#' remove_PCs(x, to_remove = 1)
#' remove_PCs(x, to_remove = NA)
#'
#' @export

remove_PCs <- \(.x, to_remove = NA, .center = TRUE, .scale = TRUE) {
  pca_result <- prcomp(.x, center = .center, scale. = .scale)
  pcs <- pca_result$x
  center <- pca_result$center
  scale <- pca_result$scale

  if (all(is.na(to_remove))) {
    return(.x)
  }

  to_remove <- unique(to_remove) # Validate to_remove indices
  if (any(to_remove > ncol(pcs))) {
    stop(
      "Some PCs to remove exceed the number of available principal components."
    )
  }

  pcs[, to_remove] <- 0
  .x <- pcs %*% t(pca_result$rotation)
  if (.scale) {
    .x <- sweep(.x, 2, scale, FUN = "*")
  }
  if (.center) {
    .x <- sweep(.x, 2, center, FUN = "+")
  }
  return(.x)
}
