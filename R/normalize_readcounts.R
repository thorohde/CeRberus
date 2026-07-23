#' Normalize read-count values
#'
#' @description
#' Converts raw read counts to log2-transformed, library-size-normalized values
#' and shifts the minimum finite value to zero. Missing values are preserved.
#'
#' @param readcounts Numeric vector of read counts.
#' @param cf1 Numeric scaling factor applied after library-size normalization.
#' @param cf2 Numeric pseudo-count added before log2 transformation.
#'
#' @return Numeric vector with the same length as `readcounts`.
#'
#' @export

#####

normalize_readcounts <- function(readcounts, cf1 = 100, cf2 = 1) {
  stopifnot(
    "readcounts must be numeric." = is.numeric(readcounts),
    "cf1 must be a single non-negative number." = length(cf1) == 1 &&
      is.finite(cf1) &&
      cf1 >= 0,
    "cf2 must be a single non-negative number." = length(cf2) == 1 &&
      is.finite(cf2) &&
      cf2 >= 0
  )

  total_counts <- sum(readcounts, na.rm = TRUE)
  if (!is.finite(total_counts) || total_counts <= 0) {
    return(rep(NA_real_, length(readcounts)))
  }

  # replaced 1e6 with 100 x length(readcounts), cf1 = 1e6, cf2 = 0.5
  x <- log2(
    (readcounts / total_counts) *
      cf1 *
      length(readcounts) +
      cf2
  ) #NA will stay NA
  if (!all(is.na(x))) {
    #not run if replicate is missing (= all NA)
    x <- x - min(x, na.rm = TRUE)
  } #smallest value is 0 regardless of cf2
  x
}

#####
