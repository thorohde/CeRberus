#' Access fields in an internal screen design
#'
#' @description
#' These methods provide read-only list-like access to fields stored in an
#' internal [`ScreenDesign-class`] object.
#'
#' @param x A `ScreenDesign` object.
#' @param name Character scalar naming a design field.
#' @param i Character scalar naming a design field.
#'
#' @return The value stored in the selected field.
#'
#' @keywords internal
#' @name extract-ScreenDesign
#' @aliases $,ScreenDesign-method
setMethod("$", "ScreenDesign", function(x, name) {
  methods::slot(x, name)
})

#' @rdname extract-ScreenDesign
#' @aliases [[,ScreenDesign,character-method
setMethod("[[", signature(x = "ScreenDesign", i = "character"), function(x, i) {
  if (length(i) != 1L) {
    stop("ScreenDesign fields must be selected one at a time.", call. = FALSE)
  }

  methods::slot(x, i)
})
