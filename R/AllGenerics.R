# All generic methods are listed here

#' @rdname addBaselineDivergence
#' @export
setGeneric("getBaselineDivergence", signature = "x", function(x, ...)
    standardGeneric("getBaselineDivergence"))

#' @rdname addBaselineDivergence
#' @export
setGeneric("addBaselineDivergence", signature = "x", function(x, ...)
    standardGeneric("addBaselineDivergence"))

#' @rdname addStepwiseDivergence
#' @export
#'
setGeneric("getStepwiseDivergence", signature = c("x"), function(x, ...)
    standardGeneric("getStepwiseDivergence"))

#' @rdname addStepwiseDivergence
#' @export
setGeneric("addStepwiseDivergence", signature = "x", function(x, ...)
    standardGeneric("addStepwiseDivergence"))

#' @rdname getBimodality
#' @export
setGeneric("getBimodality", signature = "x", function(x, ...)
    standardGeneric("getBimodality"))

#' @rdname getBimodality
#' @export
setGeneric("addBimodality", signature = "x", function(x, ...)
    standardGeneric("addBimodality"))
