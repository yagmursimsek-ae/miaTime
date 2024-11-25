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

#' @rdname addMetrics
#' @export
setGeneric("addMetrics", signature = "x", function(x, ...)
    standardGeneric("addMetrics"))
