# All generic methods are listed here

#' @rdname getBaselineDivergence
#' @export
setGeneric("getBaselineDivergence", signature = "x", function(x, ...)
    standardGeneric("getBaselineDivergence"))

#' @rdname getBaselineDivergence
#' @export
setGeneric("addBaselineDivergence", signature = "x", function(x, ...)
    standardGeneric("addBaselineDivergence"))

#' @rdname getStepwiseDivergence
#' @export
#'
setGeneric("getStepwiseDivergence", signature = c("x"), function(x, ...)
    standardGeneric("getStepwiseDivergence"))

#' @rdname getStepwiseDivergence
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

#' @rdname getShortTermChange
#' @export
setGeneric("addShortTermChange", signature = "x", function(x, ...)
    standardGeneric("addShortTermChange"))

#' @rdname getShortTermChange
#' @export
setGeneric("getShortTermChange", signature = "x", function(x, ...)
    standardGeneric("getShortTermChange"))

#' @rdname getStability
#' @export
setGeneric("getStability", signature = "x", function(x, ...)
    standardGeneric("getStability"))

#' @rdname getStability
#' @export
setGeneric("addStability", signature = "x", function(x, ...)
    standardGeneric("addStability"))
