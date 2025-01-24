#' @name addStepwiseDivergence
#' @export
#'
#' @title
#' Beta diversity between consecutive time steps
#'
#' @description
#' Calculates sample dissimilarity between consecutive time points along with
#' time difference.
#'
#' @details
#' These functions calculate time-wise divergence, meaning each sample is
#' compared to the previous i-th sample, where i is the specified time
#' interval (see \code{time.interval}). By default, the function calculates
#' divergence by comparing all samples with each other. However, it is often
#' more meaningful to calculate divergence within a specific patient or group
#' (see the \code{group} parameter).
#'
#' @return
#' \code{getStepwiseDivergence} returns \code{DataFrame} object
#' containing the sample dissimilarity and corresponding time difference between
#' samples. \code{addStepwiseDivergence}, on the other hand, returns a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object with these results in its \code{colData}.
#'
#' @inheritParams addBaselineDivergence
#'
#' @param time.interval \code{Integer scalar}. Indicates the increment between
#' time steps. By default, the function compares each sample to the
#' previous one. If you need to take every second, every third, or so, time
#' step, then increase this accordingly. (Default: \code{1L})
#
#' @examples
#' library(miaTime)
#'
#' data(hitchip1006)
#' tse <- transformAssay(hitchip1006, method = "relabundance")
#'
#' # Calculate divergence
#' tse <- addStepwiseDivergence(
#'     tse,
#'     group = "subject",
#'     time.interval = 1,
#'     time.col = "time",
#'     assay.type = "relabundance"
#'     )
#'
#' @seealso
#' \code{\link[mia:addDivergence]{mia::addDivergence()}}
#'
NULL

#' @rdname addStepwiseDivergence
#' @export
setMethod("getStepwiseDivergence", signature = c(x = "ANY"),
    function(
        x,
        time.col,
        assay.type = "counts",
        time.interval = 1L,
        group = NULL,
        method = "bray",
        ...){
        ############################# INPUT CHECK ##############################
        temp <- .check_input(
            time.col, list("character scalar"), colnames(colData(x)))
        if( !is.numeric(x[[time.col]]) ){
            stop("'time.col' must specify numeric column from colData(x)",
                call. = FALSE)
        }
        #
        .check_assay_present(assay.type, x)
        #
        temp <- .check_input(
            group, list(NULL, "character scalar"), colnames(colData(x)))
        #
        temp <- .check_input(method, list("character scalar"))
        #
        if( is.null(rownames(x)) ){
            rownames(x) <- paste0("row", seq_len(nrow(x)))
        }
        if( is.null(colnames(x)) ){
            colnames(x) <- paste0("col", seq_len(ncol(x)))
        }
        ########################### INPUT CHECK END ############################
        x <- .check_and_get_altExp(x, ...)
        # Add stepwise samples to colData
        args <- .add_reference_samples_to_coldata(
            x, time.col, group, time.interval = time.interval,
            reference.method = "stepwise", ...)
        # Create an argument list. Do not include altexp as it is already taken
        # into account.
        args <- c(
            args,
            list(assay.type = assay.type, method = method),
            list(...)[!names(list(...)) %in% c("altexp")]
        )
        # Calculate divergences
        res <- do.call(getDivergence, args)
        # Get time difference
        args[["time.col"]] <- time.col
        time_res <- do.call(.get_time_difference, args)
        # Create a DF to return
        args <- c(args, list(res = res, time_res = time_res))
        res <- do.call(.convert_divergence_to_df, args)
        return(res)
    }
)

#' @rdname addStepwiseDivergence
#' @export
setMethod("addStepwiseDivergence", signature = c(x = "SummarizedExperiment"),
    function(
        x, name = c("divergence", "time_diff", "ref_samples"), ...){
        .check_input(name, c("character vector"), length = 3L)
        # Calculate divergence
        res <- getStepwiseDivergence(x,  ...)
        # Add to colData
        res <- as.list(res) |> unname()
        x <- .add_values_to_colData(x, res, name, ...)
        return(x)
    }
)
