#' @name
#' getStability
#'
#' @export
#'
#' @title
#'
#' @description
#'
#' @details
#'
#' @references
#'
#' @return
#'
#' @inheritParams addBaselineDivergence
#'
#'
#' @param ... additional arguments.
#' \itemize{
#'   \item \code{time.interval}: \code{Integer scalar}. Indicates the increment
#'   between time steps. By default, the function compares each sample to the
#'   previous one. If you need to take every second, every third, or so, time
#'   step, then increase this accordingly. (Default: \code{1L})
#' }
#
#' @examples
#' library(miaTime)
#'
#' # Load time series data
#' data(minimalgut)
#' tse <- minimalgut
#'
#' @seealso
#' \code{\link[miaTime:getBimodality]{getBimodality()}}
NULL

#' @rdname getStability
#' @export
setMethod("getStability", signature = c(x = "SummarizedExperiment"),
    function(x, time.col, assay.type = "counts", group = NULL, ...){
        ############################## Input check #############################
        x <- .check_and_get_altExp(x, ...)
        temp <- .check_input(
          time.col, list("character scalar"), colnames(colData(x)))
        if( !is.numeric(x[[time.col]]) ){
          stop("'time.col' must specify numeric column from colData(x)",
               call. = FALSE)
        }
        .check_assay_present(assay.type, x)
        temp <- .check_input(
            group, list(NULL, "character scalar"), colnames(colData(x)))
        ########################### Input check end ############################
        # Get data into long format
        df <- meltSE(x, assay.type, add.col = c(time.col, group))
        # Calculate metrics
        df <- .calculate_stability_metrics(
            df, assay.type, time.col, group, ...)
        # Calculate stability based on calculated metrics
        res <- .calculate_stability(df, group, ...)
        # Sort the data
        res <- res[rownames(x), , drop = FALSE]
        return(res)
    }
)

################################ HELP FUNCTIONS ################################

# This function calculates the growth metrics from the data.frame.
#' @importFrom dplyr arrange group_by sym summarise mutate select
.calculate_stability_metrics <- function(
        df, assay.type, time.col, group, time.interval = 1L, ...) {
    # This following line is to suppress "no visible binding for" messages
    # in cmdcheck
    FeatureID <- .data <- value <- NULL
    #
    temp <- .check_input(time.interval, list("numeric scalar"))
    # If there are replicated samples, give warning that average is calculated
    if( anyDuplicated(df[, c("FeatureID", group, time.col)]) ){
        message("The dataset contains replicated samples. The average ",
                "abundance is calculated for each time point.")
    }
    # Sort data based on time
    df <- df |> arrange( !!sym(time.col) )
    # Group based on features and time points. If group was specified, take that
    # also into account
    if( !is.null(group) ){
        df <- df |> group_by(!!sym(group), FeatureID, !!sym(time.col))
    } else{
        df <- df |> group_by(FeatureID, !!sym(time.col))
    }
    # Calculate difference between current and previous time point and
    # difference to the reference point which is the mean value.
    # We do not have to aggregate averages as they are naturally taken into
    # account.
    df <- df |>
        # Summarize duplicated samples by taking an average
        summarise(value = mean(.data[[assay.type]], na.rm = TRUE)) |>
        # Calculate metrics
        mutate(
            diff = value - lag(value, n = time.interval),
            time_diff = .data[[time.col]] -
                lag(.data[[time.col]], n = time.interval),
            start = lag(value, n = time.interval),
            reference = mean(range(value, na.rm = TRUE)),
            reference_diff = start - reference
            ) |>
        ungroup()
    # Filter to remove NAs created by lag
    df <- df |>
        filter_at(vars(diff, time_diff, reference_diff), all_vars(!is.na(.)))
    if( nrow(df) == 0L ){
        stop("No rows left after calculating differences between oncecutive time points.", call. = FALSE)
    }
    return(df)
}

.calculate_stability <- function(df, group, ...){
    if( !is.null(group) ){
        df <- df |> group_by(!!sym(group), FeatureID)
    } else{
        df <- df |> group_by(FeatureID)
    }
    #
    df <- df |> summarise(
        stability = .calc_stability(.data, ...)
    )
    #
    if( !is.null(group) ){
        df <- df |>
            pivot_wider(names_from = group, values_from = stability)
    }
    df <- DataFrame(df, check.names = FALSE)
    rownames(df) <- df[["FeatureID"]]
    df[["FeatureID"]] <- NULL
    return(df)
}

.calc_stability <- function(.data, method = "correlation", ...){
    values <- abs( .data[["diff"]] )
    ref <- abs( .data[["reference_diff"]] )
    time <- .data[["time_diff"]]

    res <- NA
    if( method == "correlation" ){
        # If taxa is very rare, it might be that its standard deviation of
        # 'ref' is 0. Suppress the warning "the standard deviation is zero"
        res <- cor(values, ref) |> suppressWarnings()
    } else if( method == "lm" ){
        model <- lm(values ~ ref + time, data = data.frame(values, ref, time))
        res <- coef(model)
        res <- res[["ref"]]
    }
    return(res)
}
