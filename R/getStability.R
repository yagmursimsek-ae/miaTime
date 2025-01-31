#' @name
#' getStability
#'
#' @export
#'
#' @title
#' Calculate stability
#'
#' @description
#' ADD.
#'
#' @details
#' The stability is calculated by first defining reference point, \eqn{R_{f}},
#' for each feature. This is defined by taking mean of range of values
#'
#' \deqn{R_{f} = \frac{x_{f, min} + x_{f, max}}{2}}
#'
#' where \eqn{f} denotes a single feature and \eqn{x} abundance.
#'
#' Then difference between consecutive time points, \eqn{\Delta_{f, x}},
#'
#' \deqn{\Delta_{f, x} = ∣x_{f, t} - x_{f, t-1}∣}
#'
#' and difference between the previous time point and reference,
#' \eqn{\Delta_{f. R}}, are calculated:
#'
#' \deqn{\Delta_{f. R} = ∣x_{f, t-1} - R_{f}∣}
#'
#' where \eqn{t} denotes time point.
#'
#' Optionally, time difference \eqn{\Delta_{f, t}} is calculated
#'
#' \deqn{\Delta_{f, t} = t_{f, t-1} - t_{f, t-1}}
#'
#' The stability coefficient \eqn{s_{f}} is calculated with correlation
#'
#' \deqn{s_{f} = corr(\Delta_{f, x}, \Delta_{f, R})}
#'
#' or with linear model as follows
#'
#' \deqn{s_{f} = \frac{ \Delta_{f, x} - \beta_{0} - \beta_{2}*\Delta_{f, t} -
#' \epsilon_{f} }{ \Delta_{f, R} }}
#'
#' @references
#' Lahti L, et al. (2014) Tipping elements in the human intestinal ecosystem.
#' Nat Commun. doi: 10.1038/ncomms5344
#'
#' @return
#' Add here
#'
#' @inheritParams addBaselineDivergence
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
#' # Apply clr transformation
#' tse <- transformAssay(tse, method = "rclr")
#'
#' # Calculate stability single system
#' tse_sub <- tse[, tse[["StudyIdentifier"]] == "Bioreactor A"]
#' res <- getStability(tse_sub, assay.type = "rclr", time.col = "Time.hr")
#' res |> head()
#'
#' # Calculate stability for each  system simultaneously by taking time
#' # difference into account
#' res <- getStability(
#'     tse, assay.type = "rclr", time.col = "Time.hr",
#'     group = "StudyIdentifier", mode = "lm")
#' res |> head()
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
        # Calculate metrics for stability calculation
        df <- .calculate_stability_metrics(
            df, assay.type, time.col, group, ...)
        # Calculate stability based on calculated metrics
        res <- .calculate_stability(df, group, ...)
        # Sort the data so that it matches with original order
        res <- res[rownames(x), , drop = FALSE]
        return(res)
    }
)

################################ HELP FUNCTIONS ################################

# This function calculates the metrics that will be later used for calculating
# of stability.
#' @importFrom dplyr arrange group_by sym summarise mutate select
.calculate_stability_metrics <- function(
        df, assay.type, time.col, group, time.interval = 1L, ...) {
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
    # Calculate metrics for stability
    df <- df |>
        # Summarize duplicated samples by taking an average
        summarise(value = mean(.data[[assay.type]], na.rm = TRUE)) |>
        # Calculate metrics
        mutate(
            # Difference between current and previous time point
            curr_vs_prev = value - lag(value, n = time.interval),
            # Difference between current time point and feature's median
            prev_vs_ref = lag(value, n = time.interval) -
                mean(range(value, na.rm = TRUE)),
            # Time difference between consecutive time points
            time_diff = .data[[time.col]] -
                lag(.data[[time.col]], n = time.interval)
            ) |>
        ungroup()
    # Filter to remove NAs created by lag
    df <- df |> filter_at(
        vars(curr_vs_prev, prev_vs_ref, time_diff), all_vars(!is.na(.)))
    # Check that there are rows left. It might be that the lag is too large
    # for the dataset. For example, if there is only one time point, we cannot
    # calculate the difference between consecutive time points.
    if( nrow(df) == 0L ){
        stop("Cannot calculate the difference between consecutive time points ",
            "with 'ẗime.interval=", time.interval, "'.", call. = FALSE)
    }
    return(df)
}

# This function facilitates the actual calculation of stability based on the
# previously calculated metrics.
.calculate_stability <- function(df, group, ...){
    # Group based on feature so that we calculate stability for each feature.
    # If sample group was also specified, apply also it so that we calculate
    # stability for4 each feature-patient/system pair.
    if( !is.null(group) ){
        df <- df |> group_by(!!sym(group), FeatureID)
    } else{
        df <- df |> group_by(FeatureID)
    }
    # Calculate stability. We use help function to apply and control how the
    # stability is calculated.
    df <- df |> summarise(
        stability = .calc_stability(.data, ...)
    )
    # If we had sample groups, we put each sample group to own columns.
    # Otherwise we have single column which contain stability values for each
    # feature.
    if( !is.null(group) ){
        df <- df |>
            pivot_wider(names_from = group, values_from = stability)
    }
    # Convert tibble to DF and wrangle rownames.
    df <- DataFrame(df, check.names = FALSE)
    rownames(df) <- df[["FeatureID"]]
    df[["FeatureID"]] <- NULL
    return(df)
}

# This function is used to control how the stability is calculated. As an input
# we get values for single system and feature.
.calc_stability <- function(.data, mode = "correlation", ...){
    .check_input(mode, list("character scalar"), list("correlation", "lm"))
    # Get metrics
    values <- abs( .data[["curr_vs_prev"]] )
    ref <- abs( .data[["prev_vs_ref"]] )
    time <- .data[["time_diff"]]

    res <- NA
    if( mode == "correlation" ){
        # Apply simple correlation that do not take into account time
        # difference. Get the correlation between "difference between
        # consecutive time points" and "difference between previous time point
        # and reference point".
        # If taxa is very rare, it might be that its standard deviation of
        # 'ref' is 0. Suppress the warning "the standard deviation is zero"
        res <- cor(values, ref, ...) |> suppressWarnings()
    } else if( mode == "lm" ){
        # Apply linear model that takes into account time difference between
        # consecutive time points. Get how much the "difference between previous
        # time point and reference point" explains the "difference between
        # consecutive time points".
        model <- lm(values ~ ref + time, data = data.frame(values, ref, time))
        res <- coef(model)
        res <- res[["ref"]]
    }
    # The result is a single value that tells how much the distance from
    # median explains the difference between consecutive time points.
    # If the value is positive, we can say that the system is stable:
    # the abundance of taxon is predictable; with larger distance from
    # reference the changes are larger and vice versa.
    # If the value is negative, we can say that the system is unstable:
    # the magnitude fluctuations are negatively linked with the abundance level;
    # with lower abundance the changes are larger.
    return(res)
}
