#' @name
#' getShortTermChange
#'
#' @export
#'
#' @title
#' Short term changes in abundance
#'
#' @description
#' Calculates short term changes in abundance of taxa using temporal
#' abundance data.
#'
#' @details
#' This approach is used by Wisnoski NI and colleagues
#' \url{https://github.com/nwisnoski/ul-seedbank}. Their approach is based on
#' the following calculation log(present abundance/past abundance).
#' Also a compositional version using relative abundance similar to
#' Brian WJi, Sheth R et al
#' \url{https://www.nature.com/articles/s41564-020-0685-1} can be used.
#' This approach is useful for identifying short term growth behaviors of taxa.
#'
#' @references
#'
#' Ji, B.W., et al. (2020) Macroecological dynamics of gut microbiota.
#' Nat Microbiol 5, 768–775 . doi: https://doi.org/10.1038/s41564-020-0685-1
#'
#' @return
#' \code{getShortTermChange} returns \code{DataFrame} object containing
#' the short term change in abundance over time for a microbe.
#' \code{addShortTermChange}, on the other hand, returns a
#' \code{\link[SummarizedExperiment:SummarizedExperiment]{SummarizedExperiment}}
#' object with these results in its \code{metadata}.
#'
#' @inheritParams addBaselineDivergence
#'
#' @param name \code{Character scalar}. Specifies a name for storing
#' short term results. (Default: \code{"short_term_change"})
#'
#' @param ... additional arguments.
#
#' @examples
#' library(miaTime)
#'
#' # Load time series data
#' data(minimalgut)
#' tse <- minimalgut
#'
#' # Get relative abundances
#' tse <- transformAssay(tse, method = "relabundance")
#' # Calculate short term changes
#' df <- getShortTermChange(tse, assay.type = "relabundance", time.col = "Time.hr", group = "StudyIdentifier")
#'
#' # Select certain bacteria
#' select <- df[["FeatureID"]] %in% c(
#'     "Ruminococcus_bromii", "Coprococcus_catus", "Akkermansia_muciniphila")
#' df <- df[ select, ]
#'
#' # Plot results
#' library(ggplot2)
#' p <- ggplot(df, aes(x = Time.hr, y = growth_rate, colour = FeatureID)) +
#'     geom_smooth() +
#'     facet_grid(. ~ StudyIdentifier, scales = "free") +
#'     scale_y_log10()
#' p
#'
#' @seealso
#' \code{\link[mia:getStepwiseDivergence]{mia::getStepwiseDivergence()}}
NULL

#' @rdname getShortTermChange
#' @export
setMethod("addShortTermChange", signature = c(x = "SummarizedExperiment"),
    function(x, name = "short_term_change", ...){
        x <- .check_and_get_altExp(x, ...)
        args <- c(list(x = x), list(...)[!names(list(...)) %in% c("altexp")])
        # Calculate short term change
        res <- do.call(getShortTermChange, args)
        # Add to metadata
        x <- .add_values_to_metadata(x, name, res, ...)
        return(x)
    }
)

#' @rdname getShortTermChange
#' @export
setMethod("getShortTermChange", signature = c(x = "SummarizedExperiment"),
    function(
        x, time.col, assay.type = "counts", time.interval = 1L, group = NULL,
        ...){
        ############################## Input check #############################
        x <- .check_and_get_altExp(x, ...)
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
        ########################### Input check end ############################
        res <- .calculate_growth_metrics(
            x, assay.type, time.col, group, time.interval, ...)
        return(res)
    }
)

################################ HELP FUNCTIONS ################################

# wrapper to calculate growth matrix
#' @importFrom dplyr arrange group_by summarise mutate select
.calculate_growth_metrics <- function(
        x, assay.type, time.col, group, time.interval, ...) {
    # Get data in long format
    df <- meltSE(x, assay.type, add.col = c(time.col, group))
    # If there are replicated samples, give warning that average is calculated
    if( anyDuplicated(df[, c("FeatureID", group, time.col)]) ){
        warning("The dataset contains replicated samples. The average is ",
                "calculated for each time point.", call. = FALSE)
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
    df <- df |>
        # Summarize duplicated samples by taking an average
        summarise(value = mean(.data[[assay.type]], na.rm = TRUE)) |>
        # For each feature in a sample group, calculate growth metrics
        mutate(
            time_lag = !!sym(time.col) -
                lag(!!sym(time.col), n = time.interval),
            growth_diff = value - lag(value, n = time.interval),
            growth_rate = (value - lag(value, n = time.interval)) /
                lag(value, n = time.interval),
            var_abund = (value - lag(value, n = time.interval)) / time_lag
        ) |>
        # Remove value column that includes average abundances
        select(-value)
    return(df)
}
