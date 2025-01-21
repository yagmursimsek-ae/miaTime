#' @name
#' getBimodality
#'
#' @title
#' Calculate coefficient of bimodality.
#'
#' @description
#' This function calculates coefficient of bimodality for each taxa.
#'
#' @details
#' This function calculates coefficient of bimodality for each taxa. If the
#' dataset includes grouping, for instance, individual systems or patients,
#' the coefficient is calculated for each group separately.
#'
#' The coefficient of bimodality measures whether a taxon has bimodal abundance.
#' For instance, certain taxon can be high-abundant in some timepoints while
#' in others it might be rare. The coefficient can help to determine these
#' taxa.
#'
#' The coefficient of bimodality (b) is defined as follows:
#'
#' \deqn{b = \frac{1+skewness^{2}}{kurtosis+3}}
#'
#' where skewness is calculated as follows
#'
#' \deqn{skewness = \frac{\sum_{i=1}^{n}(x_{i}-\overline{x})^{3}/n}{
#' \sum_{i=1}^{n}(x_{i}-\overline{x})^{2}/n]^{3/2}}}
#'
#' and kurtosis as follows
#'
#' \deqn{kurtosis = \frac{\sum_{i=1}^{n}(x_{i}-\overline{x})^{4}/n}{
#' (\sum_{i=1}^{n}(x_{i}-\overline{x})^{2}/n)^{2}}}
#'
#' The coefficient ranges from 0-1, where 1 means high bimodality. The
#' coefficient was introduced in the paper Shade A et al. 2014, where they used
#' bimodality and abundance to determine conditionally rare taxa (CRT).
#'
#' @references
#'
#' Shade A, et al. (2014)
#' Conditionally Rare Taxa Disproportionately Contribute to Temporal Changes in
#' Microbial Diversity. doi: 10.1128/mbio.01371-14
#'
#' @return
#' \code{DataFrame} or \code{x} with results added to its \code{rowData}.
#'
#' @inheritParams addBaselineDivergence
#'
#' @param name \code{Character scalar}. Specifies a column name for storing
#' bimodality results. (Default: \code{"bimodality"})
#'
#' @param ... additional arguments.
#' \itemize{
#'   \item \code{group}: \code{Character scalar}. Specifies a name of the column
#'   from \code{colData} that identifies the grouping of the samples.
#'   (Default: \code{NULL})
#' }
#'
#' @examples
#' library(miaTime)
#' library(mia)
#'
#' data(SilvermanAGutData)
#' tse <- SilvermanAGutData
#'
#' # In this example, we are only interested on vessel 1.
#' tse <- tse[, tse[["Vessel"]] == 1]
#' tse <- transformAssay(tse, method = "relabundance")
#'
#' # Calculate bimodality
#' b <- getBimodality(tse)
#' # Determine taxa with high bimodality
#' bimodal_taxa <- names(b)[ which(b[[1]] > 0.95) ]
#'
#' # Determine taxa with abundance > 0.5%
#' abundant <- getAbundant(
#'     tse, assay.type = "relabundance", abundant.th = 0.005)
#'
#' # The detected CRT
#' crt <- intersect(bimodal_taxa, abundant)
#' head(crt)
#'
#' @seealso
#' \code{\link[mia:getAbundant]{mia::getConditionallyLowAbundant()}}
#'
NULL

#' @rdname getBimodality
#' @export
setMethod("addBimodality", signature = c(x = "SummarizedExperiment"),
    function(x, name = "bimodality", ...){
        .check_input(name, "character scalar")
        x <- .check_and_get_altExp(x, ...)
        res <- getBimodality(x, name = name, ...) |> as.list()
        x <- .add_values_to_colData(x, unname(res), names(res), MARGIN = 1L)
        return(x)
    }
)

#' @rdname getBimodality
#' @export
setMethod("getBimodality", signature = c(x = "SummarizedExperiment"),
    function(x, assay.type = "counts", ...){
        x <- .check_and_get_altExp(x, ...)
        .check_assay_present(assay.type, x)
        res <- .calculate_bimodality(x, assay.type, ...)
        return(res)
    }
)

################################ HELP FUNCTIONS ################################

# This function calculates bimodality for each group for eacg feature.
#' @importFrom dplyr group_by across all_of summarize
#' @importFrom tidyr pivot_wider
.calculate_bimodality <- function(
        tse, assay.type, name = "bimodality", group = NULL, ...){
    .check_input(name, "character scalar")
    .check_input(group, c("character scalar", "NULL"), colnames(colData(tse)))
    # Get data as long format
    df <- meltSE(tse, assay.type = assay.type, add.col = group)
    # Calculate bimodality for each group and feature
    res <- df |>
        group_by(across(all_of(c(group, "FeatureID")))) |>
        summarize(
            b = .calculate_bimodality_coefficient(.data[[assay.type]], ...),
            .groups = "drop"
            )
    # If grouping was used, put data to wide format so that columns include
    # groups and rows features.
    if( !is.null(group) ){
        res <- res |>
            pivot_wider(
                names_from = group,
                values_from = b
            )
    }
    # Put into DataFrame format and wrangle rownames
    res <- DataFrame(res, check.names = FALSE)
    rownames(res) <- res[["FeatureID"]]
    res[["FeatureID"]] <- NULL
    # Make sure that the order of data reflects the original data.
    res <- res[match(rownames(tse), rownames(res)), , drop = FALSE]
    # Rename columns by user specified name
    if( ncol(res) == 1L ){
        colnames(res) <- name
    } else{
        colnames(res) <- paste0(name, "_", colnames(res))
    }
    return(res)
}

# Function to calculate bimodality coefficient for single vector/taxon
.calculate_bimodality_coefficient <- function(x, ...) {
    skewness <- .calculate_skewness(x, ...)
    kurtosis <- .calculate_kurtosis(x, ...)
    b <- (1 + skewness^2) / (kurtosis + 3)
    return(b)
}

# Function to calculate skewness for single vector/taxon
.calculate_skewness <- function(x, ...) {
    mean_x <- mean(x, ...)
    n <- length(x)
    numerator <- sum((x - mean_x)^3) / n
    denominator <- (sum((x - mean_x)^2) / n)^(3/2)
    skewness <- numerator / denominator
    return(skewness)
}

# Function to calculate kurtosis for single vector/taxon
.calculate_kurtosis <- function(x, ...) {
    mean_x <- mean(x, ...)
    n <- length(x)
    numerator <- sum((x - mean_x)^4) / n
    denominator <- (sum((x - mean_x)^2) / n)^2
    kurtosis <- numerator / denominator
    return(kurtosis)
}
