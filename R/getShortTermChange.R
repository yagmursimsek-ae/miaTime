#' @title Short term Changes in Abundance
#'
#' @description Calculates short term changes in abundance of taxa
#'              using temporal Abundance data.
#'  
#' @param x a \code{\link{SummarizedExperiment}} object.
#' @param assay.type \code{Character scalar}. Specifies the name of assay 
#'   used in calculation. (Default: \code{"counts"})
#' @param name \code{Character scalar}. Specifies a name for storing
#' short term results. (Default: \code{"short_term_change"})
#' @param ... additional arguments.
#' 
#' 
#' @return \code{getShortTermChange} returns \code{DataFrame} object containing 
#' the short term change in abundance over time for a microbe. 
#' \code{addShortTermChange}, on the other hand, returns a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object with these results in its \code{metadata}.
#' 
#' @details This approach is used by Wisnoski NI and colleagues
#'          \url{https://github.com/nwisnoski/ul-seedbank}. Their approach is based on
#'          the following calculation log(present abundance/past abundance).
#'          Also a compositional version using relative abundance similar to
#'          Brian WJi, Sheth R et al
#'          \url{https://www.nature.com/articles/s41564-020-0685-1} can be used.
#'          This approach is useful for identifying short term growth behaviors of taxa.
#'          
#' @name addShortTermChange
#' 
#' 
#' @examples
#' 
#' # Load time series data
#' data(minimalgut)
#' tse <- minimalgut
#' 
#' short_time_labels <- c("74.5h", "173h", "438h", "434h", "390h")
#' 
#' # Subset samples by Time_label and StudyIdentifier
#' tse <- tse[, !(tse$Time_label %in% short_time_labels)]
#' tse <- tse[, (tse$StudyIdentifier == "Bioreactor A")]
#' 
#' # Get short term change
#' # Case of rarefying counts
#' tse <- transformAssay(tse, method = "relabundance")
#' getShortTermChange(tse, assay.type = relabundance, time.col = "Time.hr")
#' 
#' # Case of transforming counts
#' tse <- rarefyAssay(tse, assay.type = "counts")
#' getShortTermChange(tse, assay.type = subsampled, time.col = "Time.hr")
NULL

#' @rdname addShortTermChange
#' @export
#' @importFrom dplyr arrange as_tibble summarize "%>%"
#' @importFrom ggplot2 ggplot aes
#' @importFrom ggrepel geom_text_repel
#' @importFrom mia rarefyAssay transformAssay
#' @importFrom SummarizedExperiment colData
setMethod("getShortTermChange", signature = c(x = "SummarizedExperiment"),
    function(x, assay.type = "counts", ...){
        ############################## Input check #############################
        # Check validity of object
        if (nrow(x) == 0L){
            stop("No data available in `x` ('x' has nrow(x) == 0L.)",
                 call. = FALSE)
        }
        # Check assay.type
        .check_assay_present(assay.type, x)
        ########################### Growth Metrics ############################
        grwt <- .calculate_growth_metrics(x, assay.type = assay.type, ...)
        # Clean and format growth metrics
        grs.all <- .clean_growth_metrics(grwt, ...)
        return(grs.all)
    }
)

#' @rdname addShortTermChange
#' @export
setMethod("addShortTermChange", signature = c(x = "SummarizedExperiment"),
          function(x, assay.type = "counts", name = "short_term_change", ...){
              # Calculate short term change
              res <- getShortTermChange(x, ...)
              # Add to metadata
              x <- .add_values_to_metadata(x, name, res, ...)
              return(x)
          }
)

################################ HELP FUNCTIONS ################################

# wrapper to calculate growth matrix
.calculate_growth_metrics <- function(x, assay.type, time.col = NULL, ...) {

    # Reshape data and calculate growth metrics
    assay_data <- meltSE(x, assay.type = assay.type, 
                         add.col = time.col, row.name = "Feature_ID")
    assay_data <- assay_data %>%
        arrange( !!sym(time.col) ) %>%
        group_by(Feature_ID) %>%
        mutate(
            time_lag = !!sym(time.col) - lag( !!sym(time.col) ), 
            growth_diff =!!sym(assay.type) - lag(!!sym(assay.type)),
            growth_rate = (!!sym(assay.type) - lag(!!sym(assay.type))) / lag(!!sym(assay.type)),
            var_abund = (!!sym(assay.type) - lag(!!sym(assay.type))) / time_lag
        )
    return(assay_data)
}

.clean_growth_metrics <- function(grwt, time.col = NULL, ...) {
    # Calculate max growth
    maxgrs <- grwt %>%
        summarize(max_growth = max(growth_diff, na.rm = TRUE))
    # Merge growth data with max growth
    grs.all <- merge(grwt, maxgrs, by = "Feature_ID")
    # Add 'ismax' column indicating if the growth is the maximum
    grs.all <- grs.all %>%
        mutate(ismax = ifelse(growth_diff == max_growth, TRUE, FALSE))
    # Clean and abbreviate FeatureID names
    grs.all$Feature_IDabb <- toupper(abbreviate(grs.all$Feature_ID, 
                                               minlength = 3, 
                                               method = "both.sides"))
    # Create 'Feature.time' column combining abbreviation and time information
    grs.all$Feature_time <- paste0(grs.all$Feature_IDabb, " ", 
                                   grs.all[[time.col]], "h")
    
    return(grs.all)
}
