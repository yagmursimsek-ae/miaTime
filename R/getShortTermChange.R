#' @title Short term Changes in Abundance
#'
#' @description Calculates short term changes in abundance of taxa
#'              using temporal Abundance data.
#'  
#' @param x a \code{\link{SummarizedExperiment}} object.
#' @param assay.type \code{Character scalar}. Specifies the name of assay 
#'   used in calculation. (Default: \code{"counts"})
#' @param rarefy \code{Logical scalar}. Whether to rarefy counts.
#' (Default: \code{FALSE})
#' @param compositional \code{Logical scalar}. Whether to transform counts.
#' (Default: \code{FALSE})
#' @param depth \code{Integer scalar}. Specifies the depth used in rarefying. 
#' (Default: \code{min(assay(x, assay.type)}))
#' @param ... additional arguments.
#' 
#' 
#' @return \code{dataframe} with \code{short term change} 
#' calculations.
#' 
#' @details This approach is used by Wisnoski NI and colleagues
#'          \url{https://github.com/nwisnoski/ul-seedbank}. Their approach is based on
#'          the following calculation log(present abundance/past abundance).
#'          Also a compositional version using relative abundance similar to
#'          Brian WJi, Sheth R et al
#'          \url{https://www.nature.com/articles/s41564-020-0685-1} can be used.
#'          This approach is useful for identifying short term growth behaviors of taxa.
#'          
#' @name getShortTermChange
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
#' # Subset samples by Time_lable and StudyIdentifier
#' tse <- tse[, !(tse$Time_label %in% short_time_labels)]
#' tse <- tse[, (tse$StudyIdentifier == "Bioreactor A")]
#' 
#' # Get short term change
#' getShortTermChange(tse, rarefy = TRUE, time.col = "Time.hr")
NULL

#' @rdname getShortTermChange
#' @export
setGeneric("getShortTermChange", signature = c("x"),
    function( x, assay.type = "counts", rarefy = FALSE, compositional = FALSE, 
        depth = min(assay(x, assay.type)), ...)
        standardGeneric("getShortTermChange"))

#' @rdname getShortTermChange
#' @export
#' @importFrom dplyr arrange as_tibble summarize
#' @importFrom ggplot2 ggplot aes
#' @importFrom ggrepel geom_text_repel
#' @importFrom mia rarefyAssay transformAssay
#' @importFrom SummarizedExperiment colData
setMethod("getShortTermChange", signature = c(x = "SummarizedExperiment"),
    function(x, assay.type = "counts", rarefy = FALSE, compositional = FALSE, 
        depth = min(assay(x, assay.type)), ...){
        ############################## Input check #############################
        # Check validity of object
        if(nrow(x) == 0L){
            stop("No data available in `x` ('x' has nrow(x) == 0L.)",
                 call. = FALSE)
        }
        # Check assay.type
        .check_assay_present(assay.type, x)
        if(!.is_a_bool(rarefy)){
            stop("'rarefy' must be TRUE or FALSE.", call. = FALSE)
        }
        if(!.is_a_bool(compositional)){
            stop("'compositional' must be TRUE or FALSE.", call. = FALSE)
        }
        # Ensure that the provided depth is valid
        if ( !is.null(depth) && depth > min(assay(x, assay.type), na.rm = TRUE) ) {
            stop("Depth cannot be greater than the minimum number of counts in your data", call. = FALSE)
        
        }
        ########################### Growth Metrics ############################
        grwt <- .calculate_growth_metrics(x, assay.type = assay.type, 
                                          rarefy = rarefy, 
                                          compositional = compositional, 
                                          depth = depth, ...)
        # Clean and format growth matrics
        grs.all <- .clean_growth_metrics(grwt, ...)
        return(grs.all)
    }
)
# wrapper to calculate growth matrix
.calculate_growth_metrics <- function(x, assay.type, time.col = NULL, 
                                      rarefy, compositional, depth, ...) {
    ############################ Data Preparation ##############################
    # Initialize the filtered object based on rarefy and compositional arguments
    if (rarefy == TRUE && compositional == FALSE) {
        message("rarefy is set to TRUE, calculating short term change using counts")
        x <- rarefyAssay(x, assay.type = assay.type, depth = depth)
        assay.type <- "subsampled"
    } else if (rarefy == FALSE && compositional == FALSE) {
        message("rarefy is set to FALSE, compositional==FALSE, using raw counts")
        x <- x
    } else if (rarefy == FALSE && compositional == TRUE) {
        message("rarefy is set to FALSE, compositional==TRUE, using relative abundances")
        x <- transformAssay(x, method = "relabundance", assay.type = assay.type)
        assay.type <- "relabundance"
    } else if (rarefy == TRUE && compositional == TRUE) {
        stop("Both rarefy and compositional cannot be TRUE simultaneously", call. = FALSE)
    }
    # Reshape data and calcualte grwoth metrics
    assay_data <- meltSE(x, assay.type = assay.type, add.col = time.col)
    assay_data <- assay_data %>%
        arrange( !!sym(time.col) ) %>%
        group_by(FeatureID) %>%
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
        summarize(max.growth = max(growth_diff, na.rm = TRUE))
    # Merge growth data with max growth
    grs.all <- merge(grwt, maxgrs, by = "FeatureID")
    # Add 'ismax' column indicating if the growth is the maximum
    grs.all <- grs.all %>%
        mutate(ismax = ifelse(growth_diff == max.growth, TRUE, FALSE))
    # Clean and abbreviate FeatureID names
    grs.all$FeatureID <- gsub("_", " ", grs.all$FeatureID)
    grs.all$FeatureIDabb <- toupper(abbreviate(grs.all$FeatureID, 
                                               minlength = 3, 
                                               method = "both.sides"))
    # Create 'Feature.time' column combining abbreviation and time information
    grs.all$Feature.time <- paste0(grs.all$FeatureIDabb, " ", 
                                   grs.all[[time.col]], "h")
    
    return(grs.all)
}
