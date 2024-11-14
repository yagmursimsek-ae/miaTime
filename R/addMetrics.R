#' Add Multiple Metrics (e.g., CRT) to a \code{\link{SummarizedExperiment}} object
#' 
#' @param x a \code{\link{SummarizedExperiment}} object.
#' @param assay.type \code{Character scalar}. Specifies the name of the assay 
#'   used in calculation. (Default: \code{"counts"})
#' @param metrics \code{Character vector}. Specifies the metrics to be 
#'   calculated. For example, \code{CRT}. (Default: \code{c("CRT")})
#' @param thresholds \code{Named list}. A list of thresholds for each metric. 
#'   For example, \code{list(CRT = 2)} specifies the threshold for CRT 
#'   calculation. (Default: \code{list(CRT = 2)})
#' @param split.by \code{Character scalar}. Specifies how to split the data 
#'   for per-group or per-subject calculation. For example, \code{"subject"}. 
#'   If \code{NULL}, no splitting will occur. (Default: \code{NULL})
#' @param ... Additional arguments to pass to individual metric functions. 
#' @return A \code{\link{SummarizedExperiment}} object with added metadata 
#'   containing the calculated metrics.
#' 
NULL

#' @export
setGeneric("addMetrics", function(
    x, assay.type = "counts", 
    metrics = c("CRT"), thresholds = list(CRT = 2), split.by = NULL, ...) {
    standardGeneric("addMetrics")
})

#' @rdname addMetrics
#' @export
setMethod("addMetrics", signature = c(x = "SummarizedExperiment"),
    function(
        x, assay.type = "counts", 
        metrics = c("CRT"), thresholds = list(CRT = 2), split.by = NULL, ...){

        ############################## Input checks ############################
        # Validate 'metrics'
        if (!.is_non_empty_character(metrics)) {
            stop("'metrics' should be a non-empty character vector.", 
                 call. = FALSE)
        }
        # Check if 'split.by' exists in colData, otherwise do not split
        if (!is.null(split.by) && !split.by %in% colnames(colData(x))) {
            stop("The 'split.by' column does not exist in colData(x).", 
                 call. = FALSE)
        }

        ############################## Calculate metrics #######################
        # Split the data by subject (or other grouping variable) if required
        if (!is.null(split.by)) {
            split_data <- split(x, f = colData(x)[[split.by]])
        } else {
            split_data <- list(x)  
        }

        # Calculate metrics for each split subject/group
        results_list <- lapply(split_data, function(sub_data) {
            # Initialize a list to hold results for each metric
            metric_results <- list()
            
            # Loop through each metric and calculate it
            for (metric in metrics) {
                # Check if the metric is supported and calculate accordingly
                if (metric == "CRT") {
                    crt_results <- .calculateCRT(sub_data, thresholds[["CRT"]])
                    metric_results[["CRT"]] <- crt_results
                } else {
                    stop(paste("Metric", metric, "is not supported yet."), 
                         call. = FALSE)
                }
            }

            # Add the calculated metrics to the metadata (colData or rowData)
            for (metric in names(metric_results)) {
                colData(sub_data)[[metric]] <- metric_results[[metric]]
            }
            
            return(sub_data)
        })

        # add results to metadata
        x <- do.call(cbind, results_list)
        
        return(x)
    }
)

# Internal function to calculate Conditional Rare Taxa (CRT)
.calculateCRT <- function(data, threshold = 2) {
    # Extract the counts/abundance data
    counts <- assay(data)
    # Calculate the max and min relative abundance for each taxa/feature
    max_abundance <- apply(counts, 1, max)
    min_abundance <- apply(counts, 1, min)
    
    # Compute the max/min ratio
    ratio <- max_abundance / min_abundance
    # Identify taxa where max abundance is greater than or equal to 
    # threshold * min abundance
    crt_features <- ratio >= threshold
    
    return(crt_features)
}
