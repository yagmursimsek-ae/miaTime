#' @name addMetrics
#' @export
#' 
#' @title
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

#' @rdname addMetrics
#' @export
setMethod("addMetrics", signature = c(x = "SummarizedExperiment"),
    function(
        x, assay.type = "counts", 
        metrics = c("CRT"), thresholds = list(CRT = 2), split.by = NULL, ...){

        ############################## Input checks ############################
        # check assay
        .check_assay_present(assay.type, x)
        # check 'metrics'
        temp <- .check_input( metrics, list("character vector"))
        # Check if 'split.by' exists in colData, otherwise do not split
        if (!is.null(split.by) && !split.by %in% colnames(colData(x))) {
            stop("The 'split.by' column does not exist in colData(x).", 
                 call. = FALSE)
        }

        ########################### Grouped Calculation ########################
        # Extract assay matrix
        assay_data <- assay(x, assay.type)
        
        # Determine grouping structure
        se_list <- if (!is.null(split.by)) {
            splitOn(x, group = split.by)
        } else {
            SimpleList(All = x)
        }
        # Initialize a result data.frame to hold computed metrics
        metric_results <- lapply(se_list, function(sub_se) {
            # Extract assay matrix for the group
            assay_data <- assay(sub_se, assay.type)
            group_results <- list()
            
            # Compute each metric
            for (metric in metrics) {
                if (metric == "CRT") {
                    group_results[[metric]] <- .calculate_crt(
                        assay_data, threshold = thresholds[["CRT"]])
                } else {
                    stop(paste0("Metric '", metric, "' is not supported yet."),
                         call. = FALSE)
                }
            }
            return(group_results)
        })
        browser()
        ######################## Combine and Add to colData ####################
        # Combine results into a single data.frame
        combined_results <- do.call(rbind, 
            lapply(seq_along(metric_results), 
            function(i) {
            group_name <- names(metric_results)[i]
            group_result <- metric_results[[i]]
            data.frame(
                Sample = colnames(se_list[[i]]),
                Group = group_name,
                do.call(cbind, group_result)
            )
        }))
        
        # Sort by sample names
        combined_results <- combined_results[order(combined_results$Sample), ]
        # Add metrics to colData
        for (metric in metrics) {
            x <- .add_values_to_colData(x, 
                combined_results[[metric]], metric, ...)
        }
        
        return(x)
    }
)

# Internal function to calculate Conditional Rare Taxa (CRT)
.calculate_crt <- function(assay_data, threshold = 2) {
    # Calculate the max and min relative abundance for each taxa/feature
    max_abundance <- apply(assay_data, 1, max, na.rm = TRUE)
    min_abundance <- apply(assay_data, 1, min, na.rm = TRUE)
    
    # Compute the max/min ratio
    ratio <- max_abundance / min_abundance
    # Identify taxa where max abundance is greater than or equal to 
    # threshold * min abundance
    crt_features <- ratio >= threshold
    
    return(crt_features)
}
