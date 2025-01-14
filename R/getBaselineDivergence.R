#' @name addBaselineDivergence
#' @export
#' 
#' @title
#' Beta diversity between the baseline and later time steps 
#' 
#' @description
#' Calculates sample dissimilarity between the given baseline and other
#' time points, optionally within a group (subject, reaction chamber, or
#' similar). The corresponding time difference is returned as well.
#' 
#' @details
#' The group argument allows calculating divergence per group. If given, the
#' divergence is calculated per group.  e.g. subject, chamber, group etc.
#' Otherwise, this is done across all samples at once.
#'
#' The baseline sample(s) always need to belong to the data object i.e. they
#' can be merged into it before
#' applying this function. The reason is that they need to have comparable
#' sample data, at least some time point
#' information for calculating time differences w.r.t. baseline sample.
#' 
#' The baseline time point is by default defined as the smallest time point
#' (per group). Alternatively,
#' the user can provide the baseline vector, or a list of baseline vectors per
#' group (named list per group).
#' 
#' @return
#' \code{getBaselineDivergence} returns \code{DataFrame} object
#' containing the sample dissimilarity and corresponding time difference between
#' samples. \code{addBaselineDivergence}, on the other hand, returns a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object with these results in its \code{colData}.
#' 
#' @param x A
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object.
#' 
#' @param assay.type \code{Character scalar}. Specifies which assay values are 
#' used in the dissimilarity estimation. (Default: \code{"counts"})
#' 
#' @param group \code{Character scalar}. Specifies a name of the column from
#' \code{colData} that identifies the grouping of the samples.
#' (Default: \code{NULL})
#' 
#' @param time.col \code{Character scalar}. Specifies a name of the column from
#' \code{colData} that identifies the sampling time points for the samples.
#' 
#' @param method \code{Character scalar}. Used to calculate the dissimilarity 
#' Method is passed to the function that is specified by \code{dis.fun}. 
#' (Default: \code{"bray"})
#' 
#' @param reference \code{Character scalar}. Specifies a name of the column from
#' \code{colData} that identifies the baseline samples to be used.
#' (Default: \code{NULL})
#' 
#' @param name \code{Character scalar}. Specifies a column name for storing
#' divergence results. (Default: \code{"divergence"})
#' 
#' @param name.time \code{Character scalar}. Specifies a column name for storing
#' time differences. (Default: \code{"time_diff"})
#' 
#' @param ... Optional arguments passed into
#' \code{\link[mia:addDivergence]{mia::addDivergence()}}.
#'
#' @examples
#' library(miaTime)
#' library(mia)
#'
#' data(hitchip1006)
#' tse <- transformAssay(hitchip1006, method = "relabundance")
#' 
#' # By default, reference samples are the samples from the first timepoint
#' tse <- addBaselineDivergence(
#'     tse,
#'     group = "subject",
#'     time.col = "time",
#'     assay.type = "relabundance",
#'     method = "bray")
#' 
#' # Add reference samples to colData, if you want to specify reference
#' # samples manually
#' colData(tse)[["reference"]] <- "Sample-875"
#' tse <- addBaselineDivergence(
#'     tse,
#'     reference = "reference",
#'     group = "subject",
#'     time.col = "time",
#'     name = "divergence_from_baseline",
#'     name.time = "time_from_baseline",
#'     assay.type = "relabundance",
#'     method = "bray")
#' 
#' @seealso
#' \code{\link[mia:addDivergence]{mia::addDivergence()}}
#' 
NULL

#' @rdname addBaselineDivergence
#' @export
setMethod("getBaselineDivergence", signature = c(x = "SummarizedExperiment"),
    function(
        x,
        time.col,
        assay.type = "counts",
        reference = NULL,
        group = NULL,
        method = "bray",
        ...){
        ############################# INPUT CHECK ##############################
        x <- .check_and_get_altExp(x, ...)
        # time.col must specify numeric column from colData
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
            reference,
            list(NULL, "character scalar", "character vector"))
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
        # Add baseline samples to colData
        args <- .add_reference_samples_to_coldata(
            x, time.col, group, reference, assay.type, method,
            reference.method = "baseline", ...)
        # Create an argument list. Do not include altexp as it is already taken
        # into account.
        args <- c(
            args,
            list(assay.type = assay.type, method = method),
            list(...)[!names(list(...)) %in% c("altexp")]
        )
        # Calculate divergences
        res <- do.call(getDivergence, args)
        # Add time difference
        x <- args[["x"]]
        reference <- args[["reference"]]
        time_res <- .get_time_difference(x, time.col, reference)
        # Create a DF to return
        res <- .convert_divergence_to_df(x, res, time_res, ...)
        return(res)
    }
)

#' @rdname addBaselineDivergence
#' @export
setMethod("addBaselineDivergence", signature = c(x = "SummarizedExperiment"),
    function(x, name = "divergence", name.time = "time_diff", ...){
        # Calculate divergence
        res <- getBaselineDivergence(x, ...)
        # Add to colData
        res <- as.list(res) |> unname()
        x <- .add_values_to_colData(x, res, list(name, name.time), ...)
        return(x)
    }
)

################################ HELP FUNCTIONS ################################

# This function unifies the input of baseline samples. Despite on how the
# baseline information was provided, this function output TreeSE with baseline
# info for each sample in colData.
.add_reference_samples_to_coldata <- function(
        x, time.col, group, reference = NULL,
        ref.name = "temporal_reference_for_divergence",
        group.name = "temporal_group_for_divergence",
        time.interval = NULL,
        reference.method = "baseline",
        ...){
    #
    temp <- .check_input(
        reference,
        list(NULL, "character scalar", "character vector"))
    #
    temp <- .check_input(ref.name, list("character scalar"))
    #
    temp <- .check_input(group.name, list("character scalar"))
    #
    temp <- .check_input(time.interval, list(NULL, "numeric scalar"))
    #
    temp <- .check_input(
        reference.method, list("character scalar"),
        list("baseline", "stepwise"))
    #
    if( reference.method == "stepwise" && is.null(time.interval) ){
        stop("'time.interval' must be specified.", call. = FALSE)
    }
    # Get colData
    cd <- colData(x)
    
    # Check that group is correctly defined. It can be either NULL, a column
    # from colData or a vector that has group information for all samples.
    # If it is NULL, add group info --> all samples are in same group
    if( is.null(group) ){
        cd[[group.name]] <- rep("group", nrow(cd))
        group <- group.name
    }
    # If it is a single character value, it should specify a column from
    # colData
    is_colname <- .is_non_empty_string(group) && group %in% colnames(cd)
    # If it is a vector, then it should have values for all the samples
    is_vector <- .is_non_empty_character(group) && length(group) == nrow(cd)
    if( !(is_colname || is_vector) ){
        stop("'group' must be NULL or a single character value specifying ",
            "a column from colData(x).", call. = FALSE)
    }
    # If it was correctly defined vector, add it to colData
    if( is_vector ){
        cd[[group.name]] <- group
        group <- group.name
    }
    
    # If the reference is NULL, it means that user did not specify it.
    # Get the reference samples.
    if( is.null(reference) ){
        ref <- .get_reference_samples(
            cd, time.col, time.interval, group, reference.method)
        cd[[ref.name]] <- ref
        reference <- ref.name
    }
    # If reference was specified, check that it is specifying samples
    # correctly.
    # It can be a single character value specifying a column from colData
    # (preferred) or single character value specifying a sample.
    is_colname <- .is_non_empty_string(reference) && reference %in% colnames(cd)
    is_sample <- .is_non_empty_string(reference) && reference %in% rownames(cd)
    # Column name from colData takes precedence
    is_sample <- is_sample && !is_colname
    # It can also be a character vector. Then its length should match with
    # the length of sample or groups if "group" is specified. (At this point,
    # group cannot be NULL, because we defined it earlier if it was not
    # specified by user). Moreover, if the vector specified reference for each
    # group, it must include names that links to groups.
    is_vector_sam <- .is_non_empty_character(reference) &&
        length(reference) == nrow(cd)
    is_vector_group <- .is_non_empty_character(reference) &&
        length(reference) == length(unique(cd[[group]])) &&
        !is.null(names(reference)) && all(names(reference) %in% cd[[group]])
    # Give warning if the input was incorrect
    if( !(is_colname || is_sample || is_vector_sam ||
            is_vector_group) ){
        stop("'reference' must be NULL or a single character value specifying ",
            "a column from colData(x).", call. = FALSE)
    }
    # If the vector was for each group, extend the vector for each sample
    if( is_vector_group ){
        reference <- reference[ match(cd[[group]], names(reference)) ]
    }
    # If it was character vector or a sample name, add it to colData
    if( is_vector_sam || is_vector_group || is_sample ){
        cd[[ref.name]] <- reference
        reference <- ref.name
    }
    
    # Add modified colData back to TreeSE
    colData(x) <- cd
    # The returned value includes the TreeSE along with reference
    # column name because it might be that we have modified it.
    res <- list(x = x, reference = reference)
    return(res)
}

# This function returns the first sample for each group by default.
# Alternatively, it returns the previous ith sample for each sample in each
# group.
#' @importFrom dplyr group_by mutate arrange ungroup lag
.get_reference_samples <- function(
        df, time.col, time.interval, group, reference.method){
    rowname_col <- "temporary_rownames_column"
    reference_col <- "temporary_reference_column"
    # Store rownames and add rownames as a column
    df[[rowname_col]] <- original_order <- rownames(df)
    # Convert to data.frame and group data based on group
    df <- df |>
        as.data.frame() |>
        group_by(.data[[group]])
    
    # Determine the method and perform the respective operations
    if( reference.method == "baseline" ){
        # Find first timepoint within a group
        df <- df |>
            mutate(!!reference_col :=
                .data[[rowname_col]][which.min(.data[[time.col]])])
    } else if( reference.method == "stepwise" ){
        # For each sample, get the previous ith sample.
        # For each subject, get previous sample based on time.
        df <- df |>
            mutate(!!reference_col := lag(
                .data[[rowname_col]], n = time.interval,
                order_by = .data[[time.col]]))
    }
    # Ungroup to revert to the original structure and convert to DataFrame
    df <- df |>
        ungroup() |>
        DataFrame()
    # Put the data into original order
    rownames(df) <- df[[rowname_col]]
    df <- df[original_order, ]
    # Get only reference samples
    res <- df[[reference_col]]
    return(res)
}

# This function get time difference between a sample and its referene sample
.get_time_difference <- function(x, time.col, reference){
    # Get timepoints
    time_point <- x[[time.col]]
    # Get reference time points
    ref <- colData(x)[x[[reference]], time.col]
    # Get difference
    res <- time_point - ref
    return(res)
}

# This function converts time divergence results to DF object
.convert_divergence_to_df <- function(
        x, res, time_res, name = "divergence", name.time = "time_diff", ...){
    #
    temp <- .check_input(name, list("character scalar"))
    #
    temp <- .check_input(name.time, list("character scalar"))
    #
    df <- DataFrame(res, time_res, row.names = colnames(x))
    colnames(df) <- c(name, name.time)
    return(df)
}
