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
#' tse <- addBaselineDivergence(
#'     tse,
#'     group = "subject",
#'     time.col = "time",
#'     name = "divergence_from_baseline",
#'     name.time = "time_from_baseline",
#'     assay.type="relabundance",
#'     dis.fun = vegan::vegdist,
#'     method="bray")
#'
#' tse <- addBaselineDivergence(
#'     tse,
#'     baseline.sample = "Sample-875",
#'     group = "subject",
#'     time.col = "time",
#'     name = "divergence_from_baseline",
#'     name.time = "time_from_baseline",
#'     assay.type="relabundance",
#'     dis.fun = vegan::vegdist,
#'     method="bray")
#' 
#' @seealso
#' \code{\link[mia:addDivergence]{mia::addDivergence()}}
#' 
NULL

#' @rdname addBaselineDivergence
#' @export
setGeneric("getBaselineDivergence", signature = "x", function(x, ...)
    standardGeneric("getBaselineDivergence"))

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
        name = "divergence",
        name.time = "time_diff",
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
        temp <- .check_input(name, list("character scalar"))
        #
        temp <- .check_input(name.time, list("character scalar"))
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
        res <- .convert_divergence_to_df(x, res, time_res, name, name.time)
        return(res)
    }
)

#' @rdname addBaselineDivergence
#' @export
setGeneric("addBaselineDivergence", signature = "x", function(x, ...)
    standardGeneric("addBaselineDivergence"))

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
        x, time.col, group, reference,
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
    if( is.null(group) ){
        # If it is NULL, add group info --> all samples are in same group
        cd[[group.name]] <- rep("group", nrow(cd))
        group <- group.name
    }
    # If it is a single character value, it should specify a column from
    # colData
    is_wrong_string <- .is_non_empty_character(group) &&
        !group %in% colnames(cd)
    # If it is a vector, then it should have values for all the samples
    is_wrong_vector <- !.is_non_empty_character(group) &&
        length(group) != nrow(cd)
    if( is_wrong_string || is_wrong_vector ){
        stop("'group' must be NULL or a single character value specifying ",
            "a column from colData(x).", call. = FALSE)
    }
    # If it was correctly defined vector, add it to colData
    if( .is_non_empty_character(group) && !group %in% colnames(cd) ){
        cd[[group.name]] <- group
        group <- group.name
    }
    
    # If reference was specified, check that it is specifying samples
    # correctly.
    # It can be a single character value specifying a column from colData
    # (preferred) or single character value specifying a sample.
    is_wrong_string <- FALSE
    if( !is.null(reference) && .is_non_empty_string(reference) ){
        is_wrong_string <- !(reference %in% colnames(cd) ||
            reference %in% rownames(cd))
    }
    # It can also be a character vector. Then its length should match with
    # the length of sample or groups if "group" is specified. (At this point,
    # group cannot be NULL, because we defined it earlier if it was not
    # specified by user)
    is_wrong_vector <- FALSE
    if( !is.null(reference) && !.is_non_empty_string(reference) ){
        is_wrong_vector <- length(reference) != length(unique(cd[[group]]))
        # If the user provided a vector for each group, the vector must be named
        if( !is_wrong_vector && length(reference) != nrow(cd) &&
            is.null(names(reference)) ){
            is_wrong_vector <- TRUE
        }
        # Otherwise, we can expand the reference vector for each member of the
        # groups
        if( !is_wrong_vector && length(reference) != nrow(cd) ){
            reference <- reference[ match(cd[[group]], names(reference)) ]
        }
    }
    if( is_wrong_string || is_wrong_vector ){
        stop("'reference' must be NULL or a single character value specifying ",
            "a column from colData(x).", call. = FALSE)
    }
    # If it was character vector or a sample name, add it to colData
    if( !is.null(reference) && !(.is_non_empty_string(reference) &&
            reference %in% colnames(cd)) ){
        cd[[ref.name]] <- reference
        reference <- ref.name
    }
    
    # If the reference is now NULL, it means that user did not specify it.
    # Get the reference samples.
    if( is.null(reference) ){
        ref <- .get_reference_samples(
            cd, time.col, time.interval, group, reference.method)
        cd[[ref.name]] <- ref
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
.convert_divergence_to_df <- function(x, res, time_res, name, name.time){
    df <- DataFrame(res, time_res, row.names = colnames(x))
    colnames(df) <- c(name, name.time)
    return(df)
}
