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
        # Get time difference
        args[["time.col"]] <- time.col
        time_res <- do.call(.get_time_difference, args)
        # Create a DF to return
        args <- c(args, list(x_orig = x, res = res, time_res = time_res))
        res <- do.call(.convert_divergence_to_df, args)
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
        # If the data includes repeated timepoints, the data has multiple
        # reference samples for each sample. That is why we make sure that
        # the data is expanded accordingly.
        if( length(ref) != nrow(cd) ){
            x <- x[, names(ref)]
            cd <- cd[names(ref), ]
        } else{
            ref <- ref[ rownames(cd) ]
        }
        cd[[ref.name]] <- unname(ref)
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
    # Add original sample names (if there are duplicated timepoints, certain
    # samples are duplicated). And make column names unique.
    orig_sample_col <- "original_sample_names"
    colData(x)[[orig_sample_col]] <- colnames(x)
    colnames(x) <- make.unique(colnames(x))
    # The returned value includes the TreeSE along with reference
    # column name because it might be that we have modified it.
    res <- list(
        x = x, reference = reference, orig.sample.names = orig_sample_col)
    return(res)
}

# This function returns the first sample for each group by default.
# Alternatively, it returns the previous ith sample for each sample in each
# group.
#' @importFrom dplyr group_by mutate arrange ungroup lag
#' @importFrom tidyr unnest
.get_reference_samples <- function(
        df, time.col, time.interval, group, reference.method){
    # This following line is to suppress "no visible binding for" messages
    # in cmdcheck
    .data <- ":=" <- NULL

    rowname_col <- "temporary_rownames_column"
    reference_col <- "temporary_reference_column"
    # Give warning if the data includes replicated timepoints.
    if( anyDuplicated(df[, c(time.col, group)]) ){
        warning("The data includes duplicated timepoints. Average is applied.",
                call. = FALSE)
    }
    # Add rownames as a column
    df[[rowname_col]] <- rownames(df)
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
            mutate(!!reference_col := .get_previous_samples(
                .data, time.col, rowname_col, time.interval))
    }
    # Ungroup to revert to the original structure, make sure that
    # each cell includes single value (takes action when there are repeated
    # timepoints) and convert to DataFrame
    df <- df |>
        ungroup() |>
        unnest(cols = .data[[reference_col]]) |>
        DataFrame()
    # Get only reference samples
    res <- df[[reference_col]]
    names(res) <- df[[rowname_col]]
    return(res)
}

# This function gets df as input. The data must be already grouped if grouping
# exists. For each sample, this function finds previous timepoint. If there
# are multiple samples from previous timepoint, all are returned.
#' @importFrom dplyr lag
.get_previous_samples <- function(.data, time.col, rowname_col, time.interval){
    # Get timepoints and corresponding previous timepoints
    current_time <- unique(sort(.data[[time.col]]))
    prev_time <- lag(current_time, n = time.interval)
    prev_time <- prev_time[match(.data[[time.col]], current_time)]
    # Split sample names so that they are grouped into timepoints
    timepoint_samples <- split(.data[[rowname_col]], .data[[time.col]])
    # For each sample, assign previous samples
    prev_samples <- prev_samples[ match(prev_time, names(timepoint_samples)) ]
    prev_samples[ lengths(prev_samples) == 0L ] <- NA_character_
    return(prev_samples)
}

# This function get time difference between a sample and its reference sample
.get_time_difference <- function(x, time.col, reference, ...){
    # Get timepoints
    time_point <- x[[time.col]]
    # Get reference time points
    ref <- colData(x)[x[[reference]], time.col]
    # Get difference
    res <- time_point - ref
    return(res)
}

# This function converts time divergence results to DF object
#' @importFrom dplyr summarize
.convert_divergence_to_df <- function(
        x_orig, x, res, time_res, reference, orig.sample.names,
        name = "divergence", name.time = "time_diff", ...){
    #
    temp <- .check_input(name, list("character scalar"))
    temp <- .check_input(name.time, list("character scalar"))
    #
    df <- data.frame(res, time_res, sample = x[[orig.sample.names]])
    # If data includes replicated samples, each time point might have multiple
    # samples. We take average of these repeated time points.
    df <- df |>
        group_by(sample) |>
        summarize(res = mean(res, ...), time_res = mean(time_res, ...)) |>
        DataFrame()
    # Wrangle names
    rownames(df) <- df[["sample"]]
    df[["sample"]] <- NULL
    colnames(df) <- c(name, name.time)
    # Ensure that the order is correct, i.e., matches with the original input
    df <- df[colnames(x_orig), ]
    return(df)
}
