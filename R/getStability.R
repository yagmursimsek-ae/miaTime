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
    function(x, time.col, assay.type = "counts", reference = NULL, group = NULL,
        ...){
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
        # The reference must be NULL or column name from rowData. Additionally,
        # it can be numeric vector but then user must specify reference value
        # for each feature (or global reference).
        temp <- .check_input(
            reference,
            list(NULL, "character scalar", "numeric vector", "list"),
            supported_values = colnames(rowData(x)),
            length = c(1L, nrow(x))
            )
        ########################### Input check end ############################
        # Get data into long format
        df <- meltSE(x, assay.type, add.col = c(time.col, group))
        # If there are duplicated samples, calculate average
        df <- .summarize_duplicates(
            df, assay.type, "FeatureID", time.col, group)
        # Add reference values to df
        args <- .add_reference_for_stability(
            df, x, assay.type, reference, group)
        # Calculate metrics for stability calculation
        args <- c(args,
            list(assay.type = assay.type, time.col = time.col, group = group),
            list(...))
        df <- do.call(.calculate_stability_metrics, args)
        # Calculate stability based on calculated metrics
        res <- .calculate_stability(df, group, ...)
        # Sort the data so that it matches with original order
        res <- res[rownames(x), , drop = FALSE]
        return(res)
    }
)

################################ HELP FUNCTIONS ################################

# The dataset can include duplicated samples, i.e., multiple samples for
# each patient in certain time point. This function can be utilized to summarize
# these samples so that for each feature, there is only one data point for
# each patient, time point and feature.
#' @importFrom dplyr group_by across all_of mutate ungroup distinct
.summarize_duplicates <- function(
        df, assay.type, feature.id = NULL, time.col = NULL, group = NULL){
    # If there are replicated samples, give warning that average is calculated
    if( anyDuplicated(df[, c(feature.id, group, time.col)]) ){
        message("The dataset contains replicated samples. The average ",
                "abundance is calculated for each time point.")
        # Summarize duplicated samples by taking an average
        df <- df |>
            # Take unique set based on sample group, feature and time point
            group_by(across(all_of(c(group, feature.id, time.col)))) |>
            # Summarize them
            mutate(!!assay.type := mean(.data[[assay.type]], na.rm = TRUE)) |>
            ungroup() |>
            # Keep only unique rows
            distinct(across(all_of(c(group, feature.id, time.col))),
                .keep_all = TRUE)
    }
    return(df)
}

# This function adds reference values to the data.
.add_reference_for_stability <- function(df, x, assay.type, reference, group){
    ref_name <- "ref_values"
    # If user did not specify reference, get default
    if( is.null(reference) ){
        df <- .default_reference_for_stability(df, assay.type, group, ref_name)
    } else{
        # If user specified reference values, add them to df
        df <- .custom_reference_for_stability(df, x, reference, group, ref_name)
    }
    args <- list(df = df, reference = ref_name)
    return(args)
}

# This function calculate the default reference values and adds them to data.
#' @importFrom dplyr group_by across all_of mutate ungroup
.default_reference_for_stability <- function(df, assay.type, group, ref.name){
    # If reference was not provided, calculate default reference values,
    # i.e, median for each feature and system.
    df <- df |> group_by(across(all_of(c(group, "FeatureID")))) |>
        mutate(!!ref.name := median(.data[[assay.type]], na.rm = TRUE)) |>
        ungroup()
    return(df)
}

# This function adds user-specified reference values to the data.
.custom_reference_for_stability <- function(df, x, reference, group, ref.name){
    # If reference is a list, group must be specified
    if( is.list(reference) && is.null(group) ){
        stop("'reference' cannot be a list if 'group' is not specified. ",
            "Please use simple vector.", call. = FALSE)
    }
    # If user gave a list as reference, it must include reference for each
    # system. To map refernce values correctly, also system names must be
    # included. Moreover, the values must be numeric
    if( is.list(reference) &&
            !(all(lengths(reference) == length(unique((df[[group]])))) &&
            all(x[[group]] %in% unique(names(unlist(reference)))) &&
            all(is.numeric(unlist(reference)))
            ) ){
        stop("If 'reference' is provided as a list, all values must be ",
            "numeric. Additionally, each group must have a corresponding ",
            "reference value, and these values must be named to match the ",
            "groups specified in 'group'.", call. = FALSE)
    }

    # If reference is just single numeric value, expand it for all features
    if( .is_a_numeric(reference) ){
        reference <- rep(reference, nrow(x))
    }
    # If reference specifies a field from rowData, get the values
    if( .is_non_empty_string(reference) ){
        reference <- rowData(x)[[reference]]
    }
    # If the reference is a list, create a data.frame that can be matched with
    # the data. The data.frame will include rownames, groups and reference
    # values
    if( is.list(reference) ){
        names(reference) <- rownames(x)
        reference <- reference |> as.data.frame() |> t() |> stack() |>
            as.data.frame()
        colnames(reference) <- c("FeatureID", group, ref.name)
    } else{
        # If reference was a vector, convert it to data.frame that can be
        # added to data
        reference <- setNames(
            data.frame(rownames(x), reference) , c("FeatureID", ref.name))
    }
    # Add references to the data
    df <- dplyr::left_join(
        df, reference, by = intersect(colnames(df), colnames(reference)))
    return(df)
}

# This function calculates the metrics that will be later used for calculating
# of stability.
#' @importFrom dplyr arrange group_by across all_of mutate lag filter_at vars
#'     all_vars
.calculate_stability_metrics <- function(
        df, assay.type, time.col, reference, group, time.interval = 1L, ...){
    #
    temp <- .check_input(time.interval, list("numeric scalar"))
    # Sort data based on time
    df <- df |> arrange( !!sym(time.col) )
    # Calculate metrics for stability
    df <- df |>
        # Take unique set based on sample group, feature and time point
        group_by(across(all_of(c(group, "FeatureID")))) |>
        mutate(
            # Difference between current and previous time point
            curr_vs_prev = .data[[assay.type]] -
                lag(.data[[assay.type]], n = time.interval),
            # Difference between current time point and reference point
            prev_vs_ref = lag(.data[[assay.type]], n = time.interval) -
                .data[[reference]],
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
            "with 'time.interval=", time.interval, "'.", call. = FALSE)
    }
    return(df)
}

# This function facilitates the actual calculation of stability based on the
# previously calculated metrics.
#' @importFrom dplyr group_by across all_of mutate ungroup distinct select
#'     any_of ungroup
#' @importFrom tidyr pivot_wider
.calculate_stability <- function(df, group, calc.separately = FALSE, ...){
    temp <- .check_input(calc.separately, list("logical scalar"))
    # Group based on feature so that we calculate stability for each feature.
    # If sample group was also specified, apply also it so that we calculate
    # stability for each feature-patient/system pair.
    df <- df |> group_by(across(all_of(c(group, "FeatureID"))))
    # Calculate stability. We use help function to apply and control how the
    # stability is calculated.
    df <- df |> mutate(
        stability = .calc_stability(.data, ...)
    ) |> ungroup()

    # If user wants to calculate separately stability for values that are
    # less or greater than reference
    if( calc.separately ){
        # Create column that tells if the data point is larger than reference
        df[["larger_than_ref"]] <- ifelse(
            df[["prev_vs_ref"]] > 0, "stability_right", "stability_left")
        # For each group and feature, calculate left and right stability
        # separately
        df <- df |>
            group_by(across(all_of(c(group, "FeatureID", "larger_than_ref"))))
        df <- df |> mutate(
            stability_sep = .calc_stability(.data, ...)
        ) |> ungroup()
    }
    # Get taxa names and stability values
    df <- df |>
        select(any_of(c(
            "FeatureID", group, "stability", "larger_than_ref",
            "stability_sep"))) |> distinct()
    # If we also calculated stability separately for left and right, currently
    # they are in long format. Modify the data so that left and right values
    # are in own columns
    if( calc.separately ){
        df <- df |> pivot_wider(
            names_from = larger_than_ref, values_from = stability_sep)
    }

    # If we had sample groups, we put each sample group to own columns.
    # Otherwise we have single column which contain stability values for each
    # feature.
    if( !is.null(group) ){
        df <-  df |>
            group_by(across(all_of(c(group, "FeatureID")))) |>
            pivot_wider(
                names_from = group,
                values_from = any_of(c(
                    "stability", "stability_left", "stability_right"))
            ) |> ungroup()
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
