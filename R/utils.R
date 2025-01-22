################################### TESTING ###################################
# Methods for testing

# This function unifies input testing. The message will always be in same format
# also it makes the code simpler in main function since testing is done here.
# Borrowed from HoloFoodR.
.check_input <- function(
        variable, supported_class, supported_values = NULL, limits = NULL,
        variable_name = .get_name_in_parent(variable)){
    # Convert supported classes to character
    classes_char <- lapply(supported_class, function(class){
        if( is.null(class) ){
            class <- "NULL"
        }
        return(class)
    })
    classes_char <- unlist(classes_char)
    # Based on number of acceptable classes, the msg is different
    class_txt <- .create_msg_from_list(classes_char)
    # Create a message
    msg <- paste0("'", variable_name, "' must be ", class_txt, "." )

    # If supported values were provided
    if( !is.null(supported_values) ){
        # Convert supported values to character
        values_char <- lapply(supported_values, function(value){
            if( is.null(value) ){
                value <- "NULL"
            }
            value <- as.character(value)
            return(value)
        })
        values_char <- unlist(values_char)
        # Collapse into text
        values_txt <- paste0("'", paste(values_char, collapse = "', '"), "'")
        msg <- paste0(
            msg, " It must be one of the following options: ", values_txt)
    }

    # If limits were provided
    if( !is.null(limits) ){
        msg <- paste0(msg, " (Numeric constrains: ")
        # Add thresholds to message
        if( !is.null(limits$upper) ){
            msg <- paste0(msg, limits$upper, ">x")
        } else if(!is.null(limits$upper_include)){
            msg <- paste0(msg, limits$upper, ">=x")
        }
        if( !is.null(limits$lower) ){
            msg <- paste0(msg, "x>", limits$lower)
        } else if(!is.null(limits$lower_include)){
            msg <- paste0(msg, "x>=", limits$lower_include)
        }
        msg <- paste0(msg, ")")
    }

    # List all the input types. Run the check if the variable must be that type.
    # If correct type was found, change the result to TRUE.
    input_correct <- FALSE
    if( "NULL" %in% classes_char && is.null(variable) ){
        input_correct <- TRUE
    }
    if( "logical scalar" %in% classes_char && .is_a_bool(variable) ){
        input_correct <- TRUE
    }
    if( "logical vector" %in% classes_char && is.logical(variable) ){
        input_correct <- TRUE
    }
    if( "character scalar" %in% classes_char && .is_non_empty_string(
        variable) ){
        input_correct <- TRUE
    }
    if( "character vector" %in% classes_char && .is_non_empty_character(
        variable) ){
        input_correct <- TRUE
    }
    if( "numeric scalar" %in% classes_char && .is_a_numeric(variable) ){
        input_correct <- TRUE
    }
    if( "numeric vector" %in% classes_char && is.numeric(variable) ){
        input_correct <- TRUE
    }
    if( "integer vector" %in% classes_char && .is_integer(variable) ){
        input_correct <- TRUE
    }
    if( "integer scalar" %in% classes_char && .is_an_integer(variable) ){
        input_correct <- TRUE
    }
    if( "list" %in% classes_char && is.list(variable) && !is.data.frame(
        variable) ){
        input_correct <- TRUE
    }
    if( "data.frame" %in% classes_char && is.data.frame(variable) ){
        input_correct <- TRUE
    }
    if( "matrix" %in% classes_char && is.matrix(variable) ){
        input_correct <- TRUE
    }
    # If supported values were provided
    if( !is.null(supported_values) && !is.null(variable) ){
        # Test that if variable is in supported values
        values_correct <- lapply(supported_values, function(value){
            res <- FALSE
            if( is.null(value) && is.null(variable) || value %in% variable){
                res <- TRUE
            }
            return(res)
        })
        values_correct <- unlist(values_correct)
        # If not, then give FALSE even though class checks were correct
        if( !any(values_correct) ){
            input_correct <- FALSE
        }
    }
    # If limits were provided
    if( !is.null(limits) && !is.null(variable) ){
        if( !is.null(limits$upper) && variable >= limits$upper ){
            input_correct <- FALSE
        } else if( !is.null(
            limits$upper_include) && variable > limits$upper_include ){
            input_correct <- FALSE
        }

        if( !is.null(limits$lower) && variable <= limits$lower ){
            input_correct <- FALSE
        } else if( !is.null(
            limits$upper_include) && variable < limits$upper_include ){
            input_correct <- FALSE
        }
    }
    # Give error if variable was not correct type
    if( !input_correct ){
        stop(msg, call. = FALSE)
    }
    return(input_correct)
}

# This function creates a string from character values provided. The string
# can be used to messages. It creates a tidy list from list of values.
.create_msg_from_list <- function(classes_char, and_or = "or", ...){
    if( length(classes_char) > 2 ){
        class_txt <- paste0(
            paste(
                classes_char[seq_len(length(classes_char)-1)], collapse = ", "),
            " ", and_or, " ", classes_char[length(classes_char)])
    } else if( length(classes_char) == 2 ){
        class_txt <- paste0(
            classes_char[[1]], " ", and_or, " ", classes_char[[2]])
    } else{
        class_txt <- classes_char
    }
    return(class_txt)
}

#################### INTERNAL METHODS FROM EXTERNAL PACKAGES ###################
# internal methods loaded from other packages

.is_a_bool <- mia:::.is_a_bool
.is_non_empty_character <- mia:::.is_non_empty_character
.is_non_empty_string <- mia:::.is_non_empty_string
.is_integer <- mia:::.is_integer
.is_an_integer <- mia:::.is_an_integer
.is_a_numeric <- mia:::.is_a_numeric
.get_name_in_parent <- mia:::.get_name_in_parent
.safe_deparse <- mia:::.safe_deparse
.check_altExp_present <- mia:::.check_altExp_present
.check_assay_present <- mia:::.check_assay_present
.add_values_to_colData <- mia:::.add_values_to_colData
.check_and_get_altExp <- mia:::.check_and_get_altExp
