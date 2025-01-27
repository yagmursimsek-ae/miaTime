
test_that("getShortTermChange errors", {
    tse <- makeTSE(nrow = 100, ncol = 20)
    assayNames(tse) <- "counts"
    tse[["Time"]] <- sample(seq(1, 5), 20, replace = TRUE)

    df <- getShortTermChange(tse, time.col = "Test", assay.type = "counts") |>
        expect_error()
    df <- getShortTermChange(tse, time.col = "group", assay.type = "counts") |>
        expect_error()
    df <- getShortTermChange(tse, time.col = "group", assay.type = "test") |>
        expect_error()
    df <- getShortTermChange(tse, time.col = "Time", group = "counts") |>
        expect_error()
    df <- getShortTermChange(tse, time.col = "Time", time.interval = TRUE) |>
        expect_error()
    tse <- addShortTermChange(tse, time.col = "Time", name = TRUE) |>
        expect_error()
})

# Test with time interval 1 and no grouping
test_that("getShortTermChange calculates correctly without group", {
    tse <- makeTSE(nrow = 100, ncol = 20)
    assayNames(tse) <- "counts"
    tse[["Time"]] <- sample(seq(1, 5), 20, replace = TRUE)

    df <- getShortTermChange(tse, time.col = "Time", assay.type = "counts") |>
        expect_message()

    # Check the structure of the returned object
    expect_s4_class(df, "DataFrame")
    expect_true("time_diff" %in% colnames(df))
    expect_true("abundance_diff" %in% colnames(df))
    expect_true("growth_rate" %in% colnames(df))
    expect_true("rate_of_change" %in% colnames(df))
    expect_true("Time" %in% colnames(df))
    expect_true("FeatureID" %in% colnames(df))

    # Test 10 times with random combinations
    for( i in seq_len(10) ){
        # Get random feature and timepoint to test
        feat <- sample(rownames(tse), 1)
        time_point <- sample(tse[["Time"]], 1)
        # Get index of previous sample
        temp_df <- colData(tse)[ order(tse[["Time"]]), ]
        time_points <- unique(temp_df[["Time"]])
        prev_time <- which(time_points == time_point)-1
        ref <- rep(NA_real_, 4)
        # Calculate reference values manually if previous sample exists
        if( prev_time > 0 ){
            prev_time <- time_points[ prev_time ]
            prev_sam <- rownames(temp_df[temp_df[["Time"]] == prev_time, ])
            sam <- rownames(temp_df[temp_df[["Time"]] == time_point, ])
            # Take certain samples
            tse_curr <- tse[, sam]
            tse_prev <- tse[, prev_sam]
            # Get values
            time_diff <- time_point - prev_time
            val1 <- mean(assay(tse_curr, "counts")[feat, ], na.rm = TRUE)
            val2 <- mean(assay(tse_prev, "counts")[feat, ], na.rm = TRUE)
            ref <- c(
                time_diff,
                val1 - val2,
                (val1 - val2) / val2,
                (val1 - val2) / time_diff
            )
        }
        names(ref) <- c(
            "time_diff", "abundance_diff", "growth_rate", "rate_of_change")
        res <- df[df[["FeatureID"]] == feat & df[["Time"]] == time_point, ]
        res <- unlist(res[, names(ref)])
        #
        expect_equal(res, ref)
    }
})

# Test with time interval 2 and no grouping
test_that("getShortTermChange calculates correctly with time interval 2", {
    tse <- makeTSE(nrow = 1000, ncol = 200)
    assayNames(tse) <- "counts"
    tse[["Time"]] <- sample(seq(1, 5), 200, replace = TRUE)
    tse <- transformAssay(tse, method = "relabundance")
    df <- getShortTermChange(
        tse, time.col = "Time", assay.type = "relabundance",
        time.interval = 2) |>
        expect_message()

    # Check the structure of the returned object
    expect_s4_class(df, "DataFrame")
    expect_true("time_diff" %in% colnames(df))
    expect_true("abundance_diff" %in% colnames(df))
    expect_true("growth_rate" %in% colnames(df))
    expect_true("rate_of_change" %in% colnames(df))
    expect_true("Time" %in% colnames(df))
    expect_true("FeatureID" %in% colnames(df))

    # Test 10 times with random combinations
    for( i in seq_len(10) ){
        # Get random feature and timepoint to test
        feat <- sample(rownames(tse), 1)
        time_point <- sample(tse[["Time"]], 1)
        # Get index of previous sample
        temp_df <- colData(tse)[ order(tse[["Time"]]), ]
        time_points <- unique(temp_df[["Time"]])
        prev_time <- which(time_points == time_point)-2
        ref <- rep(NA_real_, 4)
        # Calculate reference values manually if previous sample exists
        if( prev_time > 0 ){
            prev_time <- time_points[ prev_time ]
            prev_sam <- rownames(temp_df[temp_df[["Time"]] == prev_time, ])
            sam <- rownames(temp_df[temp_df[["Time"]] == time_point, ])
            # Take certain samples
            tse_curr <- tse[, sam]
            tse_prev <- tse[, prev_sam]
            # Get values
            time_diff <- time_point - prev_time
            val1 <- mean(assay(tse_curr, "relabundance")[feat, ], na.rm = TRUE)
            val2 <- mean(assay(tse_prev, "relabundance")[feat, ], na.rm = TRUE)
            ref <- c(
                time_diff,
                val1 - val2,
                (val1 - val2) / val2,
                (val1 - val2) / time_diff
            )
        }
        names(ref) <- c(
            "time_diff", "abundance_diff", "growth_rate", "rate_of_change")
        res <- df[df[["FeatureID"]] == feat & df[["Time"]] == time_point, ]
        res <- unlist(res[, names(ref)])
        #
        expect_equal(res, ref)
    }
})

# Test with grouping
test_that("getShortTermChange calculates correctly with group", {
    tse <- makeTSE(nrow = 1000, ncol = 200)
    assayNames(tse) <- "counts"
    tse[["Time"]] <- sample(seq(1, 5), 200, replace = TRUE)
    # Remove duplicated samples
    tse <- tse[ , !duplicated(colData(tse)[, c("group", "Time")])]

    df <- getShortTermChange(
        tse, time.col = "Time", assay.type = "counts", group = "group") |>
        expect_no_message()

    # Check the structure of the returned object
    expect_s4_class(df, "DataFrame")
    expect_true("time_diff" %in% colnames(df))
    expect_true("abundance_diff" %in% colnames(df))
    expect_true("growth_rate" %in% colnames(df))
    expect_true("rate_of_change" %in% colnames(df))
    expect_true("Time" %in% colnames(df))
    expect_true("FeatureID" %in% colnames(df))

    # Test 10 times with random combinations
    for( i in seq_len(10) ){
        # Get random feature and timepoint to test
        feat <- sample(rownames(tse), 1)
        time_point <- sample(tse[["Time"]], 1)
        group <- sample(tse[["group"]], 1)
        # Get only specific group
        temp_df <- colData(tse)[ tse[["group"]] == group, ]
        # Get index of previous sample
        temp_df <- temp_df[ order(temp_df[["Time"]]), ]
        time_points <- unique(temp_df[["Time"]])
        prev_time <- which(time_points == time_point)-1
        ref <- rep(NA_real_, 4)
        # Calculate reference values manually if previous sample exists
        if( prev_time > 0 ){
            prev_time <- time_points[ prev_time ]
            prev_sam <- rownames(temp_df[temp_df[["Time"]] == prev_time, ])
            sam <- rownames(temp_df[temp_df[["Time"]] == time_point, ])
            # Take certain samples
            tse_curr <- tse[, sam]
            tse_prev <- tse[, prev_sam]
            # Get values
            time_diff <- time_point - prev_time
            val1 <- mean(assay(tse_curr, "counts")[feat, ], na.rm = TRUE)
            val2 <- mean(assay(tse_prev, "counts")[feat, ], na.rm = TRUE)
            ref <- c(
                time_diff,
                val1 - val2,
                (val1 - val2) / val2,
                (val1 - val2) / time_diff
            )
        }
        names(ref) <- c(
            "time_diff", "abundance_diff", "growth_rate", "rate_of_change")
        select <- df[["FeatureID"]] == feat & df[["Time"]] == time_point &
            df[["group"]] == group
        res <- df[select, ]
        res <- unlist(res[, names(ref)])
        #
        expect_equal(res, ref)
    }
})

# Test that getShortTimeChange and addShortTimeChange are equal
test_that("get* and add* are equal", {
    tse <- makeTSE(nrow = 1000, ncol = 200)
    assayNames(tse) <- "counts"
    tse[["Time"]] <- sample(seq(1, 5), 200, replace = TRUE)
    # Remove duplicated samples
    tse <- tse[ , !duplicated(colData(tse)[, c("Time")])]

    df <- getShortTermChange(tse, time.col = "Time") |>
        expect_no_message()
    tse <- addShortTermChange(tse, time.col = "Time", name = "test") |>
        expect_no_message()

    expect_equal(metadata(tse)[["test"]], df)
})
