# Check that input checks work
test_that("getStability input checks work",{
    tse <- makeTSE()
    assayNames(tse) <- "counts"
    tse[["time"]] <- 1L
    #
    res <- getStability(tse) |> expect_error()
    res <- getStability(tse, 1) |> expect_error()
    res <- getStability(tse, "time", group = "test") |> expect_error()
    res <- getStability(tse, "time", reference = "test") |> expect_error()
    res <- getStability(tse, "time", reference = "var1") |> expect_error()
})

# Check that stability is correctly calculated
test_that("getStability calculates correctly",{
    tse <- makeTSE(nrow = 10L, ncol = 200L)
    assay(tse, "counts", withDimnames = FALSE) <- matrix(
        rnorm(10*200), nrow = 10L, ncol = 200L)
    tse[["Time"]] <- sample(seq(1, 5), 200, replace = TRUE)
    # Remove duplicated samples
    tse <- tse[ , !duplicated(colData(tse)[, c("Time")])]
    #
    res <- getStability(tse, "Time") |> expect_no_message()
    expect_s4_class(res, "DataFrame")
    expect_equal(ncol(res), 1L)
    #
    ref <- sapply(rownames(tse), function(feat){
        # Get data for single taxon
        df <- meltSE(tse[feat, ], add.col = "Time")
        # Order data based on time
        df <- df[order(df[["Time"]]), ]
        # Calculate stability metrics
        df[["prev"]] <- c(NA, df[seq_len(nrow(df)-1), "counts"][[1]])
        df[["diff"]] <- df[["counts"]] - df[["prev"]]
        df[["ref"]] <- median(df[["counts"]])
        df[["ref_vs_prev"]] <- df[["prev"]] - df[["ref"]]
        df <- df[seq(2, nrow(df)), ]
        # Calculate stability with correlation
        temp_res <- cor(abs(df[["diff"]]), abs(df[["ref_vs_prev"]]))
        temp_res <- temp_res[[1]]
        return(temp_res)
    })
    expect_equal(res[[1]], unname(ref))
})

# Check that stability is correctly calculated with groups
test_that("getStability calculates correctly with groups",{
    tse <- makeTSE(nrow = 1000, ncol = 200)
    assay(tse, "counts", withDimnames = FALSE) <- matrix(
        rnorm(1000*200), nrow = 1000L, ncol = 200L)
    tse[["Time"]] <- sample(seq(1, 5), 200, replace = TRUE)
    # Remove duplicated samples
    tse <- tse[ , !duplicated(colData(tse)[, c("group", "Time")])]
    #
    res_all <- getStability(tse, "Time", group = "group") |> expect_no_message()
    expect_equal(ncol(res_all), length(unique(tse[["group"]])))
    #
    tse_list <- splitOn(tse, by = 2L, group = "group")
    res <- lapply(tse_list, function(x){
        getStability(x, "Time") |> expect_no_message()
    })
    res <- do.call(cbind, res)
    colnames(res) <- colnames(res_all)
    expect_equal(res_all, res)
})

# Check with duplicates
test_that("getStability calculates correctly with duplicates",{
    tse <- makeTSE(nrow = 10L, ncol = 200L)
    assay(tse, "counts", withDimnames = FALSE) <- matrix(
        rnorm(10*200), nrow = 10L, ncol = 200L)
    tse[["Time"]] <- round(runif(200, 0, 5), 1)
    # Add duplicate samples
    tse[["Time"]][1:20] <- 1
    #
    rowData(tse)[["ref_col"]] <- runif(10L, -10, 10)
    res <- getStability(tse, "Time", reference = "ref_col", mode = "lm") |>
        expect_message()
    #
    ref <- sapply(rownames(tse), function(feat){
        # Get data for single taxon
        df <- meltSE(tse[feat, ], add.col = "Time", add.row = "ref_col")
        # Order data based on time
        df <- df[order(df[["Time"]]), ]
        # Average over time points
        time_points <- unique(df[["Time"]])
        mean_vals <- sapply(time_points, function(time){
            mean(df[df[["Time"]] == time, "counts"][[1]])
        })
        df[["counts"]] <- mean_vals[match(df[["Time"]], time_points)]
        df <- unique(df[, c("counts", "ref_col", "Time")])
        # Calculate stability metrics
        df[["time_diff"]] <- c(NA, diff(df[["Time"]]))
        df[["prev"]] <- c(NA, df[seq_len(nrow(df)-1), "counts"][[1]])
        df[["diff"]] <- df[["counts"]] - df[["prev"]]
        df[["ref_vs_prev"]] <- df[["prev"]] - df[["ref_col"]]
        df <- df[seq(2, nrow(df)), ]
        # Calculate stability with linear model
        temp_res <- lm(abs(diff) ~ abs(ref_vs_prev) + time_diff, df)
        temp_res <- coef(temp_res)[[2]]
        return(temp_res)
    })
    expect_equal(res[[1]], unname(ref))
})
