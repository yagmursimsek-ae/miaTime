# Test: Basic functionality of addStepwiseDivergence
test_that("Basic functionality of addStepwiseDivergence", {
    data(hitchip1006)
    tse <- hitchip1006
    tse2 <- addStepwiseDivergence(
        tse, group = "subject", time.interval = 1, time.col = "time",
        assay.type="counts", dis.fun = vegan::vegdist, method = "bray",
        name = c("divergence", "time_diff", "ref_samples"))
    expect_equal(class(tse), class(tse2))
})

# Test: Adding new colData field with existing name generates warning
test_that("Adding new colData field with existing name generates warning", {
  data(hitchip1006)
  tse <- hitchip1006
  tse[["time_difference"]] <- NA
  expect_warning(addStepwiseDivergence(
      tse, group = "subject", time.interval = 1, time.col = "time",
      name = c("divergence", "time_difference", "ref_samples")))
})

# Test: Time difference calculation for a specific subject
test_that("Time difference calculation is correct for a specific subject", {
  data(hitchip1006)
  tse <- hitchip1006
  tse2 <- addStepwiseDivergence(
      tse, group = "subject", time.interval = 1, time.col = "time",
      assay.type="counts", dis.fun = vegan::vegdist, method = "bray",
      name = c("divergence", "time_difference", "ref_samples"))
  obs_diff <- colData(tse2)[
      which(tse2[["subject"]] == "843"), "time_difference"]
  exp_diff <- c(NA, diff(colData(tse)[
      which(tse[["subject"]] == "843"), "time"]))
  expect_equal(obs_diff, exp_diff)
})

# Test: addStepwiseDivergence with n > 1
test_that("addStepwiseDivergence with n > 1 calculates divergences correctly", {
  data(hitchip1006)
  tse <- hitchip1006
  tse2 <- addStepwiseDivergence(
      tse, group = "subject", time.interval = 2, time.col = "time",
      assay.type = "counts", dis.fun = vegan::vegdist, method = "bray",
      name = c("divergence", "time_difference", "ref_samples"))
  time_interval <- 2
  time <- colData(tse2)[which(tse2[["subject"]] == "843"), "time"]
  time_diff <- colData(tse2)[which(tse2[["subject"]] == "843"),
                             "time_difference"]
  divergence_number <- length(time) - time_interval
  divergence_calculated <- length(which(!is.na(time_diff) == TRUE))
  expect_equal(divergence_number, divergence_calculated)
})

# Test: Interval check for divergence calculation
test_that("Interval check for divergence calculation", {
  data(hitchip1006)
  tse <- hitchip1006
  tse2 <- addStepwiseDivergence(
      tse, group = "subject", time.interval = 2, time.col = "time",
      assay.type = "counts", dis.fun = vegan::vegdist, method = "bray",
      name = c("divergence", "time_difference", "ref_samples"))
  time <- colData(tse2)[which(tse2[["subject"]] == "843"), "time"]
  calculated_diff <- time[(1 + 2):length(time)] -
      time[seq_len(length(time) - 2)]
  manual_diff <- c(rep(
      NA, length(time) - length(calculated_diff)), calculated_diff)
  expect_equal(colData(tse2)[
      which(tse2[["subject"]] == "843"), "time_difference"], manual_diff)
})

# Test: Single time point results in NA divergence values
test_that("Single time point results in NA divergence values", {
  data(hitchip1006)
  tse <- hitchip1006
  tse2 <- tse[, tse[["subject"]] %in% c("900", "843", "139")]
  tse2 <- addStepwiseDivergence(
      tse2, group = "subject", time.interval = 1, time.col = "time",
      assay.type = "counts", dis.fun = vegan::vegdist, method = "bray",
      name = c("time_divergence", "time_diff", "ref_samples"))
  expect_true(all(is.na(colData(tse2)[
      which(duplicated(tse2[["subject"]]) == FALSE),
      "time_divergence"])))
})

# Test: Comparing vegan distances (bray vs euclidean)
test_that("Comparing vegan distances (bray vs euclidean)", {
  data(hitchip1006)
  tse <- hitchip1006
  tse2 <- addStepwiseDivergence(
      tse, group = "subject", time.interval = 1, time.col = "time",
      assay.type = "counts", dis.fun = vegan::vegdist, method = "bray",
      name = c("timedivergence", "timedifference", "ref_samples"))
  tse2 <- addStepwiseDivergence(
      tse2, group = "subject", time.interval = 1, time.col = "time",
      assay.type = "counts", dis.fun = vegan::vegdist, method = "euclidean",
      name = c("timedivergence2", "timedifference2", "ref_samples2"))
  expect_true(identical(tse2$timedifference, tse2$timedifference2))
  expect_true(!identical(tse2$timedivergence, tse2$timedivergence2))
})

# Test: AltExp functionality in addStepwiseDivergence
test_that("AltExp functionality in addStepwiseDivergence", {
  data(hitchip1006)
  tse <- hitchip1006
  altExp(tse, "Family") <- agglomerateByRank(tse, rank = "Family")
  tse <- addStepwiseDivergence(
      tse, group = "subject", time.interval = 1, time.col = "time",
      altexp = "Family")

  altExp(tse, "Family") <- addStepwiseDivergence(
      altExp(tse, "Family"), group = "subject", time.interval = 1,
      time.col = "time",
      name = c("timedivergence", "timedifference", "ref_samples2"))
  expect_equal(
      altExp(tse, "Family")$time_diff,
      altExp(tse, "Family")$timedifference)
  expect_equal(
      altExp(tse, "Family")$divergence,
      altExp(tse, "Family")$timedivergence)
})

# Test: getStepwiseDivergence output type
test_that("getStepwiseDivergence output type", {
  data(hitchip1006)
  tse <- hitchip1006
  divergence_result <- getStepwiseDivergence(
      tse, group = "subject", time.interval = 1, time.col = "time",
      assay.type = "counts", dis.fun = vegan::vegdist, method = "bray")
  expect_s4_class(divergence_result, "DFrame")
  expect_true(all(c("time_diff", "divergence") %in% names(divergence_result)))
})

# Test: Error if time column is missing
test_that("Error if time column is missing", {
  data(hitchip1006)
  tse <- hitchip1006
  expect_error(addStepwiseDivergence(
      tse, group = "subject", time.interval = 1, time.col = "nonexistent_time"))
})

# Test: Error if specified assay type does not exist
test_that("Error if specified assay type does not exist", {
  data(hitchip1006)
  tse <- hitchip1006
  expect_error(addStepwiseDivergence(
      tse, group = "subject", time.interval = 1, time.col = "time",
      assay.type = "nonexistent_assay"))
})

# Test: Error if group column is invalid
test_that("Error if group column is invalid", {
  data(hitchip1006)
  tse <- hitchip1006
  expect_error(addStepwiseDivergence(
      tse, group = "invalid_group", time.interval = 1, time.col = "time"))
})

# Test that time intervals calculation work
test_that(".get_reference_samples with different time intervals", {
    data(hitchip1006)
    tse <- hitchip1006
    interval_1 <- .get_reference_samples(
        colData(tse), time.col = "time", group = "subject",
        reference.method = "stepwise", time.interval = 1) |> unlist()
    interval_2 <- .get_reference_samples(
        colData(tse), time.col = "time", group = "subject",
        reference.method = "stepwise", time.interval = 2) |> unlist()
    expect_false(all(interval_1 == interval_2))
})

# Test that get* and add* gives same result
test_that(".get_reference_samples with different time intervals", {
    data(hitchip1006)
    tse <- hitchip1006
    tse <- addStepwiseDivergence(
        tse, group = "subject", time.interval = 2, time.col = "time",
        assay.type = "counts", method = "euclidean")
    res <- getStepwiseDivergence(
        tse, group = "subject", time.interval = 2, time.col = "time",
        assay.type = "counts", method = "euclidean")
    expect_equal(
        colData(tse)[, c("divergence", "time_diff", "ref_samples")], res)
})

# Basic SummarizedExperiment for testing
col_data <- DataFrame(
    time = c(0, 1, 2, 4, 10, 3),
    group = c("A", "A", "A", "B", "B", "B"),
    row.names = c("Sample1", "Sample2", "Sample3", "Sample4",
        "Sample5", "Sample6"))
count_data <- matrix(c(10, 20, 30, 40, 50, 60), ncol = 6, byrow = TRUE)
se <- SummarizedExperiment(assays = list(counts = count_data),
                           colData = col_data)

# Input validation for getStepwiseDivergence
test_that("getStepwiseDivergence input validations", {
    expect_error(getStepwiseDivergence(se, time.col = "nonexistent"))
    expect_error(getStepwiseDivergence(se, time.col = "time",
        assay.type = "unknown"))
    expect_error(getStepwiseDivergence(se, group = "nonexistent"))
    expect_error(getStepwiseDivergence(se, reference = "nonexistent"))
    expect_error(getStepwiseDivergence(se, name = "nonexistent"))
    expect_error(getStepwiseDivergence(se, name.time = "nonexistent"))
})

# Dissimilarity calculation test
test_that("getStepwiseDivergence dissimilarity calculation", {
    result <- getStepwiseDivergence(se, time.col = "time", method = "bray")
    expect_s4_class(result, "DataFrame")
    expect_true(
        all(c("divergence", "time_diff", "ref_samples") %in%
        colnames(result)))
})

# Correct time difference calculation test
test_that("getStepwiseDivergence correct time difference calculation", {
    result <- getStepwiseDivergence(se, time.col = "time", method = "bray")
    expect_true(any(is.na(result$time_diff)))
})

# addStepwiseDivergence column addition test
test_that("addStepwiseDivergence adds columns to colData", {
    se_result <- addStepwiseDivergence(se, time.col = "time", method = "bray")
    expect_true("divergence" %in% colnames(colData(se_result)))
    expect_true("time_diff" %in% colnames(colData(se_result)))
})

# Custom column naming test for addStepwiseDivergence
test_that("addStepwiseDivergence handles custom column names", {
    se_result <- addStepwiseDivergence(
        se, time.col = "time",
        name = c("custom_div", "custom_time_diff", "ref_samples"))
    # Test
    expect_true("custom_div" %in% colnames(colData(se_result)))
    expect_true("custom_time_diff" %in% colnames(colData(se_result)))
})

# Helper function: assign correct baselines
test_that(".add_reference_samples_to_coldata assigns correct baselines", {
    res <- .add_reference_samples_to_coldata(
        se, time.col = "time", group = "group")
    expect_true(
        "temporal_reference_for_divergence" %in%
        colnames(colData(res[[1]])))
})

# Reference sample assignments
test_that(".get_reference_samples stepwise", {
    stepwise <- .get_reference_samples(
        colData(se), time.col = "time", group = "group",
        reference.method = "stepwise", time.interval = 1)
    expect_equal(stepwise, c(
        NA, "Sample1", "Sample2", "Sample6", "Sample4", NA),
        check.attributes = FALSE)
})


# Test that works with different counts table
test_that("addStepwiseDivergence with multiple assay types", {
    assays(se, withDimnames = FALSE) <- list(
        counts = count_data, alt_counts = count_data * 2)
    se_result <- addStepwiseDivergence(
        se, time.col = "time", assay.type = "alt_counts")
    expect_true("divergence" %in% colnames(colData(se_result)))
})

# Test that error occurs if if method is unsupported
test_that("getStepwiseDivergence unsupported method", {
    expect_error(getStepwiseDivergence(
        se, time.col = "time", method = "unsupported"))
})

# Test that the divergence is calculated correctly for specific reference sample
test_that("addStepwiseDivergence with custom reference sample", {
    res <- getStepwiseDivergence(
        se, time.col = "time", group = "group")
    #
    se[["reference"]] <- c(NA, "Sample1", "Sample2", "Sample6", "Sample4", NA)
    time_diff <- c(NA, 1, 1, 1, 6, NA)
    ref <- getDivergence(se, reference = "reference")
    expect_equal(res[["divergence"]], ref)
    expect_equal(res[["time_diff"]], time_diff)
})

# Test that the function works correctly with replicated samples
test_that("getStepwiseDivergence with replicated samples", {
    tse <- makeTSE(nrow = 1000, ncol = 20)
    assayNames(tse) <- "counts"
    colData(tse)[["time"]] <- sample(c(1, 3, 6, 100), 20, replace = TRUE)
    res <- getStepwiseDivergence(
        tse, time.col = "time", group = "group", method = "euclidean") |>
        expect_warning()
    res <- res[, c(1, 2)]
    #
    sams <- colnames(tse)
    ref <- lapply(sams, function(sam){
        # Get data on sample
        sam_dat <- colData(tse[, sam])
        # Get its reference samples
        group_dat <- colData(tse)[tse[["group"]] == sam_dat[["group"]], ]
        time_points <- unique(sort(group_dat[["time"]]))
        ref_time <- time_points[ which(time_points == sam_dat[["time"]])-1 ]
        temp_res <- c(NA, NA)
        # Calculate distance if reference samples were found
        if( length(ref_time) > 0 ){
            # Loop through each reference sample, calculate its distance to
            # the sample, and take mean
            ref_sams <- rownames(group_dat[ group_dat[["time"]] == ref_time, ])
            ref_vals <- vapply(ref_sams, function(ref_sam){
                dist(t( assay(tse, "counts")[, c(sam, ref_sam)] ),
                    method = "euclidean")
            }, numeric(1))
            ref_vals <- mean(ref_vals)
            # Return divergence and time difference
            temp_res <- c(ref_vals, sam_dat[["time"]] - ref_time)
        }
        return(temp_res)
    })
    ref <- do.call(rbind, ref)
    ref <- DataFrame(ref)
    colnames(ref) <- colnames(res)
    rownames(ref) <- rownames(res)
    expect_equal(res, ref)
})
