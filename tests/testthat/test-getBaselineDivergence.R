# Test that the divergence and time difference is correct
test_that("addBaselineDivergence output", {
    data(hitchip1006)
    tse <- hitchip1006
    tse2 <- addBaselineDivergence(
        tse, group = "subject", time.col = "time",
        name = c("divergence_from_baseline", "time_from_baseline",
            "reference_samples"))
    # Input and output classes should match
    expect_equal(class(tse), class(tse2))
    # A subject to check time difference calculation
    time2 <- colData(tse2)[which(tse2[["subject"]] == "843"), "time"]
    time_diff_2 <- colData(tse2)[
        which(tse2[["subject"]] == "843"), "time_from_baseline"]
    expect_true( all(time2 == time_diff_2) )
    # Test divergences
    inds0 <- which(tse2[["subject"]] == "843")
    inds <- which(tse2[["subject"]] == "843")
    original.divergence <- as.matrix(
        vegan::vegdist(t(assay(tse[, inds0], "counts"))))[, 1]
    calculated.divergence <- colData(tse2)[inds, "divergence_from_baseline"]
    expect_true( all(original.divergence == calculated.divergence) )
})

# Test that the result is correct when baseline time point is not 0
test_that("Divergence in baseline other than 0", {
    data(hitchip1006)
    tse <- hitchip1006
    # Should also work when baseline is not 0
    inds <- which(tse[["subject"]] == "843")[2:5]
    tse2 <- addBaselineDivergence(
        tse[, inds], group = "subject", time.col = "time",
        name = c("divergence_from_baseline", "time_from_baseline",
            "reference_samples"))
    time2 <- tse[, inds][["time"]] - min(tse[, inds][["time"]])
    time_diff_2 <- tse2[["time_from_baseline"]]
    expect_true( all(time2 == time_diff_2) )
})

# Test that the reference work
test_that("addBaselineDivergence reference", {
    data(hitchip1006)
    tse <- hitchip1006
    # Just pick 1 subject with many time points
    # The baseline time point 0 is Sample-843
    tse <- tse[, tse[["subject"]] == "843"]
    tse2 <- addBaselineDivergence(tse, group = "subject", time.col = "time")
    # Define the baseline sample manually
    tse3 <- addBaselineDivergence(
        tse, time.col = "time", group = "subject", reference = "Sample-843",
        name = c("divergence_from_baseline", "time_from_baseline",
            "reference_samples"))
    tse4 <- addBaselineDivergence(
        tse, time.col = "time", group = "subject", reference = "Sample-1075",
        name = c("divergence_from_baseline", "time_from_baseline",
            "reference_samples"))
    # Now the times from baseline should be shifted and dissimilarities differ
    # Sample baseline when the zero time baseline is automatically checked or
    # manually set
    expect_true(all(tse2$time_from_baseline==tse3$time_from_baseline))
    # The shifted case (different, middle sample as baseline)
    expect_true(all(tse3$time_from_baseline == tse4$time_from_baseline + 0.7))

    tse5 <- addBaselineDivergence(
        tse[, tse[["subject"]] == "843"], group = "subject",
        time.col = "time",
        name = c("divergence_from_baseline", "time_from_baseline",
            "reference_samples"))
    tse6 <- addBaselineDivergence(
        tse, group = "subject", time.col = "time",
        name = c("divergence_from_baseline", "time_from_baseline",
            "reference_samples"))
    tse7 <- addBaselineDivergence(
        tse, group = "subject", time.col = "time", reference = "Sample-1075",
        name = c("divergence_from_baseline", "time_from_baseline",
            "reference_samples"))
    expect_identical(
        colData(tse5)["Sample-843", "time_from_baseline"],
        colData(tse6)["Sample-843", "time_from_baseline"])
    expect_identical(
        colData(tse5)["Sample-843", "time_from_baseline"] - 0.7,
        colData(tse7)["Sample-843", "time_from_baseline"])

    tse <- hitchip1006
    subjects <- unique(tse$subject)
    # Test with full baseline list
    baselines <- sample(colnames(tse), length(subjects))
    names(baselines) <- subjects
    baselines[names(baselines) == tse[, "Sample-843"][["subject"]]] <-
        "Sample-1075"
    tse8 <- addBaselineDivergence(
        tse, group = "subject", time.col = "time", reference = baselines,
        name = c("divergence_from_baseline", "time_from_baseline",
            "reference_samples"))
    expect_identical(
        colData(tse7)["Sample-843", "time_from_baseline"],
        colData(tse8)["Sample-843", "time_from_baseline"])
    tse[["reference_sam"]] <- baselines[ match(tse$subject, names(baselines)) ]
    res <- addBaselineDivergence(
        tse, group = "subject", time.col = "time", reference = "reference_sam",
        name = c("divergence_from_baseline", "time_from_baseline",
            "reference_samples"))
    ref <- getDivergence(tse, reference = "reference_sam")
    expect_equal(res[["divergence_from_baseline"]], ref)
})

# Test that altExp works
test_that("Test altExp", {
    data(hitchip1006)
    tse <- hitchip1006
    altExp(tse, "Family") <- agglomerateByRank(tse, rank = "Family")
    tse <- addBaselineDivergence(
        tse, group = "subject", time.col = "time", altexp = "Family")
    altExp(tse, "Family_test") <-  addBaselineDivergence(
        altExp(tse, "Family"), group = "subject", time.col = "time",
        name = c("val", "time_val", "ref"))
    # Time differences should still match
    expect_equal(
        altExp(tse, "Family")$divergence, altExp(tse, "Family_test")$val)
})

# Test that get* and add* gives same result
test_that(".get_reference_samples with different time intervals", {
    data(hitchip1006)
    tse <- hitchip1006
    tse <- addBaselineDivergence(
        tse, group = "subject", time.col = "time",
        assay.type = "counts", method = "euclidean")
    res <- getBaselineDivergence(
        tse, group = "subject", time.col = "time",
        assay.type = "counts", method = "euclidean")
    expect_equal(
        colData(tse)[, c("divergence", "time_diff", "ref_samples")], res)
})

# Basic SummarizedExperiment for testing
col_data <- DataFrame(
    time = c(0, 1, 2, 1, 2, 0),
    group = c("A", "A", "A", "B", "B", "B"),
    row.names = c("Sample1", "Sample2",
        "Sample3", "Sample4", "Sample5", "Sample6"))
count_data <- matrix(c(10, 20, 30, 40, 50, 60), ncol = 6, byrow = TRUE)
se <- SummarizedExperiment(assays = list(counts = count_data),
                           colData = col_data)

# Input validation for getBaselineDivergence
test_that("getBaselineDivergence input validations", {
    expect_error(getBaselineDivergence(se, time.col = "nonexistent"))
    expect_error(getBaselineDivergence(se, time.col = "time",
        assay.type = "unknown"))
    expect_error(getBaselineDivergence(se, group = "nonexistent"))
    expect_error(getBaselineDivergence(se, reference = "nonexistent"))
    expect_error(getBaselineDivergence(se, name = "nonexistent"))
    expect_error(getBaselineDivergence(se, name.time = "nonexistent"))
})

# Dissimilarity calculation test
test_that("getBaselineDivergence dissimilarity calculation", {
    result <- getBaselineDivergence(se, time.col = "time", method = "bray")
    expect_s4_class(result, "DataFrame")
    expect_true(all(c("divergence", "time_diff") %in% colnames(result)))
})

# Correct time difference calculation test
test_that("getBaselineDivergence correct time difference calculation", {
    result <- getBaselineDivergence(se, time.col = "time", method = "bray")
    expect_true(all(result$time_diff >= 0))
})

# addBaselineDivergence column addition test
test_that("addBaselineDivergence adds columns to colData", {
    se_result <- addBaselineDivergence(se, time.col = "time", method = "bray")
    expect_true("divergence" %in% colnames(colData(se_result)))
    expect_true("time_diff" %in% colnames(colData(se_result)))
})

# Custom column naming test for addBaselineDivergence
test_that("addBaselineDivergence handles custom column names", {
    se_result <- addBaselineDivergence(
        se, time.col = "time",
        name = c("custom_div", "custom_time_diff",
            "custom_ref_samples"))
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
test_that(".get_reference_samples baseline", {
    stepwise <- .get_reference_samples(
        colData(se), time.col = "time", group = "group",
        reference.method = "stepwise", time.interval = 1)
    expect_equal(stepwise, c(
        NA, "Sample1", "Sample2", "Sample6", "Sample4", NA))
})

# Time difference calculation
test_that(".get_time_difference calculates correct time diff", {
    reference <- c("Sample2", "Sample1", "Sample1", "Sample3", NA, "Sample4")
    se2 <- se
    colData(se2)[["ref"]] <- reference
    time_diffs <- .get_time_difference(
        se2, time.col = "time", reference = "ref")
    expect_equal(time_diffs, c(-1, 1, 2, -1, NA, -1))
})

# Convert divergence to DataFrame
test_that(".convert_divergence_to_df formats correctly", {
    col_data <- DataFrame(
        reference <- c("sample1", "sample2", "sample3", "sample4",
            "sample5", "sample6")
    )
    colnames(col_data) <- c("reference")
    divergence <- c(0.1, 0.2, 0.3, 0, NA, 2)
    time_diff <- c(0, 1, 2, 1, 0, NA)
    se <- SummarizedExperiment(
        assays = list(),
        colData = col_data
    )
    reference <- "reference"
    df <- .convert_divergence_to_df(
        se, divergence, time_diff, reference,
        name = c("test_div", "test_time_diff", "test_reference_samples"))
    expect_s4_class(df, "DataFrame")
    expect_equal(colnames(df),
        c("test_div", "test_time_diff", "test_reference_samples"))
    expect_equal(df$test_div, divergence)
    expect_equal(df$test_time_diff, time_diff)
})

# Test that works with different counts table
test_that("addBaselineDivergence with multiple assay types", {
    assays(se, withDimnames = FALSE) <- list(
        counts = count_data, alt_counts = count_data * 2)
    se_result <- addBaselineDivergence(
        se, time.col = "time", assay.type = "alt_counts")
    expect_true("divergence" %in% colnames(colData(se_result)))
})

# Test that error occurs if if method is unsupported
test_that("getBaselineDivergence unsupported method", {
    expect_error(getBaselineDivergence(
        se, time.col = "time", method = "unsupported"))
})

# Test that the divergence is calculated correctly for specific reference sample
test_that("addBaselineDivergence with custom reference sample", {
    se_result <- addBaselineDivergence(
        se, time.col = "time", reference = "Sample1")
    expect_equal(colData(se_result)["Sample1", "divergence"], 0)
})

# Test that postprocessing works with NA values
test_that(".convert_divergence_to_df with NA divergence values", {
    col_data <- DataFrame(
        reference <- c("sample1", "sample2", "sample3", "sample4",
            "sample5", "sample6")
    )
    divergence <- c(0.1, NA, 0.3, NA, 0.5, 0.6)
    time_diff <- c(0, 1, 2, 1, 0, NA)
    colnames(col_data) <- c("reference")
    se <- SummarizedExperiment(
        assays = list(),
        colData = col_data
    )
    reference <- "reference"
    df <- .convert_divergence_to_df(
        se, divergence, time_diff, reference,
        name = c("test_div", "test_time_diff", "test_reference_samples"))
    expect_s4_class(df, "DataFrame")
    expect_true(all(is.na(df$test_div[is.na(divergence)])))
})
