# Tests for .calculate_skewness (references calculated with moments::skewness)
test_that(".calculate_skewness works as expected", {
    # Symmetrical data
    expect_equal(.calculate_skewness(c(1, 2, 3)), 0)
    expect_equal(.calculate_skewness(rnorm(10000)), 0, tolerance = 0.05)
    # Positively skewed data
    expect_equal(
        round(.calculate_skewness(c(1, 1, 1, 10, 10, 10, 10, 10, 10, 10)), 2),
        -0.87)
    # Negatively skewed data
    expect_equal(
        round(.calculate_skewness(c(10, 10, 10, 10, 1, 1, 1)), 2), -0.29)
    # Constant data
    expect_equal(.calculate_skewness(rep(5, 10)), NaN)
})

# Tests for .calculate_kurtosis (references calculated with moments::kurtosis)
test_that(".calculate_kurtosis works as expected", {
    # Normalish data
    expect_equal(
        round(.calculate_kurtosis(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)), 2), 1.78)
    # High kurtosis data
    expect_equal(
        round(.calculate_kurtosis(c(1, 1, 1, 1, 50, 50, 50, 50, 50, 50)), 2),
        1.17)
    # Constant data should result in NaN
    expect_equal(.calculate_kurtosis(rep(5, 10)), NaN)
})

# Tests for .calculate_bimodality_coefficient (references calculated with
# mousetrap::bimodality_coefficient)
test_that(".calculate_bimodality_coefficient works as expected", {
    # Normalish data
    dat <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
    skew <- .calculate_skewness(dat)
    kurt <- .calculate_kurtosis(dat)
    ref <- (1+skew^2)/(kurt+3)
    expect_equal(.calculate_bimodality_coefficient(dat), ref)
    # Constant data should result in NaN
    expect_equal(.calculate_bimodality_coefficient(rep(5, 10)), NaN)
})

# Create dummy data
tse <- makeTSE(nrow = 100, ncol = 20)
assayNames(tse) <- "counts"

# Tests for .calculate_bimodality
test_that(".calculate_bimodality works as expected", {
    result <- .calculate_bimodality(tse, assay.type = "counts", group = "group")
    expect_true(is.data.frame(result) || inherits(result, "DataFrame"))
    expect_equal(ncol(result), length(unique(tse$group)))
    expect_equal(nrow(result), nrow(tse))
})

# Tests for getBimodality
test_that("getBimodality works as expected", {
    result <- getBimodality(tse, assay.type = "counts")
    expect_true(is.data.frame(result) || inherits(result, "DataFrame"))
    expect_equal(ncol(result), 1)
    expect_equal(nrow(result), nrow(tse))
    ref <- .calculate_bimodality_coefficient(assay(tse, "counts")[1, ])
    expect_equal(result[[1]][[1]], ref)
})

# Tests for addBimodality
test_that("addBimodality works as expected", {
    tse <- addBimodality(tse, assay.type = "counts", name = "testing")
    expect_true("testing" %in% colnames(rowData(tse)))
    result <- getBimodality(tse, assay.type = "counts", name = "testing")
    expect_equal(rowData(tse)[, "testing", drop = FALSE], result)
})

# Tests that error messages work
test_that("addBimodality: errors", {
    addBimodality(tse, assay.type = "test") |> expect_error()
    addBimodality(tse, group = "test") |> expect_error()
    addBimodality(tse, name = 1) |> expect_error()
    addBimodality(tse, name = c("test", "test2")) |> expect_error()
    addBimodality(tse, name = TRUE) |> expect_error()
})
