test_that("getShortTermChange", {
    library(SummarizedExperiment)
    # Load dataset
    data(minimalgut)
    tse <- minimalgut
    # Check if the function handles empty input
    empty_se <- SummarizedExperiment()
    expect_error(getShortTermChange(empty_se), 
                 "No data available in `x`")

    # Should still return a dataframe
    short_time_labels <- c("74.5h", "173h", "438h", "434h", "390h")
    # Subset samples by Time_label and StudyIdentifier 
    tse_filtered <- tse[, !(tse$Time_label %in% short_time_labels)]
    tse_filtered <- tse_filtered[, (tse_filtered$StudyIdentifier == "Bioreactor A")]
    
    expect_true(all(!(tse_filtered$Time_label %in% short_time_labels)))
    
    result <- getShortTermChange(tse_filtered, time.col = "Time.hr")
    # Expected output is a dataframe
    expect_true(is.data.frame(result))  
    expect_true("growth_diff" %in% colnames(result))
    # Test some expected properties (e.g., that growth_diff isn't all NAs)
    expect_false(all(is.na(result$growth_diff)))
})
