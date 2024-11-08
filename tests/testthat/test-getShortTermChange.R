test_that("getShortTermChange", {
    library(SummarizedExperiment)
    # Load dataset
    data(minimalgut)
    tse <- minimalgut
    # Check if the function handles empty input
    empty_se <- SummarizedExperiment()
    expect_error(getShortTermChange(empty_se), 
                 "No data available in `x`")
    # Check if assay.type argument works
    # tse_invalid <- tse
    # expect_error(
    #     getShortTermChange(tse_invalid, assay.type = "invalid_assay"),
    #     "'assay.type' must be a valid name of assays(x)"
    # )
    # Check that rarefy and compositional cannot both be TRUE
    expect_error(getShortTermChange(tse, rarefy = TRUE, compositional = TRUE, 
                                    time.col = "Time.hr"), 
                 "Both rarefy and compositional cannot be TRUE simultaneously")
    # Check if the depth argument is greater than minimum counts
    min_depth <- min(assay(tse, "counts"))
    expect_error(getShortTermChange(tse, depth = min_depth + 1, 
                                    time.col = "Time.hr"),
                 "Depth cannot be greater than the minimum number of counts in your data")
    # Check if rarefy = TRUE works
    result <- getShortTermChange(tse, rarefy = TRUE, time.col = "Time.hr")
    expect_true(is.data.frame(result))
    # Check if compositional = TRUE works
    result <- getShortTermChange(tse, compositional = TRUE, time.col = "Time.hr")
    expect_true(is.data.frame(result))  
    # Should still return a dataframe
    result <- getShortTermChange(tse, rarefy = TRUE, 
                              compositional = FALSE, 
                              time.col = "Time.hr")
    expect_true(is.data.frame(result))  
    
    short_time_labels <- c("74.5h", "173h", "438h", "434h", "390h")
    # Subset samples by Time_label and StudyIdentifier 
    tse_filtered <- tse[, !(tse$Time_label %in% short_time_labels)]
    tse_filtered <- tse_filtered[, (tse_filtered$StudyIdentifier == "Bioreactor A")]
    
    expect_true(all(!(tse_filtered$Time_label %in% short_time_labels)))
    
    result <- getShortTermChange(tse_filtered, 
                              rarefy = TRUE, time.col = "Time.hr")
    
    result <- getShortTermChange(tse_filtered, 
                              compositional = TRUE, time.col = "Time.hr")
    # Expected output is a dataframe
    expect_true(is.data.frame(result))  
    expect_true("growth_diff" %in% colnames(result))
    # Test some expected properties (e.g., that growth_diff isn't all NAs)
    expect_false(all(is.na(result$growth_diff)))
    
    min_depth <- min(assay(tse_filtered, "counts"))
    result <- getShortTermChange(tse_filtered, rarefy = TRUE, 
                                 depth = min_depth, time.col = "Time.hr")
    expect_true(is.data.frame(result)) 
    expect_error(getShortTermChange(tse_filtered, 
                                 rarefy = TRUE, 
                                 depth = min_depth + 1, time.col = "Time.hr"),
                 "Depth cannot be greater than the minimum number of counts in your data")
    
})
