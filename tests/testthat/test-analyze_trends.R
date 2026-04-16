context("Core Extraction Logics & Statistical Robustness")

test_that("run_mk extracts consistently and .mmkh_get parses fragile outputs reliably", {
    # Test internal helper directly 
    mock_vec1 <- c("new P-value" = 0.04)
    mock_vec2 <- c("new P.value" = 0.05)
    mock_vec3 <- c("Corrected P-value" = 0.01)
    
    expect_equal(EffOff:::.mmkh_get(mock_vec1, "new P[- .]?value|Corrected P"), 0.04)
    expect_equal(EffOff:::.mmkh_get(mock_vec2, "new P[- .]?value|Corrected P"), 0.05)
    expect_equal(EffOff:::.mmkh_get(mock_vec3, "new P[- .]?value|Corrected P"), 0.01)
    
    # run_mk fallback execution testing smoothly bounding against small lists unconditionally 
    res_short <- EffOff:::run_mk(c(1, 2))
    expect_true(is.na(res_short$pVal))
    expect_equal(res_short$Method, NA_character_)
})

test_that("analyze_trends throws explicit errors unconditionally elegantly on multiple basenames natively executing limits", {
    bad_df <- data.frame(
        SysTime = Sys.Date(),
        PixID = c(1, 1, 2, 2),
        Year = c(2000, 2001, 2000, 2001),
        Basename = c("PlotA", "PlotA", "PlotB", "PlotB")
    )
    
    expect_error(EffOff:::analyze_trends(bad_df), "evi_data contains multiple Basename values")
})

test_that("analyze_trends always unconditionally gracefully explicitly explicitly securely returns polymorphic list", {
    # Simple explicit bounding array mapping conditionally smoothly checking execution implicitly effectively effectively reliably safely
    good_df <- data.frame(
        PixID = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
        latitude = rep(20, 10), longitude = rep(80, 10),
        Basename = "PlotA", Year = 2000:2009, EVI = c(0.2, 0.3, 0.4, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
    )
    
    old_opts <- options(mc.cores = 1)
    old_opts <- options(mc.cores = 1)
    res <- EffOff:::analyze_trends(good_df)
    expect_type(res, "list")
    expect_true("Overall" %in% names(res))
    
    # Bounding safely tracking inherently effectively checking half-open ranges
    res_multi <- EffOff:::analyze_trends(good_df, project_start_year = 2005, project_end_year = 2008)
    expect_type(res_multi, "list")
    expect_true("Pre_Management" %in% names(res_multi))
    expect_true("Management" %in% names(res_multi))
    
    # Check boundaries
    df_pre <- res_multi$Pre_Management
    df_mgmt <- res_multi$Management
    expect_false(grepl("2005", df_pre$AnalysisDuration[1]))
    expect_true(grepl("2005", df_mgmt$AnalysisDuration[1]))
    options(old_opts)
})
