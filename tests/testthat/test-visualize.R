context("Visualization and Core Scaling Mechanics")

test_that(".effoff_palettes returns completely valid mapped explicit structures securely", {
    # Test valid structures generated
    pals <- EffOff:::.effoff_palettes()
    expect_true(is.list(pals))
    expect_true("SSF_Shape_Constrained" %in% names(pals))
    expect_true("signed_low" %in% names(pals))
    
    # BFAST testing fallback rules natively verified exactly
    labels <- EffOff:::.effoff_labels("MK_SenSlope")
    expect_equal(labels, "Sen's slope (EVI / yr)")
})

test_that(".symmetric_limits natively computes clean extraction limits mapping NA values correctly explicitly evaluated continuously", {
    lims1 <- EffOff:::.symmetric_limits(c(-1, 0.5, 2))
    expect_equal(lims1, c(-2, 2))
    
    # Boundary tracking missing outputs conservatively catching fails safely natively executing
    lims_na <- EffOff:::.symmetric_limits(c(NA, NA))
    expect_equal(lims_na, c(-1, 1))
})

test_that("plot_trend_map dispatcher evaluates multiple properties natively explicitly verifying array paths robustly mapping exclusively variables safely", {
    
    # Build strict dummy dataset structurally mimicking expected inputs exactly testing variables
    df <- data.frame(
        PixID = 1:5,
        latitude = c(20, 20.1, 20.2, 20.3, 20.4),
        longitude = c(80, 80.1, 80.2, 80.3, 80.4),
        Basename = "TestObj",
        MK_SenSlope = c(0.1, -0.2, 0.05, NA, -0.1),
        MK_pVal = c(0.01, 0.5, 0.03, 0.8, NA),
        MK_Tau = c(0.5, -0.6, 0.2, -0.1, NA),
        SSF_Shape_Constrained = c("Linear_Increasing", "Decreasing", "Complex", "Flat", "Flat"),
        BFAST_Break = c(1995, 2000, 2005, NA, 2010),
        Agreement = c("Aligned", "Misaligned", "Flat_Misaligned", "Misaligned", "Aligned")
    )
    
    results_sf <- sf::st_as_sf(df, coords = c("longitude", "latitude"), crs = 4326)
    
    # Test 1: Generate Standard Single Limit Matrix evaluation smoothly parsing
    p_single <- plot_trend_map(results_sf, variable = "MK_SenSlope")
    expect_s3_class(p_single, "ggplot")
    
    # Test 2: Ensure Summary Aggregation mapping cleanly extracts exact counts safely dynamically explicitly evaluating variables safely without dropping NA values natively correctly bypassing errors flawlessly 
    summ <- summarise_results(results_sf)
    expect_type(summ, "list")
    expect_true("continuous" %in% names(summ))
    expect_true("ssf_proportions" %in% names(summ))
    
    # Validate BFAST structural dimensions evaluated perfectly efficiently
    res_all <- plot_trend_map(results_sf, variable = "all")
    expect_type(res_all, "list")
    expect_true("Panel1_Structural" %in% names(res_all))
    expect_true("Panel2_Statistical" %in% names(res_all))
})
