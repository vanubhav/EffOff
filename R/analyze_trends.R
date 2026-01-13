#' Run BFAST Analysis
#'
#' Applies BFAST Lite to detect breakpoints in the time series.
#'
#' @param ts_data A time series object (ts).
#' @return The year of the first detected breakpoint, or NA if none found.
#' @import bfast
#' @export
run_bfast <- function(ts_data) {
    tryCatch(
        {
            bp <- bfast::bfastlite(ts_data, breaks = 1, decomp = c("stl", "stlplus"), h = 0.3, stl = "none")
            breakpoint_index <- bp$breakpoints$breakpoints[1]

            if (!is.na(breakpoint_index)) {
                breakdate_year <- bp$data_pp$time[breakpoint_index]
                # Convert decimal year to Date year
                breakdate_year_y <- floor(breakdate_year)
                return(as.integer(breakdate_year_y))
            } else {
                return(NA_integer_)
            }
        },
        error = function(e) {
            return(NA_integer_)
        }
    )
}

#' Run ShapeSelectForest Analysis
#'
#' Applies ShapeSelectForest to classify the trend shape (Flat, Decreasing, Increasing, etc.).
#'
#' @param y_val Numeric vector of the response variable (EVI).
#' @param x_pred Numeric vector of the predictor variable (Years).
#' @param flat Logical. Allow flat shape? Default TRUE.
#' @param dec Logical. Allow decreasing shape? Default TRUE.
#' @param inc Logical. Allow increasing shape? Default TRUE.
#' @return The detected shape code.
#' @import ShapeSelectForest
#' @export
run_ssf_constrained <- function(y_val, x_pred, flat = TRUE, dec = TRUE, inc = TRUE) {
    # SSF requires a matrix for y. It typically expects observations in rows (n x 1).
    y_mat <- matrix(y_val, ncol = 1)

    tryCatch(
        {
            # Suppress internal messages from coneproj/ShapeSelectForest
            utils::capture.output({
                ssf_res <- ShapeSelectForest::shape(x_pred, y_mat, "BIC",
                    get.edf0 = FALSE,
                    flat = flat, dec = dec, inc = inc,
                    jp = FALSE, invee = FALSE, vee = FALSE, db = FALSE,
                    msg = FALSE # Ensure msg is FALSE
                )
            })
            return(ssf_res$shape[1])
        },
        error = function(e) {
            # Keep error logging for now in case of other issues
            message(paste("Error in SSF:", e$message))
            return(NA_integer_)
        }
    )
}
run_ssf_unconstrained <- function(y_val, x_pred, flat = TRUE, dec = TRUE, inc = TRUE) {
    # SSF requires a matrix for y. It typically expects observations in rows (n x 1).
    y_mat <- matrix(y_val, ncol = 1)

    tryCatch(
        {
            # Suppress internal messages from coneproj/ShapeSelectForest
            utils::capture.output({
                ssf_res <- ShapeSelectForest::shape(x_pred, y_mat, "BIC",
                    get.edf0 = FALSE,
                    flat = TRUE, dec = TRUE, inc = TRUE,
                    jp = TRUE, invee = TRUE, vee = TRUE, db = TRUE,
                    msg = FALSE # Ensure msg is FALSE
                )
            })
            return(ssf_res$shape[1])
        },
        error = function(e) {
            # Keep error logging for now in case of other issues
            message(paste("Error in SSF:", e$message))
            return(NA_integer_)
        }
    )
}
#' Run Modified Mann-Kendall Test
#'
#' Applies the Modified Mann-Kendall test with Sen's Slope estimator.
#'
#' @param x Numeric vector of data.
#' @return A list containing Zc, p-value, Tau, and Sen's Slope.
#' @import modifiedmk
#' @export
run_mk <- function(x) {
    tryCatch(
        {
            mk_res <- modifiedmk::mmkh(as.numeric(x), ci = 0.95)
            return(list(
                Zc = mk_res[1],
                pVal = mk_res[2],
                Tau = mk_res[6],
                SenSlope = mk_res[7]
            ))
        },
        error = function(e) {
            return(list(Zc = NA, pVal = NA, Tau = NA, SenSlope = NA))
        }
    )
}

#' Analyze Trends for All Pixels
#'
#' Master function to apply BFAST, SSF, and MK tests to a dataset.
#'
#' @param start_year Optional integer. Start year for analysis.
#' @param end_year Optional integer. End year for analysis.
#' @return An sf object with analysis results and geometry for each pixel.
#' @export
analyze_trends <- function(evi_data, start_year = NULL, end_year = NULL) {
    # Filter by year if provided
    if (!is.null(start_year)) {
        evi_data <- evi_data %>% filter(Year >= start_year)
    }
    if (!is.null(end_year)) {
        evi_data <- evi_data %>% filter(Year <= end_year)
    }

    # Check duration
    years <- unique(evi_data$Year)
    duration <- max(years) - min(years) + 1

    run_ssf_flag <- TRUE
    if (duration < 8) {
        warning(paste("Data duration is", duration, "years. ShapeSelectForest (SSF) requires at least 8 years. SSF will be skipped."))
        run_ssf_flag <- FALSE
    }

    # Extract basename to determine AnalysisPeriod and for export
    basename_attr <- unique(evi_data$Basename)[1]

    analysis_period <- NA_character_
    if (!is.null(basename_attr)) {
        if (grepl("_AM$", basename_attr)) {
            analysis_period <- "AM"
        } else if (grepl("_AWM$", basename_attr)) {
            analysis_period <- "AWM"
        } else if (grepl("_ADM$", basename_attr)) {
            analysis_period <- "ADM"
        }
    }

    results <- list()
    pixels <- unique(evi_data$PixID)

    for (pix in pixels) {
        sub_data <- evi_data %>%
            filter(PixID == pix) %>%
            arrange(Year)
        y_val <- sub_data$EVI

        if ("Month" %in% names(sub_data)) {
            x_val <- sub_data$Year + (sub_data$Month - 1) / 12
        } else {
            x_val <- sub_data$Year
        }

        # Create TS object for BFAST
        pix_start_year <- min(sub_data$Year)

        if ("Month" %in% names(sub_data)) {
            # Monthly data
            # Ensure data is sorted by Year, Month
            sub_data <- sub_data %>% arrange(Year, Month)
            start_month <- sub_data$Month[1]
            ts_data <- ts(sub_data$EVI, start = c(pix_start_year, start_month), frequency = 12)
        } else {
            # Annual data
            ts_data <- ts(sub_data$EVI, start = pix_start_year, frequency = 1)
        }

        # BFAST
        bfast_break <- run_bfast(ts_data)

        # SSF
        if (run_ssf_flag) {
            ssf_shape_code_constrained <- run_ssf_constrained(y_val, x_val)
            ssf_shape_code_unconstrained <- run_ssf_unconstrained(y_val, x_val)
        } else {
            ssf_shape_code_constrained <- NA
            ssf_shape_code_unconstrained <- NA
        }


        # Map SSF code to description
        ssf_labels <- c(
            "1" = "Flat",
            "2" = "Decreasing",
            "3" = "One_Jump",
            "4" = "Inverted_Vee",
            "5" = "Vee",
            "6" = "Linear_Increasing",
            "7" = "Double_Jump"
        )

        if (run_ssf_flag) {
            ssf_shape_constrained <- ssf_labels[as.character(ssf_shape_code_constrained)]
            ssf_shape_unconstrained <- ssf_labels[as.character(ssf_shape_code_unconstrained)]
        } else {
            ssf_shape_constrained <- "Insufficient Data (<8 yrs)"
            ssf_shape_unconstrained <- "Insufficient Data (<8 yrs)"
        }


        # Handle unknown codes
        if (is.na(ssf_shape_constrained)) ssf_shape_constrained <- "Unknown"
        if (is.na(ssf_shape_unconstrained)) ssf_shape_unconstrained <- "Unknown"

        # MK
        mk_res <- run_mk(y_val)

        # Define end_year for the duration string
        pix_end_year <- max(sub_data$Year)

        # Get coordinates (assuming constant per pixel)
        lon <- sub_data$longitude[1]
        lat <- sub_data$latitude[1]

        results[[as.character(pix)]] <- data.frame(
            AnalysisDuration = paste0(pix_start_year, "-", pix_end_year),
            AnalysisPeriod = analysis_period,
            PixID = pix,
            longitude = lon,
            latitude = lat,
            BFAST_Break = bfast_break,
            SSF_Shape_Constrained = ssf_shape_constrained,
            SSF_Shape_Unconstrained = ssf_shape_unconstrained,
            MK_Zc = mk_res$Zc,
            MK_pVal = mk_res$pVal,
            MK_Tau = mk_res$Tau,
            MK_SenSlope = mk_res$SenSlope
        )
    }

    final_df <- bind_rows(results)

    # Convert to sf object
    results_sf <- sf::st_as_sf(final_df, coords = c("longitude", "latitude"), crs = 4326)

    # Export Results (basename_attr is already extracted)
    if (!is.null(basename_attr)) {
        # Define export directory
        export_path <- paste0("Results/", basename_attr)
        dir.create(export_path, recursive = TRUE, showWarnings = FALSE)

        # CSV Export (drop geometry)
        sf::st_write(results_sf,
            dsn = paste0(export_path, "/Trend_Results_", basename_attr, ".csv"),
            layer_options = "GEOMETRY=AS_XY",
            delete_dsn = TRUE, quiet = TRUE
        )

        # GeoJSON Export
        sf::st_write(results_sf,
            dsn = paste0(export_path, "/Trend_Results_", basename_attr, ".geojson"),
            delete_dsn = TRUE, quiet = TRUE
        )

        results_sf$Basename <- basename_attr
    }

    return(results_sf)
}
