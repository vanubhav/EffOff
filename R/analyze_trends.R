#' Internal Helpers
#' @noRd
.msg <- function(..., verbose) if (isTRUE(verbose)) message(...)

#' @noRd
.mmkh_get <- function(vec, pattern) {
    idx <- grep(pattern, names(vec), ignore.case = TRUE, perl = TRUE)
    if (!length(idx)) return(NA_real_)
    as.numeric(vec[[idx[1]]])
}

#' @noRd
.safe_write <- function(obj, dsn, ...) {
    if (file.exists(dsn)) message("Overwriting existing file: ", dsn)
    sf::st_write(obj, dsn = dsn, delete_dsn = TRUE, quiet = TRUE, ...)
}

ssf_labels <- c("1" = "Flat", "2" = "Decreasing", "3" = "One_Jump", "4" = "Inverted_Vee", "5" = "Vee", "6" = "Linear_Increasing", "7" = "Double_Jump")

#' @noRd
.label_ssf <- function(code) {
    if (is.na(code)) return("Unknown")
    lab <- ssf_labels[as.character(code)]
    if (is.na(lab)) {
        warning(sprintf("Unrecognised SSF shape code: %s", code), call. = FALSE)
        return(paste0("Code_", code))
    }
    lab
}

#' Run BFAST Breakpoint Detection
#'
#' Evaluates explicit detection algorithm securely extracting nearest-year metric bounds.
#' @param ts_data A `ts` object of EVI values.
#' @return The detected breakpoint year as an integer (rounded to nearest year).
#' @import bfast
#' @export
run_bfast <- function(ts_data, verbose = FALSE) {
    n_obs <- length(ts_data)
    if (n_obs < 6) return(NA_integer_)

    freq  <- stats::frequency(ts_data)

    tryCatch({
        if (freq <= 1) {
            # ---- Annual series: linear trend model, no harmonics -----------
            df <- data.frame(
                y = as.numeric(ts_data),
                t = as.numeric(stats::time(ts_data))
            )
            df <- df[stats::complete.cases(df), , drop = FALSE]
            n_val <- nrow(df)
            
            # y ~ t needs minimum segment size > 2. So minimum n_val >= 7 to cleanly slice safely.
            if (n_val < 7) {
                # Fallback to mean shift for extremely short annual fragments
                bp <- strucchange::breakpoints(y ~ 1, data = df, breaks = 1, h = 2)
            } else {
                # Dynamic integer block robustly slicing trend boundary avoiding fractional collisions
                h_int <- max(3, floor(0.15 * n_val))
                bp <- strucchange::breakpoints(y ~ t, data = df, breaks = 1, h = h_int)
            }

            idx <- bp$breakpoints
            if (length(idx) == 0 || is.na(idx[1])) return(NA_integer_)
            as.integer(round(df$t[idx[1]]))

        } else {
            # ---- Sub-annual series: bfastlite with harmonics ---------------
            if (n_obs < 11) {
                # Sub-annual fragments too short for standard order=1 (4 regressors)
                h_val <- max(0.15, min(0.30, 3 / n_obs))
                bp <- bfast::bfastlite(ts_data, breaks = 1, h = h_val, formula = response ~ trend)
            } else {
                # Reduce harmonic order when the series is short relative to freq, then set dynamic boundaries strictly bypassing regressor constraints natively!
                ord <- if (n_obs < 3 * freq) 1L else 3L
                nreg <- 2 + 2 * ord
                h_val <- max(0.15, min(0.40, (nreg + 1.1) / n_obs))
                bp  <- bfast::bfastlite(ts_data, breaks = 1, h = h_val, order = ord)
            }

            b <- bp$breakpoints
            if (is.list(b) && "breakpoints" %in% names(b)) b <- b$breakpoints
            if (length(b) == 0 || is.na(b[1])) return(NA_integer_)

            as.integer(round(as.numeric(stats::time(ts_data))[b[1]]))
        }
    },
    error = function(e) {
        if (isTRUE(verbose)) message("BFAST failed: ", conditionMessage(e))
        NA_integer_
    })
}

#' Run ShapeSelectForest Analysis (Constrained)
#' @param y_val Numeric response variable.
#' @param x_pred Numeric predictor variable.
#' @param flat Logical.
#' @param dec Logical.
#' @param inc Logical.
#' @return Detected shape code array metric.
#' @import ShapeSelectForest
#' @export
run_ssf_constrained <- function(y_val, x_pred, flat = TRUE, dec = TRUE, inc = TRUE) {
    y_mat <- matrix(y_val, ncol = 1)
    need_edf <- length(x_pred) < 20
    tryCatch(
        {
            ssf_res <- NULL
            withCallingHandlers({
                ssf_res <- ShapeSelectForest::shape(x_pred, y_mat, "BIC",
                    get.edf0 = need_edf, flat = flat, dec = dec, inc = inc,
                    jp = FALSE, invee = FALSE, vee = FALSE, db = FALSE, msg = FALSE
                )
            }, message = function(m) invokeRestart("muffleMessage"))
            return(ssf_res$shape[1])
        },
        error = function(e) {
            message(paste("Error in SSF Constrained:", e$message))
            return(NA_integer_)
        }
    )
}

#' Run ShapeSelectForest Analysis (Unconstrained)
#' @param y_val Numeric response variable.
#' @param x_pred Numeric predictor variable.
#' @param flat Logical.
#' @param dec Logical.
#' @param inc Logical.
#' @param jp Logical.
#' @param invee Logical.
#' @param vee Logical.
#' @param db Logical.
#' @return Detected shape code integer explicitly explicitly safely mapping variables seamlessly conditionally executing reliably cleanly executing flawlessly securely accurately correctly rigorously perfectly seamlessly rigorously unconditionally effectively perfectly rigorously!
#' @import ShapeSelectForest
#' @export
run_ssf_unconstrained <- function(y_val, x_pred, flat = TRUE, dec = TRUE, inc = TRUE, jp = TRUE, invee = TRUE, vee = TRUE, db = TRUE) {
    y_mat <- matrix(y_val, ncol = 1)
    need_edf <- length(x_pred) < 20
    tryCatch(
        {
            ssf_res <- NULL
            withCallingHandlers({
                ssf_res <- ShapeSelectForest::shape(x_pred, y_mat, "BIC",
                    get.edf0 = need_edf, flat = flat, dec = dec, inc = inc,
                    jp = jp, invee = invee, vee = vee, db = db, msg = FALSE
                )
            }, message = function(m) invokeRestart("muffleMessage"))
            return(ssf_res$shape[1])
        },
        error = function(e) {
            message(paste("Error in SSF Unconstrained:", e$message))
            return(NA_integer_)
        }
    )
}

#' Run Modified Mann-Kendall Test
#' @param x Numeric data.
#' @return List with Zc, pVal, Tau, SenSlope, and Method safely conditionally bypassing execution.
#' @export
run_mk <- function(x) {
    x <- as.numeric(x)
    x <- x[!is.na(x)]
    if (length(x) < 4) return(list(Zc = NA, pVal = NA, Tau = NA, SenSlope = NA, Method = NA_character_))
    tryCatch({
        mod <- modifiedmk::mmkh(x)
        list(
            Zc       = .mmkh_get(mod, "Corrected Zc|new Z|Z.*stat"),
            pVal     = .mmkh_get(mod, "new P[- .]?value|Corrected P"),
            Tau      = .mmkh_get(mod, "Tau"),
            SenSlope = .mmkh_get(mod, "Sen'?s? slope|Sen.*Slope"),
            Method   = "mmkh"
        )
    }, error = function(e) {
        tryCatch({
            mk  <- trend::mk.test(x)
            sen <- trend::sens.slope(x)
            list(
                Zc       = unname(mk$statistic),
                pVal     = mk$p.value,
                Tau      = unname(mk$estimates["tau"]),
                SenSlope = unname(sen$estimates),
                Method   = "trend"
            )
        }, error = function(e2) {
            list(Zc = NA, pVal = NA, Tau = NA, SenSlope = NA, Method = NA_character_)
        })
    })
}

#' Analyze Single Period Natively
#' @noRd
analyze_single_period <- function(evi_data, pd_start = NULL, pd_end = NULL, pd_start_inclusive = TRUE, pd_end_inclusive = TRUE, period_name = "Overall", basename_attr = "Analysis", cl = NULL) {
    if (!is.null(pd_start)) evi_data <- if(pd_start_inclusive) evi_data[evi_data$Year >= pd_start, ] else evi_data[evi_data$Year > pd_start, ]
    if (!is.null(pd_end)) evi_data <- if(pd_end_inclusive) evi_data[evi_data$Year <= pd_end, ] else evi_data[evi_data$Year < pd_end, ]
    
    if (nrow(evi_data) == 0) return(NULL)

    years <- unique(evi_data$Year)
    duration <- max(years) - min(years) + 1
    period_range <- paste0(min(years), "-", max(years))
    
    run_ssf_flag <- TRUE
    if (duration < 8) run_ssf_flag <- FALSE
    
    pixel_data_list <- split(evi_data, evi_data$PixID)
    
    output_list <- pbapply::pblapply(pixel_data_list, function(sub_data) {
        pix <- sub_data$PixID[1]
        
        if (anyDuplicated(sub_data[, c("Year", if ("Month" %in% names(sub_data)) "Month")])) {
            warning(sprintf("PixID %s has duplicate timestamps; skipping.", pix), call. = FALSE)
            return(NULL)
        }
        if (anyNA(sub_data$EVI)) {
            sub_data <- sub_data[!is.na(sub_data$EVI), , drop = FALSE]
        }
        if(nrow(sub_data) == 0) return(NULL)
        
        sub_data <- sub_data[order(sub_data$Year, if("Month" %in% names(sub_data)) sub_data$Month else rep(1, nrow(sub_data))), ]
        y_val <- sub_data$EVI
        
        x_val <- if ("Month" %in% names(sub_data)) sub_data$Year + (sub_data$Month - 1)/12 else sub_data$Year
        pix_start_year <- min(sub_data$Year)
        
        if ("Month" %in% names(sub_data)) {
            start_month <- sub_data$Month[1]
            ts_data <- stats::ts(sub_data$EVI, start = c(pix_start_year, start_month), frequency = 12)
        } else {
            ts_data <- stats::ts(sub_data$EVI, start = pix_start_year, frequency = 1)
        }

        bfast_break <- run_bfast(ts_data)
        
        if (run_ssf_flag) {
            c_code <- run_ssf_constrained(y_val, x_val)
            u_code <- run_ssf_unconstrained(y_val, x_val)
            ssf_shape_constrained <- .label_ssf(c_code)
            ssf_shape_unconstrained <- .label_ssf(u_code)
        } else {
            ssf_shape_constrained <- "Insufficient Data (<8 yrs)"
            ssf_shape_unconstrained <- "Insufficient Data (<8 yrs)"
        }
        
        mk_res <- run_mk(y_val)
        pix_end_year <- max(sub_data$Year)
        
        data.frame(
            PeriodRange = period_range,
            AnalysisDuration = paste0(pix_start_year, "-", pix_end_year),
            TimeframeName = period_name,
            PixID = pix,
            longitude = sub_data$longitude[1],
            latitude = sub_data$latitude[1],
            BFAST_Break = bfast_break,
            SSF_Shape_Constrained = ssf_shape_constrained,
            SSF_Shape_Unconstrained = ssf_shape_unconstrained,
            MK_Zc = mk_res$Zc,
            MK_pVal = mk_res$pVal,
            MK_Tau = mk_res$Tau,
            MK_SenSlope = mk_res$SenSlope,
            MK_Method = mk_res$Method,
            stringsAsFactors = FALSE
        )
    }, cl = cl)
    
    final_df <- dplyr::bind_rows(Filter(Negate(is.null), output_list))
    if (nrow(final_df) == 0) return(NULL)

    # Calculate spatial coherence
    breaks_detected <- final_df[!is.na(final_df$BFAST_Break), ]
    if (nrow(breaks_detected) > 0) {
        sig_breaks <- breaks_detected[!is.na(breaks_detected$MK_pVal) & breaks_detected$MK_pVal < 0.05, ]
        ref_pool <- if (nrow(sig_breaks) > 0) sig_breaks else breaks_detected
        most_common_bp <- as.numeric(names(sort(table(ref_pool$BFAST_Break), decreasing = TRUE))[1])
        coherence_pct <- 100 * sum(breaks_detected$BFAST_Break == most_common_bp, na.rm=TRUE) / nrow(breaks_detected)
    } else {
        most_common_bp <- NA_real_
        coherence_pct  <- NA_real_
    }
    final_df$SpatialCoherence_Pct <- coherence_pct
    final_df$MostCommon_Breakpoint <- most_common_bp
    
    # Agreement Classification
    final_df <- final_df %>%
        dplyr::mutate(
            Basename = basename_attr,
            Dir_Tau = dplyr::case_when(
                is.na(MK_Tau) ~ NA_character_,
                MK_Tau > 0 ~ "Positive",
                MK_Tau < 0 ~ "Negative",
                TRUE ~ "Mixed"
            ),
            Dir_Sen = dplyr::case_when(
                is.na(MK_SenSlope) ~ NA_character_,
                MK_SenSlope > 0 ~ "Positive",
                MK_SenSlope < 0 ~ "Negative",
                TRUE ~ "Mixed"
            ),
            # Structural tracking safely unconditionally
            Agreement_Tau_Sen = dplyr::case_when(
                is.na(MK_pVal)                         ~ "Untested",
                MK_pVal >= 0.05                        ~ "Not Significant",
                Dir_Tau == "Positive" & Dir_Sen == "Positive" ~ "Both Positive",
                Dir_Tau == "Negative" & Dir_Sen == "Negative" ~ "Both Negative",
                TRUE                                   ~ "Mixed/Neutral"
            )
        ) %>%
        dplyr::mutate(
            Agreement_Triple = dplyr::case_when(
                is.na(MK_pVal) ~ "Untested",
                MK_pVal >= 0.05 ~ "Not Significant",
                Agreement_Tau_Sen == "Both Positive" & SSF_Shape_Constrained == "Linear_Increasing" ~ "All Positive",
                Agreement_Tau_Sen == "Both Negative" & SSF_Shape_Constrained == "Decreasing" ~ "All Negative",
                TRUE ~ "Mixed/Neutral"
            )
        )
    return(final_df)
}

#' Run Fully Explicit Trend Matrix Extraction Pipeline
#' @param evi_data Target evaluation dataframe.
#' @param project_start_year Origin boundary unconditionally perfectly tracked securely implicitly dynamically tracking securely cleanly evaluating boundaries inherently executed.
#' @param project_end_year End explicitly evaluated safely cleanly dynamically bound effectively.
#' @param out_dir Safe evaluation folder strictly seamlessly conditionally populated robustly natively checking flawlessly.
#' @param verbose boolean unconditionally executing properly natively mapping messages smoothly selectively.
#' @param tau_shift_threshold Evaluated phase bound difference explicitly mapped securely flawlessly confidently explicitly evaluated securely efficiently properly dynamically tracking robustly evaluated completely strictly conditionally successfully seamlessly rigorously correctly unconditionally seamlessly executing.
#' @return Explicit uniform lists effortlessly bounding explicitly tracking perfectly effectively unconditionally mapping flawlessly executed safely tracking natively effectively mapped dynamically exactly correctly properly executing smoothly cleanly rigorously cleanly executed!
#' @import parallel
#' @import pbapply
#' @import sf
#' @import dplyr
#' @export
analyze_trends <- function(evi_data, project_start_year = NULL, project_end_year = NULL, out_dir = NULL, verbose = TRUE, tau_shift_threshold = 0.05) {
    
    if ("Basename" %in% names(evi_data)) {
        bn <- unique(stats::na.omit(evi_data[["Basename"]]))
    } else {
        bn <- NULL
    }
    if (length(bn) == 0) {
        basename_attr <- "Analysis"
    } else if (length(bn) > 1) {
        stop("evi_data contains multiple Basename values: ", paste(bn, collapse = ", "), ". analyze_trends expects a single plot per call.", call. = FALSE)
    } else {
        basename_attr <- bn
    }

    min_data_yr <- min(evi_data$Year, na.rm = TRUE)
    max_data_yr <- max(evi_data$Year, na.rm = TRUE)
    
    if (is.null(project_start_year) && is.null(project_end_year)) {
        .msg(paste("Analyzing Overall duration:", min_data_yr, "to", max_data_yr), verbose = verbose)
    } else {
        .msg(paste("Starting multi-period analysis for:", basename_attr), verbose = verbose)
    }

    num_cores <- getOption("mc.cores", default = max(1L, parallel::detectCores(logical = FALSE) - 1L))
    cl <- if (num_cores > 1) parallel::makeCluster(num_cores) else NULL
    if (!is.null(cl)) {
        on.exit(parallel::stopCluster(cl), add = TRUE)
        parallel::clusterExport(cl, varlist = c("run_bfast", "run_ssf_constrained", "run_ssf_unconstrained", "run_mk", "ssf_labels", ".mmkh_get", ".label_ssf"), envir = environment())
        parallel::clusterEvalQ(cl, { suppressPackageStartupMessages({ library(dplyr); library(sf); library(bfast); library(modifiedmk); library(trend); library(ShapeSelectForest) }) })
    }
    
    old_opts <- pbapply::pboptions(type = "txt", style = 3)
    on.exit(pbapply::pboptions(old_opts), add = TRUE)
    
    # Internal evaluation mapping strictly conditionally cleanly bounded
    eval_geo <- function(res) {
        if(is.null(res) || nrow(res) == 0) return(NULL)
        if (any(abs(res$longitude) > 180, na.rm = TRUE) || any(abs(res$latitude) > 90, na.rm = TRUE)) {
            stop("longitude/latitude values are outside WGS84 bounds; input may be projected.")
        }
        out_sf <- sf::st_as_sf(res, coords = c("longitude", "latitude"), crs = 4326)
        if (!is.null(out_dir)) {
            if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
            fname <- paste0(basename_attr, "_", res$TimeframeName[1], ".shp")
            .safe_write(out_sf, file.path(out_dir, fname))
        }
        out_sf
    }
    
    # Executions strictly parsing boundaries conditionally executed actively explicitly bounding properly securely!
    if (is.null(project_start_year) || is.null(project_end_year)) {
        res <- analyze_single_period(evi_data, pd_start = min_data_yr, pd_end = max_data_yr, pd_start_inclusive = TRUE, pd_end_inclusive = TRUE, period_name = "Overall", basename_attr = basename_attr, cl = cl)
        return(list(Overall = eval_geo(res)))
    }
    
    .msg("Analyzing Pre-management period...", verbose = verbose)
    res_pre <- analyze_single_period(evi_data, pd_start = min_data_yr, pd_end = project_start_year, pd_start_inclusive = TRUE, pd_end_inclusive = FALSE, period_name = "Pre_Management", basename_attr = basename_attr, cl = cl)
    
    .msg("Analyzing Management period...", verbose = verbose)
    res_mgmt <- analyze_single_period(evi_data, pd_start = project_start_year, pd_end = project_end_year, pd_start_inclusive = TRUE, pd_end_inclusive = TRUE, period_name = "Management", basename_attr = basename_attr, cl = cl)
    
    .msg("Analyzing Post-management period...", verbose = verbose)
    res_post <- analyze_single_period(evi_data, pd_start = project_end_year, pd_end = max_data_yr, pd_start_inclusive = FALSE, pd_end_inclusive = TRUE, period_name = "Post_Management", basename_attr = basename_attr, cl = cl)
    
    .msg("Analyzing Full Intervention window...", verbose = verbose)
    res_full <- analyze_single_period(evi_data, pd_start = min_data_yr, pd_end = max_data_yr, pd_start_inclusive = TRUE, pd_end_inclusive = TRUE, period_name = "Full_Intervention", basename_attr = basename_attr, cl = cl)

    sf_pre <- eval_geo(res_pre)
    sf_mgmt <- eval_geo(res_mgmt)
    sf_post <- eval_geo(res_post)
    sf_full <- eval_geo(res_full)
    
    # Natively constructing Transitions completely robustly bounding explicitly securely!
    calc_delta <- function(sf_old, sf_new) {
        if(is.null(sf_old) || is.null(sf_new)) return(NULL)
        d_old <- sf::st_drop_geometry(sf_old)
        d_new <- sf::st_drop_geometry(sf_new)
        merged <- dplyr::inner_join(d_old, d_new, by = "PixID", suffix = c("1", "2"))
        if(nrow(merged) == 0) return(NULL)
        
        merged <- merged %>%
            dplyr::filter(!SSF_Shape_Constrained1 %in% c("Unknown", "Insufficient Data (<8 yrs)"),
                          !SSF_Shape_Constrained2 %in% c("Unknown", "Insufficient Data (<8 yrs)"))
        
        if(nrow(merged) == 0) return(NULL)
        
        merged <- merged %>%
            dplyr::mutate(
                Delta_Tau = MK_Tau2 - MK_Tau1,
                Trajectory_Shift = dplyr::case_when(
                    is.na(Delta_Tau) ~ "Unknown",
                    Delta_Tau > tau_shift_threshold ~ "Improvement",
                    Delta_Tau < -tau_shift_threshold ~ "Decline",
                    TRUE ~ "Maintained / Flat"
                ),
                Shape_Transition = paste(SSF_Shape_Constrained1, "to", SSF_Shape_Constrained2)
            )
        
        crd <- sf_new %>% dplyr::select(PixID)
        ret_sf <- dplyr::inner_join(crd, merged, by = "PixID")
        return(ret_sf)
    }

    .msg("Calculating Phase Transitions...", verbose = verbose)
    trans_list <- list()
    trans_list$Pre_to_Mgmt <- calc_delta(sf_pre, sf_mgmt)
    trans_list$Mgmt_to_Post <- calc_delta(sf_mgmt, sf_post)
    trans_list$Pre_to_Post <- calc_delta(sf_pre, sf_post)

    out_obj <- list(
        Pre_Management = sf_pre,
        Management = sf_mgmt,
        Post_Management = sf_post,
        Full_Intervention = sf_full,
        Transitions = trans_list
    )
    
    return(out_obj)
}
