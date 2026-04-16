#' Summarise Trend Analysis Results
#'
#' Descriptive and domain-specific summaries for a single polygon's trend
#' analysis output at one timeframe. Mirrors the per-plot reporting from the
#' EffOff pipeline: per-pixel distributions, SSF trend proportions
#' (Decreasing / Flat / Linear_Increasing), SSF-vs-Sen alignment, spatial
#' coherence of BFAST breakpoints, and — when the timeframe is "Management" —
#' breakpoint timing relative to the reported project window.
#'
#' Note: summarises ONE sf object (one polygon, one timeframe). For
#' pre → post transitions, use summarise_transition() on a list of
#' timeframes.
#'
#' @param results_sf An `sf` object from the EffOff pipeline. Optional columns:
#'   `MK_SenSlope`, `MK_Tau`, `MK_pVal`, `SSF_Shape_Constrained`,
#'   `SSF_Shape_Unconstrained`, `BFAST_Break`, `Agreement_Tau_Sen`,
#'   `Agreement_Triple`, `Basename`, `TimeframeName`, `StYP_AV`, `EnYP_AV`.
#' @param alpha Significance threshold for MK p-value. Default 0.05.
#' @param pixel_area_m2 Landsat pixel area. Default 900 (30 m).
#' @param digits Rounding for numeric summaries. Default 4.
#'
#' @return An `effoff_summary` list.
#' @import dplyr
#' @importFrom sf st_drop_geometry st_crs
#' @importFrom tibble tibble
#' @importFrom stats median sd IQR quantile
#' @export
summarise_results <- function(results_sf,
                              alpha = 0.05,
                              pixel_area_m2 = 900,
                              digits = 4) {

    stopifnot(inherits(results_sf, "sf"))
    df <- sf::st_drop_geometry(results_sf)

    # ---- meta ---------------------------------------------------------------
    timeframe <- .safe_unique(df$TimeframeName, "TimeframeName", warn=FALSE)
    name      <- .safe_unique(df$Basename,      "Basename", warn=FALSE)

    meta <- tibble::tibble(
        name          = name,
        timeframe     = timeframe,
        n_pixels      = nrow(df),
        area_m2       = nrow(df) * pixel_area_m2,
        pixel_area_m2 = pixel_area_m2,
        crs           = sf::st_crs(results_sf)$input %||% NA_character_
    )

    # ---- continuous (all + significant-only, per Figs 3a/3b) ---------------
    cont_vars <- intersect(c("MK_SenSlope", "MK_Tau", "MK_pVal", "BFAST_Break"), names(df))
    continuous <- dplyr::bind_rows(
        lapply(cont_vars, function(v) .num_stats(df[[v]], v, digits, scope = "all"))
    )
    if ("MK_pVal" %in% names(df) && any(c("MK_SenSlope", "MK_Tau") %in% names(df))) {
        sig_df <- df[!is.na(df$MK_pVal) & df$MK_pVal < alpha, , drop = FALSE]
        sig_vars <- intersect(c("MK_SenSlope", "MK_Tau"), names(sig_df))
        if(nrow(sig_df) > 0) {
            continuous <- dplyr::bind_rows(
                continuous,
                lapply(sig_vars, function(v) .num_stats(sig_df[[v]], v, digits, scope = "significant"))
            )
        }
    }

    # ---- significance -------------------------------------------------------
    significance <- if ("MK_pVal" %in% names(df)) {
        p <- df$MK_pVal
        n_sig <- sum(p < alpha, na.rm = TRUE)
        n_tot <- sum(!is.na(p))
        tibble::tibble(
            alpha           = alpha,
            n_significant   = n_sig,
            n_tested        = n_tot,
            pct_significant = if (n_tot) round(100 * n_sig / n_tot, digits) else NA_real_
        )
    } else tibble::tibble()

    # ---- SSF proportions (3-class, manuscript framing) ---------------------
    ssf_proportions <- .ssf_proportions(df, digits)

    # ---- SSF vs Sen alignment (Section 2.2.4) ------------------------------
    alignment <- .alignment_rates(df, digits)

    # ---- Spatial coherence (Section 2.2.1) ---------------------------------
    coherence <- .spatial_coherence(df, digits)

    # ---- Management-window timing (Figures 5a–c) ---------------------------
    management_timing <- .management_timing(df, timeframe, digits)

    out <- list(
        meta              = meta,
        continuous        = continuous,
        significance      = significance,
        ssf_proportions   = ssf_proportions,
        alignment         = alignment,
        coherence         = coherence,
        management_timing = management_timing
    )
    class(out) <- c("effoff_summary", "list")
    out
}

# ---- internal helpers -------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

.safe_unique <- function(x, label, warn = TRUE) {
    if (is.null(x)) return(NA_character_)
    u <- unique(stats::na.omit(x))
    if (!length(u)) return(NA_character_)
    if (length(u) > 1 && warn)
        warning(sprintf("Multiple values for '%s'; using first.", label), call. = FALSE)
    as.character(u[1])
}

.num_stats <- function(x, name, digits, scope = "all") {
    x_clean <- x[!is.na(x)]
    if (!length(x_clean)) {
        return(tibble::tibble(variable = name, scope = scope, n = 0L, n_na = length(x),
                              mean = NA_real_, median = NA_real_, sd = NA_real_,
                              iqr = NA_real_, min = NA_real_, max = NA_real_,
                              q25 = NA_real_, q75 = NA_real_))
    }
    q <- stats::quantile(x_clean, c(0.25, 0.75), names = FALSE)
    tibble::tibble(
        variable = name, scope = scope,
        n        = length(x_clean),
        n_na     = sum(is.na(x)),
        mean     = round(mean(x_clean), digits),
        median   = round(stats::median(x_clean), digits),
        sd       = round(stats::sd(x_clean), digits),
        iqr      = round(stats::IQR(x_clean), digits),
        min      = round(min(x_clean), digits),
        max      = round(max(x_clean), digits),
        q25      = round(q[1], digits),
        q75      = round(q[2], digits)
    )
}

.ssf_proportions <- function(df, digits) {
    if (!"SSF_Shape_Constrained" %in% names(df)) return(tibble::tibble())
    three_class <- dplyr::case_when(
        df$SSF_Shape_Constrained == "Decreasing"        ~ "Decreasing",
        df$SSF_Shape_Constrained %in% c("Flat", "Complex") ~ "Flat",
        df$SSF_Shape_Constrained == "Linear_Increasing" ~ "Linear_Increasing",
        TRUE                                            ~ NA_character_
    )
    tab <- table(three_class, useNA = "no")
    tot <- sum(tab)
    if (!tot) return(tibble::tibble())
    tibble::tibble(
        class      = names(tab),
        count      = as.integer(tab),
        percentage = round(100 * as.numeric(tab) / tot, digits)
    )
}

.alignment_rates <- function(df, digits, sen_tol = 1e-4) {
    need <- c("SSF_Shape_Constrained", "MK_SenSlope")
    if (!all(need %in% names(df))) return(tibble::tibble())
    sen_dir <- dplyr::case_when(
        is.na(df$MK_SenSlope)        ~ NA_character_,
        df$MK_SenSlope >  sen_tol    ~ "positive",
        df$MK_SenSlope < -sen_tol    ~ "negative",
        TRUE                         ~ "flat"
    )
    cls <- dplyr::case_when(
        is.na(df$SSF_Shape_Constrained) | is.na(sen_dir)                             ~ NA_character_,
        df$SSF_Shape_Constrained == "Linear_Increasing" & sen_dir == "positive"      ~ "aligned",
        df$SSF_Shape_Constrained == "Decreasing"        & sen_dir == "negative"      ~ "aligned",
        df$SSF_Shape_Constrained == "Flat"              & sen_dir == "flat"          ~ "aligned",
        df$SSF_Shape_Constrained == "Flat"              & sen_dir %in% c("positive","negative") ~ "flat_misaligned",
        TRUE                                                                         ~ "misaligned"
    )
    tab <- table(cls, useNA = "no")
    tot <- sum(tab)
    tibble::tibble(
        class      = names(tab),
        count      = as.integer(tab),
        percentage = if (tot) round(100 * as.numeric(tab) / tot, digits) else NA_real_
    )
}

.spatial_coherence <- function(df, digits) {
    if (!"BFAST_Break" %in% names(df)) return(tibble::tibble())
    yrs <- df$BFAST_Break[!is.na(df$BFAST_Break)]
    if (!length(yrs)) return(tibble::tibble())
    modal <- .mode(yrs)
    tibble::tibble(
        n_pixels_with_break = length(yrs),
        modal_year          = modal,
        coherence_pct       = round(100 * sum(yrs == modal) / length(yrs), digits),
        year_range          = paste(range(yrs), collapse = "–")
    )
}

.management_timing <- function(df, timeframe, digits) {
    if (is.na(timeframe) || !grepl("Management", timeframe, ignore.case = TRUE)) return(tibble::tibble())
    need <- c("BFAST_Break", "StYP_AV", "EnYP_AV")
    if (!all(need %in% names(df))) return(tibble::tibble())

    d <- df[!is.na(df$BFAST_Break) & !is.na(df$StYP_AV) & !is.na(df$EnYP_AV), , drop = FALSE]
    if (!nrow(d)) return(tibble::tibble())

    years_from_start <- d$BFAST_Break - d$StYP_AV
    duration         <- pmax(d$EnYP_AV - d$StYP_AV, .Machine$double.eps)
    scaled_pos       <- years_from_start / duration
    within_duration  <- d$BFAST_Break >= d$StYP_AV & d$BFAST_Break <= d$EnYP_AV

    tibble::tibble(
        n_pixels_with_break     = nrow(d),
        mean_years_from_start   = round(mean(years_from_start), digits),
        median_years_from_start = round(stats::median(years_from_start), digits),
        mean_scaled_position    = round(mean(scaled_pos), digits),
        pct_within_duration     = round(100 * mean(within_duration), digits)
    )
}

.mode <- function(x) { ux <- unique(x); ux[which.max(tabulate(match(x, ux)))] }

#' @export
print.effoff_summary <- function(x, ...) {
    cat("<EffOff trend summary>\n")
    m <- x$meta
    cat(sprintf("  %s | %s\n", m$name %||% "NA", m$timeframe %||% "NA"))
    cat(sprintf("  %s pixels (~%s m^2)\n",
                format(m$n_pixels, big.mark = ","),
                format(m$area_m2,  big.mark = ",")))
    if (nrow(x$significance))
        cat(sprintf("  Significant (p < %s): %.1f%% (%s / %s)\n",
                    x$significance$alpha, x$significance$pct_significant,
                    x$significance$n_significant, x$significance$n_tested))
    if (nrow(x$ssf_proportions)) { cat("  SSF proportions:\n");          print(x$ssf_proportions, n = Inf) }
    if (nrow(x$alignment))       { cat("  SSF-vs-Sen alignment:\n");     print(x$alignment,       n = Inf) }
    if (nrow(x$coherence))
        cat(sprintf("  Spatial coherence: %.1f%% on modal year %s\n",
                    x$coherence$coherence_pct, x$coherence$modal_year))
    if (nrow(x$management_timing)) { cat("  Management-window timing:\n"); print(x$management_timing, n = Inf) }
    invisible(x)
}
