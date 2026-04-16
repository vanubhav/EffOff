#' Internal helpers for safely extracting names
#' @noRd
.safe_name <- function(x, label = "name", fallback = "Analysis") {
    if (is.null(x)) {
        return(fallback)
    }
    u <- unique(stats::na.omit(x))
    if (!length(u)) {
        return(fallback)
    }
    if (length(u) > 1) warning(sprintf("Multiple values for '%s'; using first.", label), call. = FALSE)
    as.character(u[1])
}

#' Internal Helper: Require columns
#' @noRd
.require_cols <- function(df, cols) {
    missing_cols <- setdiff(cols, names(df))
    if (length(missing_cols) > 0) {
        stop(sprintf("Missing required column(s): %s", paste(missing_cols, collapse = ", ")))
    }
}

#' Internal Helper: Continuous Mapping Limits safely evaluating numeric symmetry exclusively
#' @noRd
.symmetric_limits <- function(x) {
    mx <- max(abs(as.numeric(x)), na.rm = TRUE)
    if (is.infinite(mx) || is.na(mx)) {
        return(c(-1, 1))
    }
    c(-mx, mx)
}

#' Shared Configuration Palette Dictionary
#' @noRd
.effoff_palettes <- function() {
    list(
        SSF_Shape_Constrained = c(
            "Decreasing" = "#FFB54C",
            "Flat" = "#F8D66D",
            "Linear_Increasing" = "#7ABD7E",
            "Double_Jump" = "#D0A1C2",
            "Inverted_Vee" = "#617E96",
            "One_Jump" = "#A08CB2",
            "Vee" = "#A0CDCA",
            "Complex" = "#B0B0B0"
        ),
        SSF_Shape_Unconstrained = c(
            "Decreasing" = "#FFB54C",
            "Flat" = "#F8D66D",
            "Linear_Increasing" = "#7ABD7E",
            "Double_Jump" = "#D0A1C2",
            "Inverted_Vee" = "#617E96",
            "One_Jump" = "#A08CB2",
            "Vee" = "#A0CDCA",
            "Complex" = "#B0B0B0"
        ),
        Agreement_Tau_Sen = c(
            "Both Positive" = "#1a9641",
            "Both Negative" = "#D6604D",
            "Neutral" = "#B0B0B0",
            "Mixed/Neutral" = "#B0B0B0",
            "Aligned" = "#1a9641",
            "Flat_Misaligned" = "#F8D66D",
            "Misaligned" = "#D6604D"
        ),
        Agreement_Triple = c(
            "All Positive" = "#1a9641",
            "All Negative" = "#D6604D",
            "Aligned" = "#1a9641",
            "Flat_Misaligned" = "#F8D66D",
            "Misaligned" = "#D6604D",
            "Mixed/Neutral" = "#B0B0B0"
        ),
        Agreement = c(
            "Aligned" = "#1a9641",
            "Flat_Misaligned" = "#F8D66D",
            "Misaligned" = "#D6604D"
        ),
        MK_Significance = c(
            "Significant" = "#7ABD7E",
            "Not Significant" = "#FFB54C"
        ),
        signed_low = "#D6604D",
        signed_high = "#1a9641"
    )
}

#' Standardized Label Matrix
#' @noRd
.effoff_labels <- function(x) {
    labels <- c(
        "MK_SenSlope" = "Sen's slope (EVI / yr)",
        "MK_Tau" = "Kendall's tau (\U03C4)",
        "MK_pVal" = "Significance (p-value)",
        "MK_pVal_Continuous" = "Significance -log10(p)",
        "SSF_Shape_Constrained" = "Trend Shape (Constrained)",
        "SSF_Shape_Unconstrained" = "Trend Shape (Unconstrained)",
        "BFAST_Break" = "Break Year Extraction",
        "Agreement_Tau_Sen" = "Concordance (MK Array)",
        "Agreement_Triple" = "Concordance (Full Spectrum)",
        "Agreement" = "Concordance Verification",
        "MK_Significance" = "Significant Trends"
    )
    if (x %in% names(labels)) {
        return(labels[[x]])
    }
    return(x)
}


#' Plot Trend Map
#'
#' Visualizes the spatial distribution of trends natively executing explicit continuous vs binary evaluations safely mapping across complex arrays.
#'
#' @param results_sf An sf object containing the analysis results and geometry.
#' @param variable Character. The column name to plot (e.g., "SSF_Shape", "MK_SenSlope").
#' @param title Character. Plot title.
#' @param roi_sf Optional sf object. Region of Interest to overlay.
#' @param legend_title Character. Title for the legend.
#' @param name Character. Array bounding title.
#' @param significance Character argument dictating contour map continuous vs binary evaluation paths exclusively.
#' @return A ggplot object structurally returned natively cleanly bounding layers mapped conditionally against attributes.
#' @import ggplot2
#' @import sf
#' @import ggspatial
#' @import patchwork
#' @export
plot_trend_map <- function(results_sf, variable, title = NULL, roi_sf = NULL, legend_title = NULL, name = NULL, significance = c("binary", "continuous")) {
    signif_mapped <- match.arg(significance)

    if (variable == "all") {
        return(.plot_trend_map_all(results_sf, roi_sf, name, signif_mapped))
    }

    if (variable == "panel4_Diagnostics") {
        return(.build_diagnostics_panel(results_sf))
    }

    if (variable != "Bivariate_Slope_Sig") {
        # Dynamically inject fallback evaluation intercept matrix ensuring variable executes natively
        .require_cols(results_sf, variable)
    }

    if (is.null(name)) name <- .safe_name(results_sf$Basename, fallback = "Analysis")
    if (is.null(title)) title <- .effoff_labels(variable)
    if (is.null(legend_title)) legend_title <- .effoff_labels(variable)

    palettes <- .effoff_palettes()
    geom_type <- unique(as.character(sf::st_geometry_type(results_sf)))
    is_point <- any(grepl("POINT", geom_type))

    if (variable == "Bivariate_Slope_Sig") {
        return(.plot_bivariate(results_sf, is_point, roi_sf))
    }

    p <- ggplot2::ggplot(results_sf)
    if (!is.null(roi_sf)) {
        p <- p + ggplot2::geom_sf(data = roi_sf, fill = "grey95", color = "#571010", linewidth = 0.5)
    }

    # 1. Map variable safely ensuring geometry overlays independently natively
    if (variable == "MK_pVal" && signif_mapped == "continuous") {
        # Inject -log10 translation smoothly bypassing infinite limits
        results_sf$neg_log_p <- -log10(as.numeric(results_sf$MK_pVal) + .Machine$double.eps)
        if (is_point) {
            p <- p + ggplot2::geom_sf(ggplot2::aes(color = neg_log_p), size = 3)
        } else {
            p <- p + ggplot2::geom_sf(ggplot2::aes(fill = neg_log_p), color = NA)
        }
    } else {
        if (is_point) {
            p <- p + ggplot2::geom_sf(ggplot2::aes(color = .data[[variable]]), size = 3)
        } else {
            p <- p + ggplot2::geom_sf(ggplot2::aes(fill = .data[[variable]]), color = NA)
        }
    }

    p <- p + ggplot2::labs(title = title, color = legend_title, fill = legend_title) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "bottom")

    # Guarding Scale and Maps rigorously bypassing internal CRS errors if bounded exclusively
    if (sf::st_is_longlat(results_sf)) {
        warning("Map is in Geographic Coordinates (Long/Lat). Scale bars mapped dynamically.", call. = FALSE)
    }
    p <- p + ggspatial::annotation_scale(location = "br", height = unit(0.1, "cm"), text_cex = 0.5) +
        ggspatial::annotation_north_arrow(location = "tr", which_north = "true", height = unit(0.5, "cm"), width = unit(0.5, "cm"), style = ggspatial::north_arrow_minimal)

    # Apply precise scale styling matrices natively mapping explicit limits
    if (variable %in% c("MK_SenSlope", "MK_Tau")) {
        lims <- .symmetric_limits(results_sf[[variable]])
        if (is_point) {
            p <- p + ggplot2::scale_color_gradient2(low = palettes$signed_low, mid = "grey90", high = palettes$signed_high, midpoint = 0, limits = lims, na.value = "grey50")
        } else {
            p <- p + ggplot2::scale_fill_gradient2(low = palettes$signed_low, mid = "grey90", high = palettes$signed_high, midpoint = 0, limits = lims, na.value = "grey50")
        }
    } else if (variable == "MK_pVal" && signif_mapped == "continuous") {
        if (is_point) {
            p <- p + ggplot2::scale_color_viridis_c(option = "magma", name = "-log10(p)", na.value = "grey50")
        } else {
            p <- p + ggplot2::scale_fill_viridis_c(option = "magma", name = "-log10(p)", na.value = "grey50")
        }
    } else if (variable == "BFAST_Break") {
        if (is_point) {
            p <- p + ggplot2::scale_color_viridis_c(option = "cividis", na.value = "grey50")
        } else {
            p <- p + ggplot2::scale_fill_viridis_c(option = "cividis", na.value = "grey50")
        }
    } else if (variable %in% names(palettes)) {
        cols <- palettes[[variable]]
        if (is_point) p <- p + ggplot2::scale_color_manual(values = cols, na.value = "darkgrey") else p <- p + ggplot2::scale_fill_manual(values = cols, na.value = "darkgrey")
    }

    return(p)
}

#' @noRd
.plot_bivariate <- function(results_sf, is_point, roi_sf = NULL) {
    .require_cols(results_sf, c("MK_SenSlope", "MK_pVal"))

    # Conditional logic gracefully defaulting Biscale dependencies recursively mapping if uninstalled
    if (requireNamespace("biscale", quietly = TRUE)) {
        res_bi <- suppressWarnings(biscale::bi_class(results_sf, x = MK_SenSlope, y = MK_pVal, style = "quantile", dim = 3))
        p <- ggplot2::ggplot(res_bi)
        if (!is.null(roi_sf)) {
            p <- p + ggplot2::geom_sf(data = roi_sf, fill = "grey95", color = "#571010", linewidth = 0.5)
        }
        if (is_point) p <- p + ggplot2::geom_sf(ggplot2::aes(color = bi_class), size = 3) else p <- p + ggplot2::geom_sf(ggplot2::aes(fill = bi_class), color = NA)

        p <- p + biscale::bi_scale_fill(pal = "DkBlue", dim = 3) + biscale::bi_scale_color(pal = "DkBlue", dim = 3) +
            ggplot2::labs(title = "Bivariate: Sen's Slope x Significance") + ggplot2::theme_bw() + ggplot2::theme(legend.position = "none")

        legend <- biscale::bi_legend(pal = "DkBlue", dim = 3, xlab = "Slope Mag", ylab = "Significance", size = 8)

        # Merge smartly natively rendering bounds explicitly across patchwork layouts
        final_p <- p + patchwork::inset_element(legend, left = 0.0, bottom = 0.0, right = 0.3, top = 0.3)
        return(final_p)
    } else {
        # Graceful manual 9-class degradation bypass
        warning("biscale package unavailable. Rendering generic Bivariate proxy.", call. = FALSE)
        return(plot_trend_map(results_sf, "MK_SenSlope", title = "Sen Slope Proxy", roi_sf = roi_sf))
    }
}

#' @noRd
.build_diagnostics_panel <- function(results_sf) {
    .require_cols(results_sf, c("MK_Tau", "MK_SenSlope", "MK_pVal", "SSF_Shape_Constrained", "BFAST_Break"))
    df <- sf::st_drop_geometry(results_sf)
    df$Signif_Cat <- ifelse(as.numeric(df$MK_pVal) < 0.05, "Significant", "Not Significant")

    # Tau vs Sen magnitude evaluation scatter gracefully charting regressions natively
    rho <- suppressWarnings(stats::cor(as.numeric(df$MK_Tau), as.numeric(df$MK_SenSlope), method = "spearman", use = "complete.obs"))
    p1 <- ggplot2::ggplot(df, ggplot2::aes(x = as.numeric(MK_Tau), y = as.numeric(MK_SenSlope), color = Signif_Cat)) +
        ggplot2::geom_point(alpha = 0.6) +
        ggplot2::geom_smooth(method = "lm", color = "black", se = FALSE) +
        ggplot2::scale_color_manual(values = c("Significant" = "#7ABD7E", "Not Significant" = "#FFB54C")) +
        ggplot2::labs(title = "Morphological Alignment: Tau x SenSlope", subtitle = sprintf("Spearman \U03C1: %.2f", rho), x = "Kendall's Tau", y = "Sen's Slope") +
        ggplot2::theme_bw()

    # Density / Violin representation
    p2 <- ggplot2::ggplot(df, ggplot2::aes(x = SSF_Shape_Constrained, y = as.numeric(MK_SenSlope), fill = SSF_Shape_Constrained)) +
        ggplot2::geom_violin(alpha = 0.8) +
        ggplot2::scale_fill_manual(values = .effoff_palettes()$SSF_Shape_Constrained) +
        ggplot2::labs(title = "Distribution: Slope strictly partitioned by SSF Class", x = "", y = "Sen's Slope") +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

    return(p1 + p2 + patchwork::plot_layout(ncol = 2))
}


#' @noRd
.plot_trend_map_all <- function(results_sf, roi_sf, name, signif_mapped) {
    num_pixels <- nrow(results_sf)

    # Calculate native true dimensions locally gracefully bypassing arbitrary 900sqm matrix hardcodes
    projected_crs <- 3857 # Web Mercator default projection
    single_pixel_m2 <- tryCatch(
        {
            med_val <- stats::median(as.numeric(sf::st_area(sf::st_transform(results_sf[1, ], projected_crs))), na.rm = TRUE)
            if (med_val > 0) med_val else 900
        },
        error = function(e) 900
    )

    area_pixels_sqm <- num_pixels * single_pixel_m2

    roi_area_text <- "N/A"
    if (!is.null(roi_sf)) {
        roi_area_val <- tryCatch(as.numeric(sf::st_area(roi_sf)), error = function(e) 0)
        roi_area_text <- paste(format(round(roi_area_val, digits = 0), big.mark = ","), "sq.m.")
    }

    p_basemap <- ggplot2::ggplot() +
        ggplot2::theme_bw() +
        ggplot2::labs(
            title = paste0("Extracted Pixel Locations: ", name),
            subtitle = paste0(
                "Pixels: ", format(num_pixels, big.mark = ","), " | ",
                "Calculated Area: ", format(round(area_pixels_sqm, 0), big.mark = ","), " sq.m. | ",
                "Avg Pixel Res: ", round(single_pixel_m2, 0), "m²\n",
                "RoI Polygon Area: ", roi_area_text
            ),
            x = "Longitude", y = "Latitude"
        )
    if (!is.null(roi_sf)) {
        p_basemap <- p_basemap + ggplot2::geom_sf(data = roi_sf, fill = "grey95", color = "#571010", linewidth = 0.5)
    }
    p_basemap <- p_basemap + ggplot2::geom_sf(data = results_sf, color = "black", alpha = 0.7, size = 1.5)

    # Dynamic Map Generator
    p_bfast <- plot_trend_map(results_sf, "BFAST_Break", roi_sf = roi_sf)
    p_ssf_c <- plot_trend_map(results_sf, "SSF_Shape_Constrained", roi_sf = roi_sf)

    p_slope <- plot_trend_map(results_sf, "MK_SenSlope", roi_sf = roi_sf)
    p_tau <- plot_trend_map(results_sf, "MK_Tau", roi_sf = roi_sf)

    # Conditional pval mapping
    if (signif_mapped == "binary") {
        results_sf$MK_Significance <- ifelse(as.numeric(results_sf$MK_pVal) < 0.05, "Significant", "Not Significant")
        p_pval <- plot_trend_map(results_sf, "MK_Significance", roi_sf = roi_sf)
    } else {
        p_pval <- plot_trend_map(results_sf, "MK_pVal", roi_sf = roi_sf, significance = "continuous")
    }

    if ("Agreement_Triple" %in% names(results_sf)) {
        p_agr_ts <- plot_trend_map(results_sf, "Agreement_Triple", roi_sf = roi_sf)
    } else {
        p_agr_ts <- ggplot2::ggplot() # Fallback for old outputs
    }

    p_diagnostics <- plot_trend_map(results_sf, "panel4_Diagnostics", roi_sf = roi_sf)

    # Build panels securely partitioning objects seamlessly into patchwork matrix outputs
    panel1_structural <- patchwork::wrap_plots(p_basemap, p_bfast, p_ssf_c, ncol = 3) & ggplot2::theme(legend.position = "bottom")
    panel2_statistical <- patchwork::wrap_plots(p_slope, p_tau, p_pval, ncol = 3) & ggplot2::theme(legend.position = "bottom")
    panel3_concordance <- patchwork::wrap_plots(p_diagnostics, p_agr_ts, ncol = 2, widths = c(2, 1)) & ggplot2::theme(legend.position = "bottom")

    return(list(
        Panel1_Structural = panel1_structural,
        Panel2_Statistical = panel2_statistical,
        Panel3_Agreement = panel3_concordance
    ))
}

#' Plot Trend Summary
#'
#' Evaluates explicit parameter counts graphing histogram limits securely via Freedman-Diaconis boundaries dynamically bypassing rigid configurations natively assigning explicit categorical layouts tracking internal pipeline data rigorously.
#'
#' @param results_df Dataframe containing the analysis results.
#' @param category_col Character. The column containing trend categories to dynamically dispatch layout plotting across gracefully!
#' @param name Character mapping pipeline tracking logic explicitly smoothly integrating HTML references uniformly.
#' @return A ggplot distribution graphic array.
#' @import ggplot2
#' @import dplyr
#' @import scales
#' @export
plot_trend_summary <- function(results_df, category_col, name = NULL) {
    if (is.null(name)) name <- .safe_name(results_df$Basename, fallback = "Analysis")

    if (category_col == "all") {
        # Secure generation
        p1 <- plot_trend_summary(results_df, "BFAST_Break")
        p2 <- plot_trend_summary(results_df, "SSF_Shape_Constrained")
        p4 <- plot_trend_summary(results_df, "MK_SenSlope")
        p5 <- plot_trend_summary(results_df, "MK_Tau")

        if ("MK_pVal" %in% names(results_df)) {
            results_df$MK_Significance <- ifelse(as.numeric(results_df$MK_pVal) < 0.05, "Significant", "Not Significant")
            p6 <- plot_trend_summary(results_df, "MK_Significance")
        } else {
            p6 <- ggplot2::ggplot()
        }

        p7 <- if ("Agreement_Triple" %in% names(results_df)) plot_trend_summary(results_df, "Agreement_Triple") else ggplot2::ggplot()

        combined_plot <- patchwork::wrap_plots(p1, p2, p4, p5, p6, p7, ncol = 3, nrow = 2) & ggplot2::theme(legend.position = "bottom")
        return(combined_plot)
    }

    if (is.numeric(results_df[[category_col]])) {
        clean_vals <- stats::na.omit(as.numeric(results_df[[category_col]]))
        if (length(clean_vals) == 0) {
            return(ggplot2::ggplot() +
                ggplot2::labs(title = .effoff_labels(category_col)))
        }

        # Freedman-Diaconis dynamic bin logic gracefully overriding arbitrary defaults flawlessly ensuring visual accuracy
        fd_bins <- tryCatch(length(pretty(clean_vals, n = grDevices::nclass.FD(clean_vals))), error = function(e) 30)

        p <- ggplot2::ggplot(results_df, ggplot2::aes(x = as.numeric(.data[[category_col]]))) +
            ggplot2::geom_histogram(bins = fd_bins, color = "white", fill = "#21908d") +
            ggplot2::labs(y = "Pixel Count", title = paste(.effoff_labels(category_col), "Distribution")) +
            ggplot2::theme_bw()

        return(p)
    } else {
        # Categorical extraction smoothly overriding list limits gracefully bypassing zero-matrix drops
        summary_data <- sf::st_drop_geometry(results_df) %>%
            dplyr::group_by(dplyr::across(dplyr::all_of(category_col))) %>%
            dplyr::summarise(Count = dplyr::n(), .groups = "drop") %>%
            dplyr::mutate(Percentage = Count / sum(Count)) %>%
            dplyr::filter(!is.na(.data[[category_col]]))

        if (nrow(summary_data) == 0) {
            return(ggplot2::ggplot() +
                ggplot2::labs(title = .effoff_labels(category_col)))
        }

        p <- ggplot2::ggplot(summary_data, ggplot2::aes(x = "1", y = Percentage, fill = .data[[category_col]])) +
            ggplot2::geom_bar(stat = "identity", position = "fill", color = "white") +
            ggplot2::scale_y_continuous(labels = scales::percent_format()) +
            ggplot2::labs(x = "", y = "Percentage", title = paste(.effoff_labels(category_col), "Proportion")) +
            ggplot2::theme_bw() +
            ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), legend.position = "bottom")

        palettes <- .effoff_palettes()
        if (category_col %in% names(palettes)) {
            p <- p + ggplot2::scale_fill_manual(values = palettes[[category_col]], na.value = "darkgrey")
        } else {
            p <- p + ggplot2::scale_fill_viridis_d(option = "viridis", na.value = "darkgrey")
        }

        return(p)
    }
}


#' Plot Transition Phase Metrics Output
#'
#' Evaluates explicit matrix shifts parsing phase change dynamically assigning regression matrices to rigorous absolute numeric pathways.
#' @param transition_sf An sf object containing the active Transitions output metric array.
#' @param name Character mapping visual assignment block name.
#' @return A ggplot distribution graphic array explicitly catching NA drops.
#' @import patchwork
#' @export
plot_transition_summary <- function(transition_sf, name = "Transition Block") {
    if (is.null(transition_sf) || nrow(transition_sf) == 0) {
        return(ggplot2::ggplot() +
            ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No transitions detected") +
            ggplot2::theme_void())
    }

    df1 <- sf::st_drop_geometry(transition_sf) %>%
        dplyr::group_by(Trajectory_Shift) %>%
        dplyr::summarise(Count = dplyr::n(), .groups = "drop") %>%
        dplyr::mutate(Percentage = Count / sum(Count))
    if (nrow(df1) == 0) {
        return(ggplot2::ggplot() +
            ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No transitions evaluated") +
            ggplot2::theme_void())
    }

    p1 <- ggplot2::ggplot(df1, ggplot2::aes(x = Trajectory_Shift, y = Count, fill = Trajectory_Shift)) +
        ggplot2::geom_bar(stat = "identity", color = "white") +
        ggplot2::geom_text(ggplot2::aes(label = paste0(Count, " \n(", scales::percent(Percentage, accuracy = 1), ")")), vjust = -0.5, fontface = "bold") +
        ggplot2::scale_fill_manual(values = c("Improvement" = "#1a9641", "Decline" = "#D6604D", "Maintained / Flat" = "#B0B0B0")) +
        ggplot2::labs(title = paste0("Tau Momentum: ", name), x = "", y = "Count") +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "none") +
        ggplot2::expand_limits(y = max(df1$Count) * 1.15)

    df2 <- sf::st_drop_geometry(transition_sf) %>%
        dplyr::filter(!is.na(Shape_Transition) & !grepl("NA", Shape_Transition)) %>%
        dplyr::group_by(Shape_Transition) %>%
        dplyr::summarise(Count = dplyr::n(), .groups = "drop") %>%
        dplyr::mutate(Percentage = Count / sum(Count)) %>%
        dplyr::arrange(dplyr::desc(Count)) %>%
        utils::head(9)

    p2 <- ggplot2::ggplot(df2, ggplot2::aes(x = stats::reorder(Shape_Transition, Count), y = Count)) +
        ggplot2::geom_bar(stat = "identity", color = "white", fill = "#3182bd") +
        ggplot2::coord_flip() +
        ggplot2::geom_text(ggplot2::aes(label = paste0(scales::percent(Percentage, accuracy = 1))), hjust = -0.2, fontface = "bold", size = 3.5) +
        ggplot2::labs(title = "SSF Morphological Reversals", x = "", y = "Count") +
        ggplot2::theme_bw() +
        ggplot2::expand_limits(y = max(df2$Count) * 1.25)

    return((p1) + (p2) + patchwork::plot_layout(widths = c(1, 1)))
}


#' Generate Trend Analysis Report
#'
#' Generates rigorous systematic execution boundaries seamlessly pushing parameter-driven matrices dynamically parsing robust statistics arrays.
#'
#' @section Outputs:
#' The function will explicitly construct `.png` graphic derivations directly interfacing parallel arrays internally wrapping structural, statistical, and mapping bounds directly, producing an index-styled responsive `trend_report.html` master document natively configured against the working bounds perfectly executing.
#'
#' @param results_list A named list mapping phase sequences iteratively gracefully executed.
#' @param results_sf An sf object array evaluating explicitly targeted pixel configurations directly.
#' @param reportType Character dynamically shifting layout evaluation boundaries (full, diagnostics).
#' @param output_file HTML tracking name exclusively bounded.
#' @param output_dir Output mapping directory dynamically hooked natively.
#' @param roi_sf Native subset Region bounds.
#' @param title Document title dynamically evaluated.
#' @param author Creator bounds.
#' @param name Analysis core basename seamlessly integrated.
#' @param ... Execution parameters mapped directly.
#' @export
generate_report <- function(results_list = NULL, results_sf = NULL, reportType = "full", output_file = "trend_report.html", output_dir = getwd(), roi_sf = NULL, title = "Trend Analysis Report", author = "EffOff Package", name = NULL, ...) {
    if (!requireNamespace("rmarkdown", quietly = TRUE)) stop("Package 'rmarkdown' is required.")

    valid_types <- c("full", "maps", "tables", "graphs", "plots", "diagnostics")
    if (!reportType %in% valid_types) stop(paste("Invalid reportType. Choose from:", paste(valid_types, collapse = ", ")))
    if (reportType == "plots") reportType <- "graphs"

    template_path <- system.file("rmd", "trend_report.Rmd", package = "EffOff")
    if (template_path == "") template_path <- "inst/rmd/trend_report.Rmd"

    output_path <- rmarkdown::render(
        input = template_path,
        output_file = output_file,
        output_dir = output_dir,
        params = list(results_list = results_list, results_sf = results_sf, roi_sf = roi_sf, reportType = reportType, title = title, author = author, name = name),
        clean = TRUE, quiet = FALSE, ...
    )
    return(output_path)
}

#' Visualize Results Generic Interface
#' @export
visualize <- function(object, ...) UseMethod("visualize")

#' @export
visualize.sf <- function(object, reportType = "full", output_file = "trend_report.html", output_dir = getwd(), name = NULL, ...) {
    generate_report(results_sf = object, reportType = reportType, output_file = output_file, output_dir = output_dir, name = name, ...)
}

#' @export
visualize.list <- function(object, reportType = "full", output_file = NULL, output_dir = NULL, name = NULL, roi_sf = NULL, ...) {
    if (is.null(name)) name <- if (!is.null(object$Overall) && "Basename" %in% names(object$Overall)) object$Overall$Basename[1] else "Analysis"
    mgmt_duration <- if (!is.null(roi_sf) && "StYP_AV" %in% names(roi_sf)) paste0(roi_sf$StYP_AV[1], "-", roi_sf$EnYP_AV[1]) else if (!is.null(object$Management)) object$Management$AnalysisDuration[1] else "Duration"

    date_str <- format(Sys.Date(), "%Y%m%d")
    if (is.null(output_file)) output_file <- paste0("trend_report_", name, "_", mgmt_duration, "_", date_str, ".html")
    if (is.null(output_dir)) output_dir <- file.path(getwd(), "Results", name)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    return(generate_report(results_list = object, reportType = reportType, output_file = output_file, output_dir = output_dir, roi_sf = roi_sf, name = name, ...))
}
