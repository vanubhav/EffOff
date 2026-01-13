#' Plot Trend Map
#'
#' Visualizes the spatial distribution of trends.
#'
#' @param results_sf An sf object containing the analysis results and geometry.
#' @param variable Character. The column name to plot (e.g., "SSF_Shape", "MK_SenSlope").
#' @param title Character. Plot title.
#' @param roi_sf Optional sf object. Region of Interest to overlay.
#' @param legend_title Character. Title for the legend.
#' @return A ggplot object.
#' @import ggplot2
#' @import sf
#' @import ggspatial
#' @import patchwork
#' @export
plot_trend_map <- function(results_sf, variable, title = "Trend Map", roi_sf = NULL, legend_title = NULL, name = NULL) {
    # Determine name: use argument if provided, otherwise check attribute, else default
    if (is.null(name)) {
        name <- attr(results_sf, "basename")
        if (is.null(name)) name <- "Analysis"
    }

    # Default legend title to variable name if not provided
    if (is.null(legend_title)) legend_title <- variable

    # Ensure roi_sf is an sf object if provided
    if (!is.null(roi_sf) && !inherits(roi_sf, "sf")) {
        roi_sf <- sf::st_as_sf(roi_sf)
    }

    # Handle batch plotting
    if (variable == "all") {
        # Create Binary P-value column
        if ("MK_pVal" %in% names(results_sf)) {
            results_sf$MK_Significance <- ifelse(results_sf$MK_pVal < 0.05, "Significant", "Not Significant")
        } else {
            results_sf$MK_Significance <- NA
        }

        # Panel 1:Map, BFAST, SSF Constrained, SSF Unconstrained

        # Calculate stats for subtitle
        num_pixels <- nrow(results_sf)
        # Assuming 30m resolution -> 900 sq.m per pixel
        area_pixels_sqm <- num_pixels * 900


        if (!is.null(roi_sf)) {
            roi_area_val <- sf::st_area(roi_sf)
            # Format nicely
            roi_area_text <- paste(format(round(as.numeric(roi_area_val), digits = 4), big.mark = ","), "sq.m.")
        }

        p_basemap <- ggplot() +
            geom_sf(data = roi_sf, fill = "grey95", color = "#571010", linewidth = 0.5) +
            theme_bw() +
            labs(
                title = paste0("Extracted Pixel Locations for ", name),
                subtitle = paste0(
                    "No. Pixels: ", format(num_pixels, big.mark = ","),
                    " | Pixel Area: ", format(area_pixels_sqm, big.mark = ","), " sq.m.",
                    " | RoI Area: ", roi_area_text
                ),
                x = "Longitude",
                y = "Latitude"
            ) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

        # Add ROI if available
        if (!is.null(roi_sf)) {
            p_basemap <- p_basemap +
                geom_sf(data = roi_sf, fill = "grey95", color = "#571010", linewidth = 0.5)
        }

        # Add Pixel Points
        # results_sf contains the geometry of pixels (POINTS)
        p_basemap <- p_basemap +
            geom_sf(data = results_sf, color = "black", alpha = 0.7, size = 1.5)


        p_bfast <- plot_trend_map(results_sf, "BFAST_Break", title = "BFAST Breakpoint", roi_sf = roi_sf, legend_title = "Year")
        p_ssf_c <- plot_trend_map(results_sf, "SSF_Shape_Constrained", title = "SSF Constrained", roi_sf = roi_sf, legend_title = "Shape")
        p_ssf_u <- plot_trend_map(results_sf, "SSF_Shape_Unconstrained", title = "SSF Unconstrained", roi_sf = roi_sf, legend_title = "Shape")

        panel0 <- p_basemap
        panel1 <- patchwork::wrap_plots(p_bfast, p_ssf_c, p_ssf_u, ncol = 3, nrow = 1)

        # Panel 2: Sen's Slope, Tau, P-value (Binary)
        p_slope <- plot_trend_map(results_sf, "MK_SenSlope", title = "Sen's Slope", roi_sf = roi_sf, legend_title = "Slope")
        p_tau <- plot_trend_map(results_sf, "MK_Tau", title = "Mann-Kendall Tau", roi_sf = roi_sf, legend_title = "Tau")
        p_pval <- plot_trend_map(results_sf, "MK_Significance",
            title = "Significance (p < 0.05)",
            subtitle = "Sen's Slope and MK Tau", roi_sf = roi_sf, legend_title = "Significance"
        )

        panel2 <- patchwork::wrap_plots(p_slope, p_tau, p_pval,
            ncol = 3, nrow = 1
        )

        # Export Plots to PNG
        export_path <- paste0("Results/", name)
        dir.create(export_path, recursive = TRUE, showWarnings = FALSE)

        ggplot2::ggsave(paste0("Basemap_", name, ".png"),
            path = export_path,
            plot = panel0, width = 297, height = 210, units = "mm", bg = "transparent", dpi = 600
        )
        ggplot2::ggsave(paste0("Structural_", name, ".png"),
            path = export_path,
            plot = panel1, width = 297, height = 210, units = "mm", bg = "transparent", dpi = 600
        )
        ggplot2::ggsave(paste0("Statistical_", name, ".png"),
            path = export_path,
            plot = panel2, width = 297, height = 210, units = "mm", bg = "transparent", dpi = 600
        )

        return(list(
            Panel0_Basemap = panel0,
            Panel1_Structural = panel1,
            Panel2_Statistical = panel2
        ))
    }

    # Check geometry type to decide between fill (polygons) and color (points)
    geom_type <- unique(as.character(sf::st_geometry_type(results_sf)))
    is_point <- any(grepl("POINT", geom_type))

    p <- ggplot(results_sf)

    # Add ROI as base layer if provided
    if (!is.null(roi_sf)) {
        p <- p + geom_sf(data = roi_sf, fill = "grey95", color = "#571010", linewidth = 0.5)
    }

    if (is_point) {
        p <- p + geom_sf(aes(color = .data[[variable]]), size = 1.5) +
            labs(title = title, color = legend_title)
    } else {
        p <- p + geom_sf(aes(fill = .data[[variable]]), color = NA) +
            labs(title = title, fill = legend_title)
    }

    p <- p +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        ggspatial::annotation_scale(location = "br", height = unit(0.1, "cm"), text_cex = 0.5) +
        ggspatial::annotation_north_arrow(
            location = "tl", which_north = "true",
            height = unit(0.5, "cm"), width = unit(0.5, "cm"),
            style = ggspatial::north_arrow_fancy_orienteering
        )

    # Apply custom scales based on variable
    if (variable %in% c("MK_SenSlope", "MK_Tau")) {
        if (is_point) {
            p <- p + scale_color_gradient2(low = "brown", mid = "#F8D66D", high = "darkgreen", midpoint = 0)
        } else {
            p <- p + scale_fill_gradient2(low = "brown", mid = "#F8D66D", high = "darkgreen", midpoint = 0)
        }
    } else if (variable %in% c("SSF_Shape_Constrained", "SSF_Shape_Unconstrained")) {
        ssf_colors <- c(
            "Decreasing" = "#FFB54C",
            "Flat" = "#F8D66D",
            "Linear_Increasing" = "#7ABD7E",
            "Double_Jump" = "#D0A1C2",
            "Inverted_Vee" = "#617E96",
            "One_Jump" = "#A08CB2",
            "Vee" = "#A0CDCA"
        )
        if (is_point) {
            p <- p + scale_color_manual(values = ssf_colors)
        } else {
            p <- p + scale_fill_manual(values = ssf_colors)
        }
    } else if (variable == "MK_Significance") {
        sig_colors <- c(
            "Significant" = "darkgreen",
            "Not Significant" = "brown"
        )
        if (is_point) {
            p <- p + scale_color_manual(values = sig_colors, na.value = "darkgrey")
        } else {
            p <- p + scale_fill_manual(values = sig_colors, na.value = "darkgrey")
        }
    } else {
        # Dynamically add scale based on variable type and geometry
        if (is.numeric(results_sf[[variable]])) {
            if (is_point) {
                p <- p + scale_color_viridis_c(option = "viridis")
            } else {
                p <- p + scale_fill_viridis_c(option = "viridis")
            }
        } else {
            if (is_point) {
                p <- p + scale_color_viridis_d(option = "viridis")
            } else {
                p <- p + scale_fill_viridis_d(option = "viridis")
            }
        }
    }


    return(p)
}

#' Plot Trend Summary
#'
#' Creates a stacked bar chart of trend categories.
#'
#' @param results_df Dataframe containing the analysis results.
#' @param category_col Character. The column containing trend categories.
#' @return A ggplot object.
#' @import ggplot2
#' @import dplyr
#' @import scales
#' @export
plot_trend_summary <- function(results_df, category_col, name = NULL) {
    # Determine name
    if (is.null(name)) {
        name <- attr(results_df, "basename")
        if (is.null(name)) name <- "Analysis"
    }

    # Handle batch plotting
    if (category_col == "all") {
        # Create Binary P-value column if needed
        if ("MK_pVal" %in% names(results_df)) {
            results_df$MK_Significance <- ifelse(results_df$MK_pVal < 0.05, "Significant", "Not Significant")
        } else {
            results_df$MK_Significance <- NA
        }

        # Generate individual plots
        p1 <- plot_trend_summary(results_df, "BFAST_Break") + labs(title = "BFAST Breakpoint")
        p2 <- plot_trend_summary(results_df, "SSF_Shape_Constrained") + labs(title = "SSF 'Constrained'")
        p3 <- plot_trend_summary(results_df, "SSF_Shape_Unconstrained") + labs(title = "SSF 'Unconstrained'")
        p4 <- plot_trend_summary(results_df, "MK_SenSlope") + labs(title = "Sen's Slope")
        p5 <- plot_trend_summary(results_df, "MK_Tau") + labs(title = "Mann-Kendall Tau")
        p6 <- plot_trend_summary(results_df, "MK_Significance") + labs(title = "Significance (p < 0.05)")

        # Combine into a grid
        combined_plot <- patchwork::wrap_plots(p1, p2, p3, p4, p5, p6, ncol = 3, nrow = 2)

        # Export Plot to PNG
        export_path <- paste0("Results/", name)
        dir.create(export_path, recursive = TRUE, showWarnings = FALSE)

        ggplot2::ggsave(paste0("Summary_", name, ".png"),
            path = export_path,
            plot = combined_plot, width = 297, height = 420, units = "mm", bg = "transparent", dpi = 600
        )

        return(combined_plot)
    }

    # Check if variable is numeric
    if (is.numeric(results_df[[category_col]])) {
        # For numeric variables, plot a histogram
        if (category_col %in% c("MK_SenSlope", "MK_Tau")) {
            ggplot(results_df, aes(x = .data[[category_col]], fill = after_stat(x))) +
                geom_histogram(bins = 30, color = "white") +
                scale_fill_gradient2(low = "brown", mid = "#F8D66D", high = "darkgreen") +
                labs(y = "Count", title = "Distribution") +
                theme_bw()
        } else {
            ggplot(results_df, aes(x = .data[[category_col]])) +
                geom_histogram(bins = 30, fill = "#21908d", color = "white") + # Viridis teal
                labs(y = "Count", title = "Distribution") +
                theme_bw()
        }
    } else {
        # For categorical variables, plot a stacked bar chart
        summary_data <- results_df %>%
            group_by(across(all_of(category_col))) %>%
            summarise(Count = n()) %>%
            mutate(Percentage = Count / sum(Count))

        p <- ggplot(summary_data, aes(x = "1", y = Percentage, fill = .data[[category_col]])) +
            geom_bar(stat = "identity", position = "fill") +
            scale_y_continuous(labels = scales::percent_format()) +
            labs(x = "", y = "Percentage", title = "Trend Distribution") +
            theme_bw() +
            theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

        if (category_col %in% c("SSF_Shape_Constrained", "SSF_Shape_Unconstrained")) {
            ssf_colors <- c(
                "Decreasing" = "#FFB54C",
                "Flat" = "#F8D66D",
                "Linear_Increasing" = "#7ABD7E",
                "Double_Jump" = "#D0A1C2",
                "Inverted_Vee" = "#617E96",
                "One_Jump" = "#A08CB2",
                "Vee" = "#A0CDCA"
            )
            p <- p + scale_fill_manual(values = ssf_colors)
        } else if (category_col == "MK_Significance") {
            sig_colors <- c(
                "Significant" = "darkgreen",
                "Not Significant" = "brown"
            )
            p <- p + scale_fill_manual(values = sig_colors, na.value = "darkgrey")
        } else {
            p <- p + scale_fill_viridis_d(option = "viridis")
        }

        return(p)
    }
}

#' Generate Trend Analysis Report
#'
#' Generates an HTML report summarizing the vegetation trend analysis results.
#'
#' @param results_sf An sf object containing the analysis results.
#' @param reportType Character. The type of report to generate. Options: "full" (default), "maps", "tables", "graphs".
#' @param output_file Character. The name of the output HTML file.
#' @param output_dir Character. The directory to save the report. Default is current working directory.
#' @param roi_sf Optional sf object. Region of Interest to overlay.
#' @param title Character. Report title. Default is "Trend Analysis Report".
#' @param author Character. Report author. Default is "EffOff Package".
#' @param ... Additional arguments passed to rmarkdown::render.
#' @return The path to the generated HTML file.
#' @export
generate_report <- function(results_sf,
                            reportType = "full",
                            output_file = "trend_report.html",
                            output_dir = getwd(),
                            roi_sf = NULL,
                            title = "Trend Analysis Report",
                            author = "EffOff Package",
                            name = NULL,
                            ...) {
    if (!requireNamespace("rmarkdown", quietly = TRUE)) {
        stop("Package 'rmarkdown' is required to generate reports. Please install it.")
    }

    # Validate reportType
    valid_types <- c("full", "maps", "tables", "graphs", "plots")
    if (!reportType %in% valid_types) {
        stop(paste("Invalid reportType. Choose from:", paste(valid_types, collapse = ", ")))
    }
    if (reportType == "plots") reportType <- "graphs"

    # Find template
    template_path <- system.file("rmd", "trend_report.Rmd", package = "EffOff")
    if (template_path == "") {
        # If installed package doesn't have it (yet), check local inst/rmd for dev
        template_path <- file.path("inst", "rmd", "trend_report.Rmd")
        if (!file.exists(template_path)) {
            # Fallback to check if we are in the package root
            template_path <- file.path("trend_report.Rmd")
            if (!file.exists(template_path)) {
                stop("Report template 'trend_report.Rmd' not found.")
            }
        }
    }

    # Render report
    output_path <- rmarkdown::render(
        input = template_path,
        output_file = output_file,
        output_dir = output_dir,
        params = list(
            roi_sf = roi_sf,
            reportType = reportType,
            title = title,
            author = author,
            name = name
        ),
        clean = TRUE,
        quiet = FALSE,
        ...
    )

    return(output_path)
}

#' Visualize Results
#'
#' Generic function to visualize or report on analysis results.
#'
#' @param object An object to visualize (typically an sf object with results).
#' @param ... Additional arguments passed to specific methods.
#' @export
visualize <- function(object, ...) {
    UseMethod("visualize")
}

#' @describeIn visualize Generate a report for sf objects.
#' @param reportType Character. Type of report ("full", "maps", "tables", "graphs").
#' @param output_file Output filename.
#' @param output_dir Output directory.
#' @export
visualize.sf <- function(object, reportType = "full", output_file = "trend_report.html", output_dir = getwd(), name = NULL, ...) {
    generate_report(results_sf = object, reportType = reportType, output_file = output_file, output_dir = output_dir, name = name, ...)
}
