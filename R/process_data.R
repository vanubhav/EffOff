#' Clean and Process EVI Data
#'
#' This function cleans the raw EVI data extracted from GEE. It filters invalid values,
#' prioritizes sensors (LS8 > LS5 > LS7) to handle overlaps, and fills missing dates.
#' It also interpolates missing EVI values and carries forward coordinate information.
#' Finally, it calculates seasonal and annual means.
#'
#' @param evi_data A dataframe containing the raw EVI data with columns: PixID, Date, Sensor, EVI, longitude, latitude.
#' @param start_year Numeric. The start year of the analysis. Data before this year (used for gap filling) will be removed after interpolation. Default is 1990.
#' @param wet_months Numeric vector of months considered 'wet'. Default is 6:9.
#' @param dry_months Numeric vector of months considered 'dry'. Default is c(1:5, 10:12).
#' @return A list containing three dataframes: AM (Annual Mean), ADM (Annual Dry Season Mean), and AWM (Annual Wet Season Mean).
#' @import dplyr
#' @import tidyr
#' @import zoo
#' @import lubridate
#' @import stringr
#' @export
clean_evi_data <- function(evi_data, start_year = 1990, wet_months = 6:9, dry_months = c(1:5, 10:12)) {
    # Ensure Date format
    evi_data$Date <- as.Date(evi_data$Date)

    # Filter invalid EVI values
    evi_data <- evi_data %>%
        dplyr::filter(EVI >= -1 & EVI <= 1)

    # Standardize PixID
    # Format PixID to be "Pix_00001" style for consistency
    evi_data <- evi_data %>%
        dplyr::mutate(
            PixID = stringr::str_pad(as.character(PixID), 5, pad = "0"),
            PixID = paste0("Pix_", PixID)
        )

    # Add Year and Month
    evi_data <- evi_data %>%
        dplyr::mutate(
            Year = lubridate::year(Date),
            Month = lubridate::month(Date)
        )

    # Note: We do NOT filter by start_year here yet.
    # We keep the buffer data (e.g., Dec of prev year) for interpolation.

    # Sensor Prioritization: LS8 > LS5 > LS7
    # If multiple sensors have data for the same month/pixel, keep the highest priority one.
    evi_data <- evi_data %>%
        dplyr::mutate(Priority = dplyr::case_when(
            Sensor == "LC08" ~ 3,
            Sensor == "LT05" ~ 2,
            Sensor == "LE07" ~ 1,
            TRUE ~ 0
        )) %>%
        dplyr::group_by(PixID, Year, Month) %>%
        dplyr::arrange(dplyr::desc(Priority)) %>%
        dplyr::slice(1) %>%
        dplyr::ungroup()

    # Fill months
    # Create a complete grid of PixID and Year-Month
    full_grid <- expand.grid(
        PixID = unique(evi_data$PixID),
        Date = seq(
            min(evi_data$Date, na.rm = TRUE),
            max(evi_data$Date, na.rm = TRUE),
            by = "month"
        )
    ) %>%
        dplyr::mutate(
            Year = lubridate::year(Date),
            Month = lubridate::month(Date)
        )

    # Join grid + interpolate EVI
    # Carry forward longitude and latitude
    evi_filled <- full_grid %>%
        dplyr::left_join(
            evi_data %>% dplyr::select(PixID, Year, Month, EVI, Sensor, longitude, latitude),
            by = c("PixID", "Year", "Month")
        ) %>%
        dplyr::group_by(PixID) %>%
        dplyr::arrange(Date) %>%
        dplyr::mutate(
            # interpolate EVI
            EVI = zoo::na.approx(EVI, na.rm = FALSE),
            # Fill missing Sensor with 'Generated_Mean'
            Sensor = tidyr::replace_na(Sensor, "Generated_Mean"),
            # carry coords downward and upward
            longitude = zoo::na.locf(longitude, na.rm = FALSE),
            longitude = zoo::na.locf(longitude, na.rm = FALSE, fromLast = TRUE),
            latitude = zoo::na.locf(latitude, na.rm = FALSE),
            latitude = zoo::na.locf(latitude, na.rm = FALSE, fromLast = TRUE)
        ) %>%
        dplyr::ungroup()

    # Final Filter: Remove data before the start_year
    # This removes the buffer month(s) used for interpolation
    evi_filled <- evi_filled %>%
        dplyr::filter(Year >= start_year)

    # Capture basename attribute
    basename_attr <- attr(evi_data, "basename")

    # Re-attach attribute to filtered data before passing to calculate_means
    if (!is.null(basename_attr)) {
        attr(evi_filled, "basename") <- basename_attr
    }

    # Calculate and return seasonal means
    means_list <- calculate_means(evi_filled, wet_months = wet_months, dry_months = dry_months)

    # Re-attach attribute to each dataframe in the list, appending the season type
    if (!is.null(basename_attr)) {
        for (season in names(means_list)) {
            attr(means_list[[season]], "basename") <- paste0(basename_attr, "_", season)
        }
    }

    return(means_list)
}

#' Calculate Seasonal and Annual Means
#'
#' Calculates Annual Means (AM), Annual Dry Season Means (ADM), and Annual Wet Season Means (AWM).
#'
#' @param cleaned_data A dataframe output from `clean_evi_data`.
#' @param wet_months Numeric vector of months considered 'wet'. Default is 6:9.
#' @param dry_months Numeric vector of months considered 'dry'. Default is c(1:5, 10:12).
#' @return A list containing three dataframes: AM, ADM, AWM.
#' @note At least 2 months are required for each season to calculate a valid mean.
#' @examples
#' \dontrun{
#' # Use custom months: Wet (May-Oct), Dry (Nov-Apr)
#' means <- calculate_means(
#'     cleaned_data,
#'     wet_months = 5:10,
#'     dry_months = c(11, 12, 1, 2, 3, 4)
#' )
#' }
#' @export
calculate_means <- function(cleaned_data, wet_months = 6:9, dry_months = c(1:5, 10:12)) {
    # Validation: Ensure at least 2 months per season
    if (length(wet_months) < 2) {
        stop("At least 2 months are required for the wet season.")
    }
    if (length(dry_months) < 2) {
        stop("At least 2 months are required for the dry season.")
    }

    # Inform user about the months being used
    message("Calculating seasonal means with the following definitions:")
    message(paste("Wet Season Months:", paste(wet_months, collapse = ", ")))
    message(paste("Dry Season Months:", paste(dry_months, collapse = ", ")))

    # Determine name for files
    name <- attr(cleaned_data, "basename")
    if (is.null(name)) name <- "Analysis"

    # Create export directory
    export_path <- paste0("Results/", name)
    dir.create(export_path, recursive = TRUE, showWarnings = FALSE)

    # Helper to export CSV and GeoJSON
    export_dataset <- function(data, prefix, suffix) {
        filename_base <- paste0(prefix, "_", suffix)

        # CSV
        write.csv(data, file.path(export_path, paste0(filename_base, ".csv")), row.names = FALSE)

        # GeoJSON (requires valid coordinates)
        if (nrow(data) > 0) {
            # Convert to sf
            data_sf <- sf::st_as_sf(data, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)

            sf::st_write(data_sf,
                dsn = file.path(export_path, paste0(filename_base, ".geojson")),
                delete_dsn = TRUE, quiet = TRUE
            )
        }
    }

    # Annual Mean
    am <- cleaned_data %>%
        dplyr::group_by(PixID, Year) %>%
        dplyr::summarise(
            EVI = mean(EVI, na.rm = TRUE),
            longitude = dplyr::first(longitude),
            latitude = dplyr::first(latitude),
            .groups = "drop"
        )
    export_dataset(am, "Data_AM", name)

    # Wet Season Mean
    awm <- cleaned_data %>%
        dplyr::filter(Month %in% wet_months) %>%
        dplyr::group_by(PixID, Year) %>%
        dplyr::summarise(
            EVI = mean(EVI, na.rm = TRUE),
            longitude = dplyr::first(longitude),
            latitude = dplyr::first(latitude),
            .groups = "drop"
        )
    export_dataset(awm, "Data_AWM", name)

    # Dry Season Mean
    adm <- cleaned_data %>%
        dplyr::filter(Month %in% dry_months) %>%
        dplyr::group_by(PixID, Year) %>%
        dplyr::summarise(
            EVI = mean(EVI, na.rm = TRUE),
            longitude = dplyr::first(longitude),
            latitude = dplyr::first(latitude),
            .groups = "drop"
        )
    export_dataset(adm, "Data_ADM", name)

    return(list(AM = am, ADM = adm, AWM = awm))
}
