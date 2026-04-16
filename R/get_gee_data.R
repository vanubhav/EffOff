utils::globalVariables(c("PixID", "Date", "EVI", "Sensor", "coord_id"))

#' Extract EVI Data from Google Earth Engine
#'
#' This function extracts Enhanced Vegetation Index (EVI) data from Landsat 5, 7, and 8
#' for a given spatial feature (sf object). It harmonizes the data across sensors
#' and returns the raw pixel values.
#'
#' @param sf_object An sf object containing the region of interest (polygons).
#' @param start_year Numeric. The start year for data extraction. Default is 1990.
#' @param end_year Numeric. The end year for data extraction. Default is current year.
#' @param scale Numeric. The scale in meters for the reduction. Default is 30.
#' @param file_name Character. Optional custom name for the output file and task. Default is NULL.
#' @return A dataframe containing the extracted EVI data.
#' @details
#' This function uses `rgee::ee_as_sf(..., via = "drive")` to handle large datasets.
#' **IMPORTANT:**
#' *   You must have authenticated Google Drive with `rgee` (run `ee_Initialize(drive = TRUE)`).
#' *   This package requires access to Google Drive and other Google services via the `googledrive` and `tidyverse` packages for automated file handling.
#' *   An active internet connection is required for this function to communicate with Google Earth Engine and Google Drive.
#'
#' **WARNING:** Extracting pixel-level data for large regions or long time series can generate
#' massive datasets. If the process fails or takes too long, consider splitting your region
#' into smaller polygons. The program will periodically check for the completed file in your
#' Google Drive and download the file when it becomes available.
#' @import rgee
#' @import sf
#' @import dplyr
#' @import googledrive
#' @importFrom dplyr %>%
#' @export
extract_evi <- function(sf_object, start_year = 1990, end_year = as.integer(format(Sys.Date(), "%Y")), scale = 30, file_name = NULL) {
    # Check if rgee is initialized
    if (!rgee::ee_Initialize(drive = TRUE, quiet = TRUE)) {
        stop("rgee not initialized. Please run ee_Initialize(drive = TRUE) first.")
    }

    # Ensure sf object is 2D (remove Z/M dimensions if present)
    sf_object <- sf::st_zm(sf_object, drop = TRUE, what = "ZM")

    # Ensure geometry is valid
    sf_object <- sf::st_make_valid(sf_object)

    # Convert sf object to ee object
    roi <- rgee::sf_as_ee(sf_object)

    # Define time range
    # Start 1 month early (December of previous year) for better gap filling
    prev_year <- start_year - 1
    start_date <- paste0(prev_year, "-12-01")
    end_date <- paste0(end_year, "-12-31")

    # Helper functions

    # Apply scaling factors
    apply_scale_factors <- function(image) {
        opticalBands <- image$select("SR_B.")$multiply(0.0000275)$add(-0.2)
        thermalBands <- image$select("ST_B.*")$multiply(0.00341802)$add(149.0)
        return(image$addBands(opticalBands, NULL, TRUE)$addBands(thermalBands, NULL, TRUE))
    }

    # Calculate EVI
    calc_evi <- function(image) {
        evi <- image$expression(
            "2.5 * ((NIR - Red) / (NIR + 6 * Red - 7.5 * Blue + 1))",
            list(
                NIR = image$select("NIR"),
                Red = image$select("Red"),
                Blue = image$select("Blue")
            )
        )$rename("EVI")
        return(image$addBands(evi))
    }

    # Cloud masking function (Standard Landsat QA)
    mask_clouds <- function(image) {
        qa <- image$select("QA_PIXEL")
        # Bits 3 (cloud) and 4 (cloud shadow)
        cloud_mask <- qa$bitwiseAnd(bitwShiftL(1, 3))$eq(0)
        shadow_mask <- qa$bitwiseAnd(bitwShiftL(1, 4))$eq(0)
        mask <- cloud_mask$And(shadow_mask)
        return(image$updateMask(mask))
    }

    # Rename bands for harmonization
    rename_oli <- function(image) {
        image$select(
            c("SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B6", "SR_B7", "QA_PIXEL"),
            c("Blue", "Green", "Red", "NIR", "SWIR1", "SWIR2", "QA_PIXEL")
        )
    }

    rename_etm <- function(image) {
        image$select(
            c("SR_B1", "SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B7", "QA_PIXEL"),
            c("Blue", "Green", "Red", "NIR", "SWIR1", "SWIR2", "QA_PIXEL")
        )
    }

    # Warning about data volume
    message("WARNING: Extracting pixel-level data for large regions or long time series can generate massive datasets.")
    message("If the process fails or takes too long, consider splitting your region into smaller polygons.")

    tryCatch(
        {
            # Define Filters
            # CLOUD_COVER < 50
            # GEOMETRIC_RMSE_MODEL < 10
            # IMAGE_QUALITY == 9 (or IMAGE_QUALITY_OLI == 9)

            common_filters <- rgee::ee$Filter$And(
                rgee::ee$Filter$bounds(roi),
                rgee::ee$Filter$date(start_date, end_date),
                rgee::ee$Filter$lt("CLOUD_COVER", 50),
                rgee::ee$Filter$lt("GEOMETRIC_RMSE_MODEL", 10)
            )

            # Collections
            l8 <- rgee::ee$ImageCollection("LANDSAT/LC08/C02/T1_L2")$
                filter(common_filters)$
                filter(rgee::ee$Filter$eq("IMAGE_QUALITY_OLI", 9))$
                map(apply_scale_factors)$
                map(rename_oli)$
                map(mask_clouds)$
                map(calc_evi)

            l7 <- rgee::ee$ImageCollection("LANDSAT/LE07/C02/T1_L2")$
                filter(common_filters)$
                filter(rgee::ee$Filter$eq("IMAGE_QUALITY", 9))$
                map(apply_scale_factors)$
                map(rename_etm)$
                map(mask_clouds)$
                map(calc_evi)

            l5 <- rgee::ee$ImageCollection("LANDSAT/LT05/C02/T1_L2")$
                filter(common_filters)$
                filter(rgee::ee$Filter$eq("IMAGE_QUALITY", 9))$
                map(apply_scale_factors)$
                map(rename_etm)$
                map(mask_clouds)$
                map(calc_evi)

            # Merge collections
            merged_col <- l8$merge(l7)$merge(l5)$sort("system:time_start")

            # Function to sample pixels from each image
            sample_pixels <- function(img) {
                date <- img$date()$format("YYYY-MM-dd")
                sensor <- img$get("SENSOR_ID")

                # Add pixel coordinates
                img <- img$addBands(rgee::ee$Image$pixelLonLat())

                # Sample pixels within the ROI
                # We extract EVI, latitude, longitude
                # geometries = TRUE might help avoid some internal rgee/reticulate conversion issues
                # by keeping the object as a standard FeatureCollection with geometry
                samples <- img$select(c("EVI", "latitude", "longitude"))$sample(
                    region = roi,
                    scale = scale,
                    geometries = TRUE,
                    dropNulls = TRUE
                )

                # Add metadata to each feature
                return(samples$map(function(f) {
                    f$set("Date", date)$set("Sensor", sensor)
                }))
            }

            # Flatten the collection of collections into a single FeatureCollection
            fc <- merged_col$map(sample_pixels)$flatten()

            # Define names for Drive task and Local file
            if (!is.null(file_name)) {
                drive_name <- file_name
                local_name <- file_name
            } else {
                timestamp <- format(Sys.time(), "%d%m%Y_%H%M")
                drive_name <- paste0("ROI_EVI_", timestamp)
                local_name <- NULL # Will set to TaskID later
            }

            # Create local export directory if it doesn't exist
            export_dir <- file.path(getwd(), "EffOff_exports")
            if (!dir.exists(export_dir)) {
                dir.create(export_dir, recursive = TRUE)
            }

            message("Starting Google Drive export task...")

            # Create export task
            task <- rgee::ee_table_to_drive(
                collection = fc,
                description = drive_name,
                folder = "EffOff_files",
                fileFormat = "CSV"
            )

            # Start task
            task$start()

            # Get Task ID immediately
            task_id <- task$id
            message(paste("Task started. ID:", task_id))

            # Finalize local filename
            if (is.null(local_name)) {
                local_name <- task_id
            }

            output_file <- file.path(export_dir, paste0(local_name, ".csv"))

            # Inform user about the process
            message("########################################################################")
            message("This task has started and is running on the Google Earth Engine servers.")
            message("Once completed, the file will be transferred to a Google Drive folder, and also downloaded to a local folder in your working directory.")
            message("Please check the Earth Engine Task Manager for progress if the file is not generated in the Google Drive folder. The program will periodically check for the completed file in your Google Drive and download the file when it becomes available.")
            message("Larger regions will take longer times to process. Processing times will also depend on the number of parallel queries to your Google Earth Engine account.")
            message("########################################################################")

            # Download file
            message("Downloading data from Google Drive using googledrive package...")

            # Retry loop to find the file (handling potential Drive indexing latency)
            found_file <- NULL
            attempts <- 0
            max_attempts <- 60 # Wait up to 10 minutes (60 * 10s)

            while (is.null(found_file) && attempts < max_attempts) {
                attempts <- attempts + 1
                message(sprintf("Attempt %d of %d to locate file in 'EffOff_files' folder (waiting 10s)...", attempts, max_attempts))

                # 1. Get the folder ID for "EffOff_files"
                # We assume it exists because GEE created it, or we create/find it
                # But GEE export should have created it if it didn't exist.
                # However, to be safe, we search for the folder first.
                folder_res <- googledrive::drive_find(
                    pattern = "EffOff_files",
                    type = "folder",
                    n_max = 1
                )

                if (nrow(folder_res) == 0) {
                    # If folder not found yet, wait and retry (GEE might be slow creating it)
                    Sys.sleep(10)
                    next
                }

                folder_id <- folder_res$id[1]

                # 2. List files INSIDE that folder
                # Sort by modifiedTime desc to get the latest file first
                folder_files <- googledrive::drive_ls(path = googledrive::as_id(folder_id), orderBy = "modifiedTime desc")

                if (nrow(folder_files) > 0) {
                    # Look for exact match (with or without .csv)
                    target_names <- c(drive_name, paste0(drive_name, ".csv"))
                    exact_matches <- folder_files[folder_files$name %in% target_names, ]

                    if (nrow(exact_matches) > 0) {
                        # Take the first exact match (which is the latest due to sorting)
                        found_file <- exact_matches[1, ]
                    } else {
                        # Fallback: check if any file starts with the filename (split files)
                        partial_matches <- folder_files[grepl(paste0("^", drive_name), folder_files$name), ]
                        if (nrow(partial_matches) > 0) {
                            message("Exact filename match not found, but found files matching pattern. Using the latest one.")
                            found_file <- partial_matches[1, ]
                        }
                    }
                }

                if (is.null(found_file)) {
                    Sys.sleep(10) # Wait before retrying
                }
            }

            if (is.null(found_file)) {
                stop(paste0("Could not find the exported file '", drive_name, "' in 'EffOff_files' folder after ", max_attempts, " attempts. Please check your Drive manually."))
            }

            message(paste("Found file:", found_file$name, "(ID:", found_file$id, ")"))

            # Download
            googledrive::drive_download(
                file = googledrive::as_id(found_file$id),
                path = output_file,
                overwrite = TRUE,
                verbose = TRUE
            )

            message(paste("Data successfully extracted and saved to:", output_file))
            message("The coordinates of pixels are in the EPSG:4326 Geodetic Coordinate System")
            message("Each pixel covers a 900 sq.m. area (30m resolution)")

            # Read the CSV
            evi_df <- read.csv(output_file)


            if (nrow(evi_df) == 0) {
                warning("No EVI data found for the specified region and time period.")
                return(data.frame(PixID = integer(), Date = character(), EVI = numeric(), Sensor = character()))
            }

            # Convert to sf object to assign PixID based on coordinates
            # remove = FALSE keeps the coords in the dataframe
            evi_sf <- sf::st_as_sf(evi_df, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)

            # Assign PixID based on unique coordinates
            coords <- sf::st_coordinates(evi_sf)
            evi_sf$coord_id <- paste(coords[, 1], coords[, 2], sep = "_")

            # Create a mapping from coord_id to PixID
            unique_coords <- unique(evi_sf$coord_id)
            coord_map <- setNames(seq_along(unique_coords), unique_coords)

            evi_sf$PixID <- coord_map[evi_sf$coord_id]

            # Convert back to dataframe and select columns
            evi_final <- evi_sf %>%
                sf::st_drop_geometry() %>%
                dplyr::select("PixID", "Date", "EVI", "Sensor", "longitude", "latitude") %>%
                dplyr::arrange(PixID, Date)

            # Add Basename column
            evi_final$Basename <- local_name

            return(evi_final)
        },
        error = function(e) {
            message("Error during GEE extraction:")
            message(e$message)
            return(NULL)
        }
    )
}
