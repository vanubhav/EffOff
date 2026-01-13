# Test Script for EffOff

# 1. Load Package
library(devtools)
load_all(".")

# 2. Check Package Name
print(paste("Package Name in Desc:", packageDescription("EffOff")$Package))
print(paste("Title:", packageDescription("EffOff")$Title))

# 3. Simulate Data for Report
# Create a dummy sf object mimicking the results
library(sf)
library(dplyr)

# Create a grid
grid <- st_make_grid(st_as_sfc(st_bbox(c(xmin = 0, xmax = 1, ymin = 0, ymax = 1))), n = c(10, 10))
results_sf <- st_sf(ID = 1:100, geometry = grid)

# Add dummy data
set.seed(123)
results_sf$BFAST_Break <- sample(2000:2020, 100, replace = TRUE)
results_sf$SSF_Shape_Constrained <- sample(c("linearIncreasing", "flat", "decreasing"), 100, replace = TRUE)
results_sf$MK_SenSlope <- rnorm(100, 0, 0.1)
results_sf$MK_pVal <- runif(100, 0, 0.1)

# Create dummy ROI
roi_sf <- st_as_sf(st_as_sfc(st_bbox(c(xmin = -0.1, xmax = 1.1, ymin = -0.1, ymax = 1.1))))
roi_sf$State <- "TestState"
roi_sf$Division <- "TestDiv"
roi_sf$Range <- "TestRange"
roi_sf$StYP_AV <- "2000"
roi_sf$EnYP_AV <- "2010"
roi_sf$uid <- "12345"

# 4. Generate Report (Full)
output_file <- "test_report_full.html"
report_path_full <- visualize(results_sf, reportType = "full", output_file = output_file, roi_sf = roi_sf)
print(paste("Full Report generated at:", report_path_full))

# 5. Generate Report (Maps)
output_file_map <- "test_report_maps.html"
report_path_map <- visualize(results_sf, reportType = "maps", output_file = output_file_map, roi_sf = roi_sf)
print(paste("Map Report generated at:", report_path_map))
