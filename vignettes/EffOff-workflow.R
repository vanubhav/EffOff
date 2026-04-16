## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)


## ----setup--------------------------------------------------------------------
library(EffOff)
library(sf)
library(ggplot2)


## ----eval=FALSE---------------------------------------------------------------
# # Load shapefile
# roi <- st_read("path/to/your/shapefile.shp")
# 
# # Initialize GEE (requires authentication)
# # rgee::ee_Initialize()
# 
# # Extract full-history EVI data
# # This will pull data from 1990 to present by default.
# evi_data <- extract_evi(roi)


## ----eval=FALSE---------------------------------------------------------------
# # Clean and fill gaps
# cleaned_data <- clean_evi_data(evi_data)
# 
# # Calculate Annual and Seasonal Means
# means <- calculate_means(cleaned_data)
# annual_means <- means$AM
# dry_means <- means$ADM
# wet_means <- means$AWM


## ----eval=FALSE---------------------------------------------------------------
# # Analyze trends on Annual Means
# results_list <- analyze_trends(annual_means, project_start_year = 2005, project_end_year = 2012)
# 
# # View the generated periods
# names(results_list)
# # e.g., "Pre_Management", "Management", "Post_Management", "Full_Intervention", "Overall", "Transitions"


## ----eval=FALSE---------------------------------------------------------------
# # Plot trend map for the entire Management period
# mgmt_sf <- results_list$Management
# plot_trend_map(mgmt_sf, variable = "SSF_Shape_Constrained", title = "Vegetation Trend Shapes (Management Phase)")
# 
# # You can analyze indicator agreement maps directly:
# plot_trend_map(mgmt_sf, variable = "Agreement_Triple", title = "Full Agreement Between Tau, Sen, and SSF")
# 
# # You can also generate a robust HTML breakdown for all phases
# # by passing the list to the visualize wrapper:
# visualize(results_list)

