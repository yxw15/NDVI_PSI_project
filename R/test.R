# Load required libraries
library(tidyverse)
library(spdep)
library(sf)
library(ggplot2)

# Set working directory and load the data
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")
load("results/Data/All_Species_Quantiles_PSI_TDiff.RData")

# Subset the data for species "Pine"
pine_data <- subset(all_results_df, species == "Pine")

# Convert your Pine data to an sf object (assuming 'x' is longitude and 'y' is latitude)
pine_data_sf <- st_as_sf(pine_data, coords = c("x", "y"), crs = 4326)

# Transform to a projected CRS (e.g., UTM zone appropriate for your region; here, EPSG:32633 is used as an example)
pine_data_proj <- st_transform(pine_data_sf, crs = 31467)

# Extract the projected coordinates (in meters)
coords_proj <- st_coordinates(pine_data_proj)

# Define the distance threshold as 5000 meters (5 km)
threshold <- 5000

# Create a neighbor list for all points within 5 km
nb <- dnearneigh(coords_proj, 0, threshold)

# Convert the neighbor list into spatial weights (using row-standardized weights)
listw <- nb2listw(nb, style = "W", zero.policy = TRUE)

# Perform Moran's I tests for the three parameters
moran_quantiles <- moran.test(pine_data$Quantiles, listw, zero.policy = TRUE)
moran_soil <- moran.test(pine_data$soil_water_potential, listw, zero.policy = TRUE)
moran_transp <- moran.test(pine_data$transpiration_deficit, listw, zero.policy = TRUE)

# Print the Moran's I test results for each parameter
cat("Moran's I test for Quantiles:\n")
print(moran_quantiles)
cat("\nMoran's I test for Soil Water Potential:\n")
print(moran_soil)
cat("\nMoran's I test for Transpiration Deficit:\n")
print(moran_transp)
