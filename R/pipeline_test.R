# Load required libraries
library(terra)
library(ggplot2)
library(dplyr)
library(tidyr)

# Set working directory
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

source("R/functions_NDVI.R")
source("R/functions_PSI.R")
source("R/functions_TDiff.R")

# Define input files and parameters
input_file <- "../TreeSpeciesGermany/TreeSpeciesMODIS.tif"
output_folder <- "results/Species_Maps"
NDVI <- "../WZMAllDOYs/Quantiles_209.nc"
start_date <- "2003-01-01"
month_day <- "07-28"
depth <- 100
years <- 2003:2024

# Helper function to create directories if missing
ensure_directory <- function(dir_path) {
  if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE)
}

# Ensure output directories exist
ensure_directory(output_folder)

# Define species file paths
species_list <- list(
  Beech = list(
    psi_nc = "../Allan_Yixuan/PSImean_AMJJA_8days_Bu_bfv_20032024_compressed.nc",
    tdiff_nc = "../Allan_Yixuan/TDiffsum_AMJJA_8days_Bu_bfv20032024_compressed.nc",
    mask = "results/Species_Maps/Beech_mask.tif",
    output_dir = "results/Beech_100"
  ),
  Oak = list(
    psi_nc = "../Allan_Yixuan/PSImean_AMJJA_8days_Ei_bfv_20032024_compressed.nc",
    tdiff_nc = "../Allan_Yixuan/TDiffsum_AMJJA_8days_Ei_bfv20032024_compressed.nc",
    mask = "results/Species_Maps/Oak_mask.tif",
    output_dir = "results/Oak_100"
  ),
  Spruce = list(
    psi_nc = "../Allan_Yixuan/PSImean_AMJJA_8days_Fi_bfv_20032024_compressed.nc",
    tdiff_nc = "../Allan_Yixuan/TDiffsum_AMJJA_8days_Fi_bfv20032024_compressed.nc",
    mask = "results/Species_Maps/Spruce_mask.tif",
    output_dir = "results/Spruce_100"
  ),
  Pine = list(
    psi_nc = "../Allan_Yixuan/PSImean_AMJJA_8days_Ki_bfv_20032024_compressed.nc",
    tdiff_nc = "../Allan_Yixuan/TDiffsum_AMJJA_8days_Ki_bfv20032024_compressed.nc",
    mask = "results/Species_Maps/Pine_mask.tif",
    output_dir = "results/Pine_100"
  )
)

# Function to process a single species
process_species <- function(species_name, paths) {
  cat(paste("Processing", species_name, "...\n"))
  
  ensure_directory(paths$output_dir)
  
  # Convert PSI and TDiff NetCDFs to rasters
  psi_df <- transfer_psi_to_df(paths$psi_nc, start_date)
  save_psi_raster(psi_df, month_day, depth, paths$output_dir)

  tdiff_df <- transfer_tdiff_to_df(paths$tdiff_nc, start_date)
  save_tdiff_raster(tdiff_df, month_day, paths$output_dir)

  # Process NDVI, PSI, TDiff for each year
  species_results <- lapply(years, function(year) {
    process_NDVI_PSI_TDiff_species_year(NDVI, species_name, paths$mask, year, paths$output_dir)
  })
  
  # Combine yearly results
  species_df <- bind_rows(species_results)
  
  # Save data
  save_path <- file.path(paths$output_dir, "Quantiles_PSI_TDiff_df_100.RData")
  save(species_df, file = save_path)
  cat(paste("Saved data for", species_name, "to", save_path, "\n"))
  
  return(species_df)
}

# Process all species
all_results_df <- bind_rows(lapply(names(species_list), function(species) {
  process_species(species, species_list[[species]])
}))

all_results_df <- rbind(Oak, Beech, Spruce, Pine)

# Save combined dataset
output_data_dir <- "results/Data"
ensure_directory(output_data_dir)

final_results_path <- file.path(output_data_dir, "All_Species_Quantiles_PSI_TDiff_100.RData")
save(all_results_df, file = final_results_path)
cat("Final dataset saved.\n")

cat("All data preprocessing completed successfully.\n")
