# Load required libraries
library(terra)
library(ggplot2)
library(dplyr)
library(tidyr)

# Set working directory
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

# Source required R functions
source("R/functions_NDVI.R")
source("R/functions_PSI.R")
source("R/functions_TDiff.R")

# Define global parameters
start_date <- "2003-01-01"
depth <- 50
years <- 2003:2024

# Month-specific settings (renamed to avoid duplicate names)
months_config <- list(
  April  = list(month_day = "04-23", NDVI = "../WZMAllDOYs/Quantiles_113.nc"),
  May    = list(month_day = "05-25", NDVI = "../WZMAllDOYs/Quantiles_145.nc"),
  June   = list(month_day = "06-26", NDVI = "../WZMAllDOYs/Quantiles_177.nc"),
  July   = list(month_day = "07-28", NDVI = "../WZMAllDOYs/Quantiles_209.nc"),
  August = list(month_day = "08-29", NDVI = "../WZMAllDOYs/Quantiles_241.nc")
)

# Species-specific settings
species_config <- list(
  Beech = list(
    psi_nc = "../Allan_Yixuan/PSImean_AMJJA_8days_Bu_bfv_20032024_compressed.nc",
    tdiff_nc = "../Allan_Yixuan/TDiffsum_AMJJA_8days_Bu_bfv20032024_compressed.nc",
    mask = "results/Species_Maps/Beech_mask.tif"
  ),
  Oak = list(
    psi_nc = "../Allan_Yixuan/PSImean_AMJJA_8days_Ei_bfv_20032024_compressed.nc",
    tdiff_nc = "../Allan_Yixuan/TDiffsum_AMJJA_8days_Ei_bfv20032024_compressed.nc",
    mask = "results/Species_Maps/Oak_mask.tif"
  ),
  Spruce = list(
    psi_nc = "../Allan_Yixuan/PSImean_AMJJA_8days_Fi_bfv_20032024_compressed.nc",
    tdiff_nc = "../Allan_Yixuan/TDiffsum_AMJJA_8days_Fi_bfv20032024_compressed.nc",
    mask = "results/Species_Maps/Spruce_mask.tif"
  ),
  Pine = list(
    psi_nc = "../Allan_Yixuan/PSImean_AMJJA_8days_Ki_bfv_20032024_compressed.nc",
    tdiff_nc = "../Allan_Yixuan/TDiffsum_AMJJA_8days_Ki_bfv20032024_compressed.nc",
    mask = "results/Species_Maps/Pine_mask.tif"
  )
)

# Utility function to ensure a directory exists
ensure_directory <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  }
}

# Function to process each species for a given month
process_species <- function(species_name, species_paths, NDVI_file, month_day) {
  message(sprintf("Processing %s for %s...", species_name, month_day))
  
  # Ensure output directory exists
  species_paths$output_dir <- file.path(species_paths$output_dir)
  ensure_directory(species_paths$output_dir)
  
  # Convert PSI NetCDF to dataframe and save as raster
  psi_df <- tryCatch({
    transfer_psi_to_df(species_paths$psi_nc, start_date)
  }, error = function(e) {
    stop(sprintf("Error in PSI conversion for %s: %s", species_name, e))
  })
  save_psi_raster(psi_df, month_day, depth, species_paths$output_dir)
  
  # Convert TDiff NetCDF to dataframe and save as raster
  tdiff_df <- tryCatch({
    transfer_tdiff_to_df(species_paths$tdiff_nc, start_date)
  }, error = function(e) {
    stop(sprintf("Error in TDiff conversion for %s: %s", species_name, e))
  })
  save_tdiff_raster(tdiff_df, month_day, species_paths$output_dir)
  
  # Process NDVI, PSI, and TDiff for each year and combine results
  species_results <- lapply(years, function(year) {
    process_NDVI_PSI_TDiff_species_year(NDVI_file, species_name, species_paths$mask, year, species_paths$output_dir)
  })
  species_df <- bind_rows(species_results)
  
  # Save the species-specific combined data
  save_path <- file.path(species_paths$output_dir, paste0("NDVI_PSI_TDiff_df_", gsub("-", "", month_day), ".RData"))
  save(species_df, file = save_path)
  message(sprintf("Saved data for %s to %s", species_name, save_path))
  
  return(species_df)
}

# Function to process data for an entire month across all species
process_month <- function(month_name, months_config, species_config) {
  message(sprintf("\n##### Starting processing for %s #####", month_name))
  
  month_day <- months_config[[month_name]]$month_day
  NDVI_file <- months_config[[month_name]]$NDVI
  
  # Process data for each species within the month
  species_results_list <- lapply(names(species_config), function(species) {
    sp_paths <- species_config[[species]]
    sp_paths$output_dir <- file.path("results_monthly", month_name, species)
    process_species(species, sp_paths, NDVI_file, month_day)
  })
  
  # Combine all species data and add a month column
  month_results_df <- bind_rows(species_results_list)
  month_results_df$month <- month_name
  
  # Save the monthly combined dataset
  output_data_dir <- file.path("results", "Data", month_name)
  ensure_directory(output_data_dir)
  final_results_path <- file.path(output_data_dir, paste0("All_Species_NDVI_PSI_TDiff_", month_name, ".RData"))
  save(month_results_df, file = final_results_path)
  message(sprintf("Final dataset for %s saved to %s", month_name, final_results_path))
  message(sprintf("##### Processing Data for %s completed #####\n", month_name))
  
  return(month_results_df)
}

# Process each month and collect the results
processed_months <- lapply(names(months_config), function(m) {
  process_month(m, months_config, species_config)
})

# Combine data from all months into one data frame
All_Species_Quantiles_PSI_TDiff_year_month <- bind_rows(processed_months)

# Save the combined dataset
final_combined_path <- file.path("results", "Data", "All_Species_Quantiles_PSI_TDiff_year_month.RData")
save(All_Species_Quantiles_PSI_TDiff_year_month, file = final_combined_path)
message("\n🎉 All Data for all months completed successfully! 🎉")

### combine all dataframes ###
# Define species and months
species_list <- c("Beech", "Oak", "Pine", "Spruce")
months_list <- c("April", "May", "June", "July", "August")

# Initialize an empty list to store data
all_species_data <- list()

# Loop through each species and month to load and process the data
for (species in species_list) {
  for (month in months_list) {
    # Construct the directory path
    dir_path <- paste0("results_monthly/", month, "/", species, "/")
    
    # Get all RData files in the directory
    rdata_files <- list.files(path = dir_path, pattern = "NDVI_PSI_TDiff_df_.*\\.RData$", full.names = TRUE)
    
    # Check if at least one file is found
    if (length(rdata_files) > 0) {
      # Load the first available file (assuming there's only one per month per species)
      load(rdata_files[1])
      
      # Ensure the dataframe is named consistently (assuming it's loaded as species_df)
      species_df$species <- species
      species_df$month <- month
      
      # Store in the list
      all_species_data <- append(all_species_data, list(species_df))
      
      message(paste("Loaded:", rdata_files[1]))
    } else {
      message(paste("No file found in:", dir_path))
    }
  }
}

# Combine all dataframes into one
if (length(all_species_data) > 0) {
  final_df <- do.call(rbind, all_species_data)
  
  # Save the combined dataframe
  save(final_df, file = "results/Data/All_species_month_year_Quantiles_PSI_TDiff.RData")
  
  print("All species data successfully combined and saved.")
} else {
  print("No data found for any species and month.")
}

message("\n🎉 All Data for all months completed successfully! 🎉")
