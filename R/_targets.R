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

# Define month-specific settings
month_list <- list(
  April  = list(month_day = "04-23", NDVI = "../WZMAllDOYs/Quantiles_113.nc"),
  May    = list(month_day = "05-25", NDVI = "../WZMAllDOYs/Quantiles_145.nc"),
  June   = list(month_day = "06-26", NDVI = "../WZMAllDOYs/Quantiles_177.nc"),
  July   = list(month_day = "07-28", NDVI = "../WZMAllDOYs/Quantiles_209.nc"),
  August = list(month_day = "08-29", NDVI = "../WZMAllDOYs/Quantiles_241.nc")
)

# Define species-related parameters
species_list <- list(
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

# Helper function to create directories
ensure_directory <- function(dir_path) {
  if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
}

# Function to process each species
process_species <- function(species_name, paths, NDVI, month_day) {
  cat(paste("Processing", species_name, "for", month_day, "...\n"))
  
  ensure_directory(paths$output_dir)  # Ensure output directory exists
  
  # Convert PSI and TDiff NetCDFs to rasters
  psi_df <- tryCatch(transfer_psi_to_df(paths$psi_nc, start_date), 
                     error = function(e) stop("Error in PSI conversion:", e))
  save_psi_raster(psi_df, month_day, depth, paths$output_dir)
  
  tdiff_df <- tryCatch(transfer_tdiff_to_df(paths$tdiff_nc, start_date), 
                       error = function(e) stop("Error in TDiff conversion:", e))
  save_tdiff_raster(tdiff_df, month_day, paths$output_dir)
  
  # Process NDVI, PSI, and TDiff for each year
  species_results <- lapply(years, function(year) {
    process_NDVI_PSI_TDiff_species_year(NDVI, species_name, paths$mask, year, paths$output_dir)
  })
  
  # Save results
  species_df <- bind_rows(species_results)
  save_path <- file.path(paths$output_dir, paste0("NDVI_PSI_TDiff_df_", gsub("-","",month_day), ".RData"))
  save(species_df, file = save_path)
  cat(paste("Saved data for", species_name, "to", save_path, "\n"))
  
  return(species_df)
}

# Process each month
for (month in names(month_list)) {
  cat("\n##### Starting processing for", month, "#####\n")
  
  month_day <- month_list[[month]]$month_day
  NDVI_file <- month_list[[month]]$NDVI
  
  month_species_results <- lapply(names(species_list), function(species) {
    sp_paths <- species_list[[species]]
    sp_paths$output_dir <- file.path("results_monthly", month, species)
    process_species(species, sp_paths, NDVI_file, month_day)
  })
  
  month_results_df <- bind_rows(month_species_results)
  
  # Save combined dataset for the month
  output_data_dir <- file.path("results/Data", month)
  ensure_directory(output_data_dir)
  final_results_path <- file.path(output_data_dir, paste0("All_Species_NDVI_PSI_TDiff_", month, ".RData"))
  save(month_results_df, file = final_results_path)
  cat("Final dataset for", month, "saved.\n")
  
  # Generate plots
  if (nrow(month_results_df) > 0) {
    figures_dir <- file.path("results/Figures", month)
    ensure_directory(figures_dir)
    
    plot_functions <- list(
      plot_series_NDVI_PSI_same_month,
      plot_species_NDVI_PSI,
      plot_correlation_NDVI_PSI_avg,
      plot_correlation_NDVI_PSI_species,
      plot_mean_box_NDVI_PSI,
      plot_mean_box_NDVI_TDiff,
      plot_mean_box_TDiff_PSI,
      plot_box_NDVI_PSI,
      plot_box_TDiff_PSI,
      plot_box_NDVI_TDiff,
      plot_box_merge_NDVI_PSI,
      plot_box_merge_TDiff_PSI,
      plot_box_merge_NDVI_TDiff,
      plot_NDVI_PSIbin_slope,
      plot_TDiff_PSIbin_slope,
      plot_NDVI_TDiff_slope
    )
    
    plot_paths <- paste0(figures_dir, "/", c(
      "series_NDVI_PSI_same_month.png",
      "species_NDVI_PSI.png",
      "correlation_NDVI_PSI_avg.png",
      "correlation_NDVI_PSI_species.png",
      "mean_box_NDVI_PSI.png",
      "mean_box_NDVI_TDiff.png",
      "mean_box_TDiff_PSI.png",
      "box_NDVI_PSI.png",
      "box_TDiff_PSI.png",
      "box_NDVI_TDiff.png",
      "box_merge_NDVI_PSI.png",
      "box_merge_TDiff_PSI.png",
      "box_merge_NDVI_TDiff.png",
      "NDVI_PSIbin_slope.png",
      "TDiff_PSIbin_slope.png",
      "NDVI_TDiff_slope.png"
    ))
    
    mapply(function(f, p) if (exists(deparse(substitute(f)))) f(month_results_df, p), plot_functions, plot_paths)
    
    cat("All plots for", month, "generated successfully.\n")
  } else {
    cat("No data available for plotting for", month, ".\n")
  }
  
  cat("##### Processing for", month, "completed #####\n")
}

cat("\n🎉 All analyses for all months completed successfully! 🎉\n")
