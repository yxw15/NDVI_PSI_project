setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

library(ncdf4)
library(terra)
library(dplyr)
library(reshape2)

# --- Config ---
dates_start <- "2003-01-01"
years       <- 2003:2024

rew_files <- list(
  Beech  = "../Allan_Yixuan/RELAWATmean_AMJJA_8days_Bu_bfv20032024_compressed.nc",
  Oak    = "../Allan_Yixuan/RELAWATmean_AMJJA_8days_Ei_bfv20032024_compressed.nc",
  Spruce = "../Allan_Yixuan/RELAWATmean_AMJJA_8days_Fi_bfv20032024_compressed.nc",
  Pine   = "../Allan_Yixuan/RELAWATmean_AMJJA_8days_Ki_bfv20032024_compressed.nc"
)

months_config <- list(
  April  = list(day="04-23"),
  May    = list(day="05-25"),
  June   = list(day="06-26"),
  July   = list(day="07-28"),
  August = list(day="08-29")
)

# Directory template: results_monthly_REW/April/Beech/
out_dir_template <- "results_monthly_REW/{month}/{species}"


# Ensure a directory exists (recursively)
ensure_directory <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive=TRUE)
}

# Read NetCDF REW, return melted data frame (no depth)
rew_nc_to_df <- function(nc_file) {
  nc <- nc_open(nc_file)
  time <- ncvar_get(nc, "time")  # time dimension length
  x    <- ncvar_get(nc, "x")     # x dimension length
  y    <- ncvar_get(nc, "y")     # y dimension length
  rew  <- ncvar_get(nc, "relawat") # [time, y, x]
  dates <- as.Date(dates_start) + time
  nc_close(nc)
  
  # Check actual array shape
  arr_shape <- dim(rew)
  cat("Shape: ", paste(arr_shape, collapse=", "), "\n")
  cat("time length: ", length(dates), "; y length: ", length(y), "; x length: ", length(x), "\n")
  
  # If order is [time, y, x], melt as such:
  df <- reshape2::melt(rew, varnames=c("time_idx", "y_idx", "x_idx"), value.name="relawat")
  df$time <- as.character(dates[df$time_idx])
  df$y <- y[df$y_idx]
  df$x <- x[df$x_idx]
  df <- df[,c("time", "y", "x", "relawat")]
  df <- na.omit(df)
  df
}


# Save rasters for the month-day closest to the target (per year)
save_rew_rasters_for_monthday <- function(df, month, month_day, species, out_dir) {
  ensure_directory(out_dir)
  for (yr in years) {
    # Desired date (e.g., "2003-04-23")
    target_date <- as.Date(sprintf("%d-%s", yr, month_day))
    available_dates <- unique(df$time)
    # Find the date in the NetCDF closest to the target date
    closest_date <- available_dates[which.min(abs(as.Date(available_dates) - target_date))]
    pts <- df %>% filter(time == closest_date)
    if (nrow(pts) == 0) next
    r <- rast(pts[,c("x","y","relawat")], type="xyz")
    crs(r) <- "EPSG:31467"
    file_out <- file.path(out_dir, sprintf("rew_%d.tif", yr))
    writeRaster(r, file_out, overwrite=TRUE)
  }
}

# --- MAIN LOOP ---
for (species in names(rew_files)) {
  cat("Processing:", species, "\n")
  nc_file <- rew_files[[species]]
  df <- rew_nc_to_df(nc_file)
  for (month in names(months_config)) {
    month_day <- months_config[[month]]$day
    out_dir <- gsub("\\{month\\}", month, out_dir_template)
    out_dir <- gsub("\\{species\\}", species, out_dir)
    cat("  Month:", month, "| Day:", month_day, "\n")
    save_rew_rasters_for_monthday(df, month, month_day, species, out_dir)
  }
}

# --- CHECK COMPLETENESS OF RASTER OUTPUTS ---
years <- 2003:2024

all_complete <- TRUE # flag for overall completeness

for (month in names(months_config)) {
  for (species in names(rew_files)) {
    out_dir <- gsub("\\{month\\}", month, out_dir_template)
    out_dir <- gsub("\\{species\\}", species, out_dir)
    
    # Get expected filenames
    expected_files <- file.path(out_dir, sprintf("rew_%d.tif", years))
    # Check which files exist
    files_exist <- file.exists(expected_files)
    
    if (!all(files_exist)) {
      all_complete <- FALSE
      cat("\nMISSING FILES for", species, "in", month, ":\n")
      cat(expected_files[!files_exist], sep="\n")
    }
  }
}

if (all_complete) {
  cat("\n✅ All expected .tif files are present for all months and species!\n")
} else {
  cat("\n⚠️ Some .tif files are missing. See above for details.\n")
}

###
species_maps_files <- list(
  Beech  = "species_map_MODIS/Beech.tif",
  Oak    = "species_map_MODIS/Oak.tif",
  Spruce = "species_map_MODIS/Spruce.tif",
  Pine   = "species_map_MODIS/Pine.tif"
)

months_config <- list(
  July   = list(day = "07-28", NDVI = "../WZMAllDOYs/Quantiles_209.nc"),
  August = list(day = "08-29", NDVI = "../WZMAllDOYs/Quantiles_241.nc")
)

PSI_folder <- "results_monthly"
REW_folder <- "results_monthly_REW"

psi.file.July.Beech <- "July/Beech/psi_2003.tif"
rew.file.July.Beech <- "July/Beech/rew_2003.tif"

library(terra)
library(dplyr)

years <- 2003:2024
species_names <- names(species_maps_files)
months <- names(months_config)

get_file_path <- function(folder, month, species, prefix, year) {
  file.path(folder, month, species, sprintf("%s_%d.tif", tolower(prefix), year))
}

results_list <- list()

for (species in species_names) {
  cat("Processing species:", species, "\n")
  # Load species mask
  species_mask <- rast(species_maps_files[[species]])
  
  for (month in months) {
    cat("  Month:", month, "\n")
    ndvi_path <- months_config[[month]]$NDVI
    ndvi_raster <- rast(ndvi_path)
    
    for (year in years) {
      cat("    Year:", year, "\n")
      # Build file paths
      psi_path <- get_file_path(PSI_folder, month, species, "psi", year)
      rew_path <- get_file_path(REW_folder, month, species, "rew", year)
      
      # Check existence
      if (!file.exists(psi_path) | !file.exists(rew_path)) {
        cat("      Missing PSI or REW file for", species, month, year, "\n")
        next
      }
      
      psi_rast <- rast(psi_path)
      rew_rast <- rast(rew_path)
      
      # Project PSI and REW to NDVI CRS and extent
      psi_proj <- project(psi_rast, ndvi_raster)
      rew_proj <- project(rew_rast, ndvi_raster)
      
      # Stack: NDVI, PSI, REW
      stack_ras <- c(ndvi_raster, psi_proj, rew_proj)
      names(stack_ras) <- c("NDVI", "PSI", "REW")
      
      # Mask by species map
      stack_masked <- mask(stack_ras, species_mask)
      
      df <- as.data.frame(stack_masked, xy=TRUE, na.rm=TRUE)
      if (nrow(df) == 0) next
      
      df$species <- species
      df$month <- month
      df$year <- year
      
      results_list[[length(results_list) + 1]] <- df
    }
  }
}

final_df <- bind_rows(results_list)

# Columns: x, y, NDVI, PSI, REW, species, month, year
head(final_df)

