# 🎉 Starting NDVI-PSI-TDiff rootzone processing script...

# 📚 Load required libraries
library(terra)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ncdf4)

# 📂 Set working directory
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

# 🔧 Global parameters
dates_start <- "2003-01-01"
years <- 2003:2024

# 🗓️ Month config
months_config <- list(
  April  = list(day = "04-23", NDVI = "../WZMAllDOYs/Quantiles_113.nc"),
  May    = list(day = "05-25", NDVI = "../WZMAllDOYs/Quantiles_145.nc"),
  June   = list(day = "06-26", NDVI = "../WZMAllDOYs/Quantiles_177.nc"),
  July   = list(day = "07-28", NDVI = "../WZMAllDOYs/Quantiles_209.nc"),
  August = list(day = "08-29", NDVI = "../WZMAllDOYs/Quantiles_241.nc")
)

# 🌳 Species config
species_config <- list(
  Beech  = list(psi = "../ALLAN_PIA_SoilMoisture/PSIM/PSImean_rootzone_AMJJA_8days_Bu_bfv_20032024_compressed.nc",
               tdiff = "../ALLAN_PIA_SoilMoisture/TDIFF/TDiffsum_AMJJA_8days_Bu_bfv20032024_compressed.nc",
                mask = "species_map_MODIS/Beech.tif"),
  Oak    = list(psi = "../ALLAN_PIA_SoilMoisture/PSIM/PSImean_rootzone_AMJJA_8days_Ei_bfv_20032024_compressed.nc",
                tdiff = "../ALLAN_PIA_SoilMoisture/TDIFF/TDiffsum_AMJJA_8days_Ei_bfv20032024_compressed.nc",
                mask = "species_map_MODIS/Oak.tif"),
  Spruce = list(psi = "../ALLAN_PIA_SoilMoisture/PSIM/PSImean_rootzone_AMJJA_8days_Fi_bfv_20032024_compressed.nc",
                tdiff = "../ALLAN_PIA_SoilMoisture/TDIFF/TDiffsum_AMJJA_8days_Fi_bfv20032024_compressed.nc",
                mask = "species_map_MODIS/Spruce.tif"),
  Pine   = list(psi = "../ALLAN_PIA_SoilMoisture/PSIM/PSImean_rootzone_AMJJA_8days_Ki_bfv_20032024_compressed.nc",
                tdiff = "../ALLAN_PIA_SoilMoisture/TDIFF/TDiffsum_AMJJA_8days_Ki_bfv20032024_compressed.nc",
                mask = "species_map_MODIS/Pine.tif")
)

# Utility: ensure directory exists
ensure_directory <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

# -----------------------------------------------------------------------------
# Function: transfer_psi_to_df
# -----------------------------------------------------------------------------
transfer_psi_to_df <- function(nc_file) {
  message("🌊 Loading PSI NetCDF: ", nc_file)
  nc <- nc_open(nc_file)
  time_vals <- ncvar_get(nc, "time")
  x_vals <- ncvar_get(nc, "x")
  y_vals <- ncvar_get(nc, "y")
  psi_array <- ncvar_get(nc, "psi")
  nc_close(nc)
  
  dates <- as.Date(dates_start) + time_vals
  dimnames(psi_array) <- list(x = x_vals, y = y_vals, time = as.character(dates))
  
  df <- as.data.frame.table(psi_array, responseName = "soil_water_potential", stringsAsFactors = FALSE) %>%
    mutate(
      x = as.numeric(as.character(x)),
      y = as.numeric(as.character(y)),
      time = as.Date(time)
    ) %>%
    drop_na(soil_water_potential)
  
  message("✅ PSI df ready: ", nrow(df), " rows.")
  return(df)
}

# -----------------------------------------------------------------------------
# Function: save_psi_raster
# -----------------------------------------------------------------------------
save_psi_raster <- function(df, month_day, out_dir) {
  message("🌊 Saving PSI rasters for ", month_day)
  ensure_directory(out_dir)
  sub <- df %>% filter(format(time, "%m-%d") == month_day)
  
  for (yr in years) {
    yr_pts <- sub %>% filter(format(time, "%Y") == yr)
    if (nrow(yr_pts) == 0) next
    r <- rast(yr_pts[, c("x", "y", "soil_water_potential")], type = "xyz")
    crs(r) <- "epsg:31467"
    path <- file.path(out_dir, sprintf("psi_%d.tif", yr))
    writeRaster(r, path, overwrite = TRUE)
    message("   🔹 PSI year ", yr, " saved.")
  }
  message("✅ PSI rasters done for ", month_day)
}

# -----------------------------------------------------------------------------
# Function: transfer_tdiff_to_df
# -----------------------------------------------------------------------------
transfer_tdiff_to_df <- function(nc_file) {
  message("🔥 Loading TDiff NetCDF: ", nc_file)
  nc <- nc_open(nc_file)
  time_vals <- ncvar_get(nc, "time")
  x_vals <- ncvar_get(nc, "x")
  y_vals <- ncvar_get(nc, "y")
  tdiff_array <- ncvar_get(nc, "tdiff")
  nc_close(nc)
  
  dates <- as.Date(dates_start) + time_vals
  dimnames(tdiff_array) <- list(x = x_vals, y = y_vals, time = as.character(dates))
  
  df <- as.data.frame.table(tdiff_array, responseName = "transpiration_deficit", stringsAsFactors = FALSE) %>%
    mutate(
      x = as.numeric(as.character(x)),
      y = as.numeric(as.character(y)),
      time = as.Date(time)
    ) %>%
    drop_na(transpiration_deficit)
  
  message("✅ TDiff df ready: ", nrow(df), " rows.")
  return(df)
}

# -----------------------------------------------------------------------------
# Function: save_tdiff_raster
# -----------------------------------------------------------------------------
save_tdiff_raster <- function(df, month_day, out_dir) {
  message("🔥 Saving TDiff rasters for ", month_day)
  ensure_directory(out_dir)
  sub <- df %>% filter(format(time, "%m-%d") == month_day)
  
  for (yr in years) {
    yr_pts <- sub %>% filter(format(time, "%Y") == yr)
    if (nrow(yr_pts) == 0) next
    r <- rast(yr_pts[, c("x", "y", "transpiration_deficit")], type = "xyz")
    crs(r) <- "epsg:31467"
    path <- file.path(out_dir, sprintf("tdiff_%d.tif", yr))
    writeRaster(r, path, overwrite = TRUE)
    message("   🔹 TDiff year ", yr, " saved.")
  }
  message("✅ TDiff rasters done for ", month_day)
}

# -----------------------------------------------------------------------------
# Function: process_NDVI_PSI_TDiff_species_year
# -----------------------------------------------------------------------------
process_NDVI_PSI_TDiff_species_year <- function(NDVI_file, species, year, out_dir) {
  message("🌿 [", species, "] Combining year ", year)
  
  ### 1️⃣ read the multi‐layer NDVI file, pick the layer for this year
  ndvi_all  <- rast(NDVI_file)
  layer_idx <- year - min(years) + 1
  ndvi_year <- ndvi_all[[layer_idx]]
  
  ### 2️⃣ mask it to the species extent
  mask_r  <- rast(species_config[[species]]$mask)
  ndvi_m  <- mask(ndvi_year, mask_r)
  
  ### 3️⃣ **rename** that one band to "Quantiles" ###
  names(ndvi_m) <- "Quantiles"
  
  ### 4️⃣ load your PSI & TDiff rasters for this year, reproject & mask
  psi_r   <- project(rast(file.path(out_dir, sprintf("psi_%d.tif", year))), mask_r, method = "bilinear")
  tdiff_r <- project(rast(file.path(out_dir, sprintf("tdiff_%d.tif", year))), mask_r, method = "bilinear")
  
  ### 5️⃣ save the masked rasters (optional)
  writeRaster(ndvi_m,  file.path(out_dir, sprintf("NDVI_%d_mask.tif",    year)), overwrite = TRUE)
  writeRaster(psi_r,   file.path(out_dir, sprintf("PSI_%d_mask.tif",     year)), overwrite = TRUE)
  writeRaster(tdiff_r, file.path(out_dir, sprintf("TDiff_%d_mask.tif",   year)), overwrite = TRUE)
  
  ### 6️⃣ stack the three single‐band rasters and convert to a data.frame
  combo <- c(ndvi_m, psi_r, tdiff_r)
  df    <- as.data.frame(combo, xy = TRUE, na.rm = TRUE)
  
  ### 7️⃣ add metadata columns
  df$year    <- year
  df$species <- species
  
  message("   ✅ Combined data rows: ", nrow(df))
  return(df)
}


# -----------------------------------------------------------------------------
# Function: process_species
# -----------------------------------------------------------------------------
process_species <- function(species, month_name, redo_ndvi = FALSE) {
  message(sprintf("🎬 [%s] -> %s", month_name, species))
  sp_cfg <- species_config[[species]]
  mo_cfg <- months_config[[month_name]]
  
  out_dir <- file.path("results_monthly_rootzone", month_name, species)
  ensure_directory(out_dir)
  rdata_fp <- file.path(out_dir, sprintf("NDVI_PSI_TDiff_%s_rootzone.RData", gsub("-", "", mo_cfg$day)))
  
  if (file.exists(rdata_fp) && !redo_ndvi) {
    message("🎯 Already processed: ", species, " - skipping.")
    load(rdata_fp)
    return(species_df)
  }
  
  # Only recompute PSI/TDiff if files don't exist
  psi_r_path <- file.path(out_dir, sprintf("psi_%d.tif", min(years)))
  if (!file.exists(psi_r_path)) {
    psi_df <- transfer_psi_to_df(sp_cfg$psi)
    save_psi_raster(psi_df, mo_cfg$day, out_dir)
  }
  
  tdiff_r_path <- file.path(out_dir, sprintf("tdiff_%d.tif", min(years)))
  if (!file.exists(tdiff_r_path)) {
    tdiff_df <- transfer_tdiff_to_df(sp_cfg$tdiff)
    save_tdiff_raster(tdiff_df, mo_cfg$day, out_dir)
  }
  
  # Always redo NDVI and recombine
  lst <- lapply(years, process_NDVI_PSI_TDiff_species_year,
                NDVI_file = mo_cfg$NDVI, species = species, out_dir = out_dir)
  species_df <- bind_rows(lst)
  save(species_df, file = rdata_fp)
  message("📦 Saved species RData: ", rdata_fp)
  return(species_df)
}

# -----------------------------------------------------------------------------
# Function: process_month
# -----------------------------------------------------------------------------
process_month <- function(month_name, redo_ndvi = FALSE) {
  message(sprintf("🗓️ Processing month %s", month_name))
  ensure_directory(file.path("results", "Data", "rootzone"))
  
  # Redo NDVI for each species
  dfs <- lapply(names(species_config), process_species, month_name = month_name, redo_ndvi = redo_ndvi)
  
  month_df <- bind_rows(dfs) %>% mutate(month = month_name)
  out_dir <- file.path("results", "Data", "rootzone")
  save(month_df, file = file.path(out_dir, sprintf("AllSpecies_%s_rootzone.RData", month_name)))
  message("💾 Saved month data: ", month_name)
  return(month_df)
}


# -----------------------------------------------------------------------------
# Main loop
# -----------------------------------------------------------------------------
all_months <- lapply(names(months_config), function(m) process_month(m, redo_ndvi = FALSE))
combined <- bind_rows(all_months) %>% mutate(rootzone = TRUE)
save(combined, file = file.path("results", "Data", "AllSpecies_AllMonths_rootzone.RData"))
message("🎉 Rootzone processing complete! 🎉")
