# 📚 Load required libraries
library(terra)
library(dplyr)

# 📂 Set working directory
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

# 🔧 Global parameters
years       <- 2003:2024
months_cfg  <- list(
  April  = list(day="04-23", NDVI="../WZMAllDOYs/Quantiles_113.nc"),
  May    = list(day="05-25", NDVI="../WZMAllDOYs/Quantiles_145.nc"),
  June   = list(day="06-26", NDVI="../WZMAllDOYs/Quantiles_177.nc"),
  July   = list(day="07-28", NDVI="../WZMAllDOYs/Quantiles_209.nc"),
  August = list(day="08-29", NDVI="../WZMAllDOYs/Quantiles_241.nc")
)
species_cfg <- list(
  Beech  = list(
    psi   = "../ALLAN_PIA_SoilMoisture/PSIM/PSImean_rootzone_AMJJA_8days_Bu_bfv_20032024_compressed.nc",
    tdiff = "../ALLAN_PIA_SoilMoisture/TDIFF/TDiffsum_AMJJA_8days_Bu_bfv20032024_compressed.nc",
    mask  = "species_map_MODIS/Beech.tif"
  ),
  Oak    = list(
    psi   = "../ALLAN_PIA_SoilMoisture/PSIM/PSImean_rootzone_AMJJA_8days_Ei_bfv_20032024_compressed.nc",
    tdiff = "../ALLAN_PIA_SoilMoisture/TDIFF/TDiffsum_AMJJA_8days_Ei_bfv20032024_compressed.nc",
    mask  = "species_map_MODIS/Oak.tif"
  ),
  Spruce = list(
    psi   = "../ALLAN_PIA_SoilMoisture/PSIM/PSImean_rootzone_AMJJA_8days_Fi_bfv_20032024_compressed.nc",
    tdiff = "../ALLAN_PIA_SoilMoisture/TDIFF/TDiffsum_AMJJA_8days_Fi_bfv20032024_compressed.nc",
    mask  = "species_map_MODIS/Spruce.tif"
  ),
  Pine   = list(
    psi   = "../ALLAN_PIA_SoilMoisture/PSIM/PSImean_rootzone_AMJJA_8days_Ki_bfv_20032024_compressed.nc",
    tdiff = "../ALLAN_PIA_SoilMoisture/TDIFF/TDiffsum_AMJJA_8days_Ki_bfv20032024_compressed.nc",
    mask  = "species_map_MODIS/Pine.tif"
  )
)

# -----------------------------------------------------------------------------
# Function: process one species × one month
# -----------------------------------------------------------------------------
process_species_month <- function(species, month_name) {
  mo    <- months_cfg[[month_name]]
  sp    <- species_cfg[[species]]
  maskr <- rast(sp$mask)
  
  # read all three as multi-layer SpatRasters
  ndvi  <- rast(mo$NDVI)
  psi   <- rast(sp$psi)     # variable 'psi'
  tdiff <- rast(sp$tdiff)   # variable 'tdiff'
  
  # extract their time vectors as Dates
  ndates <- as.Date(time(ndvi))
  pdates <- as.Date(time(psi))
  tdates <- as.Date(time(tdiff))
  
  out_list <- list()
  
  for (yr in years) {
    target_date <- as.Date(sprintf("%d-%s", yr, mo$day))
    i_ndvi <- which(ndates == target_date)
    i_psi  <- which(pdates == target_date)
    i_tdiff<- which(tdates == target_date)
    
    if (length(i_ndvi)==0 || length(i_psi)==0 || length(i_tdiff)==0) next
    
    # single‐layer rasters for this date
    r_ndvi <- ndvi[[i_ndvi]]
    r_psi  <- psi[[i_psi ]]
    r_td   <- tdiff[[i_tdiff]]
    
    # reproject & mask
    r_ndvi <- project(r_ndvi, maskr, method="near") %>% mask(maskr)
    r_psi  <- project(r_psi,  maskr, method="bilinear") %>% mask(maskr)
    r_td   <- project(r_td,   maskr, method="bilinear") %>% mask(maskr)
    
    # stack & to data.frame
    stk <- c(r_ndvi, r_psi, r_td)
    names(stk) <- c("Quantiles","soil_water_potential","transpiration_deficit")
    
    df <- as.data.frame(stk, xy=TRUE, na.rm=TRUE)
    df$year    <- yr
    df$species <- species
    df$month   <- month_name
    
    out_list[[as.character(yr)]] <- df
  }
  
  bind_rows(out_list)
}

# -----------------------------------------------------------------------------
# Loop over all species & months, combine into one big data.frame
# -----------------------------------------------------------------------------
all_dfs <- list()
for (month in names(months_cfg)) {
  for (sp in names(species_cfg)) {
    message("Processing ", sp, " in ", month)
    df <- process_species_month(sp, month)
    all_dfs[[paste(sp, month, sep = "_")]] <- df
  }
}
combined_df <- bind_rows(all_dfs)

# 👀 Inspect
head(combined_df)
save(combined, file = file.path("results", "Data", "AllSpecies_AllMonths_rootzone_terra.RData"))