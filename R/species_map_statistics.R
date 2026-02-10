# ===============================
# Setup
# ===============================
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

library(terra)

# ===============================
# BWI statistics (ha → km²)
# ===============================
bwi <- data.frame(
  Species = c("Oak", "Beech", "Spruce", "Pine"),
  Area_ha = c(
    1265470,  # Oak (Eiche)
    1819006,  # Beech (Buche)
    2296866,  # Spruce (Fichte)
    2396433   # Pine (Kiefer)
  )
)

bwi$Area_km2 <- bwi$Area_ha * 0.01

# ===============================
# Raster-based area calculation
# ===============================
pixel_size_m   <- 231.25
pixel_area_km2 <- (pixel_size_m^2) / 1e6

rasters <- list(
  Oak    = rast("species_map_MODIS/Oak.tif"),
  Beech  = rast("species_map_MODIS/Beech.tif"),
  Spruce = rast("species_map_MODIS/Spruce.tif"),
  Pine   = rast("species_map_MODIS/Pine.tif")
)

calc_area_km2 <- function(r) {
  n_pixels <- global(r == 1, "sum", na.rm = TRUE)[1, 1]
  n_pixels * pixel_area_km2
}

raster_area <- data.frame(
  Species = names(rasters),
  Area_km2_Raster = sapply(rasters, calc_area_km2)
)

# ===============================
# Combine & compare
# ===============================
comparison <- merge(
  bwi[, c("Species", "Area_km2")],
  raster_area,
  by = "Species"
)

names(comparison)[2] <- "Area_km2_BWI"

comparison$Difference_km2 <- comparison$Area_km2_Raster - comparison$Area_km2_BWI
comparison$Difference_pct <- 100 * comparison$Difference_km2 / comparison$Area_km2_BWI

# NEW: raster area as % of BWI
comparison$Raster_pct_of_BWI <- 
  100 * comparison$Area_km2_Raster / comparison$Area_km2_BWI

# ===============================
# Print results
# ===============================
print(comparison, row.names = FALSE)


# total forest area (ha)
area_total_ha <- 10971061

# BWI species areas (ha)
bwi_ha <- data.frame(
  Species = c("Oak", "Beech", "Spruce", "Pine"),
  Area_ha = c(1265470, 1819006, 2296866, 2396433)
)

# percentage of total forest area
bwi_ha$Percentage_total <- 100 * bwi_ha$Area_ha / area_total_ha

# optional: also in km²
bwi_ha$Area_km2 <- bwi_ha$Area_ha * 0.01

print(bwi_ha)


# ============================================================
# Version 2: Reproject rasters to a metric CRS and compute area
# ============================================================

# Choose an equal-area CRS for Europe (recommended for area stats)
crs_area <- "EPSG:3035"  # ETRS89 / LAEA Europe (meters)

# Function: project raster to meters, count ones, compute km²
calc_area_projected_km2 <- function(r, crs_area) {
  r_m <- project(r, crs_area, method = "near")  # keep classes (0/1) intact
  n_pixels <- global(r_m == 1, "sum", na.rm = TRUE)[1, 1]
  pixel_area_km2 <- prod(res(r_m)) / 1e6       # res() in meters, so m² -> km²
  n_pixels * pixel_area_km2
}

# Compute raster areas (km²) after reprojection
raster_area_proj <- data.frame(
  Species = names(rasters),
  Area_km2_Raster_proj = sapply(rasters, calc_area_projected_km2, crs_area = crs_area)
)

# Join with BWI km² and compute percentage
comparison_proj <- merge(
  bwi[, c("Species", "Area_km2")],
  raster_area_proj,
  by = "Species"
)

names(comparison_proj)[2] <- "Area_km2_BWI"
comparison_proj$Raster_pct_of_BWI_proj <-
  100 * comparison_proj$Area_km2_Raster_proj / comparison_proj$Area_km2_BWI

# Print
print(comparison_proj[, c("Species", "Area_km2_BWI", "Area_km2_Raster_proj", "Raster_pct_of_BWI_proj")],
      row.names = FALSE)

