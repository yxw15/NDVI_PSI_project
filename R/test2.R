setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

library(ncdf4)
library(reshape2)
library(ggplot2)
library(terra)

months_config <- list(
  April  = list(month_day = "04-23", NDVI = "../WZMAllDOYs/Quantiles_113.nc"),   # 🌱
  May    = list(month_day = "05-25", NDVI = "../WZMAllDOYs/Quantiles_145.nc"),   # ☀️
  June   = list(month_day = "06-26", NDVI = "../WZMAllDOYs/Quantiles_177.nc"),   # 🌻
  July   = list(month_day = "07-28", NDVI = "../WZMAllDOYs/Quantiles_209.nc"),   # 🌳
  August = list(month_day = "08-29", NDVI = "../WZMAllDOYs/Quantiles_241.nc")    # 🍂
)

# Species-specific settings
species_config <- list(
  Beech = list(
    psi_nc   = "../Allan_Yixuan/PSImean_AMJJA_8days_Bu_bfv_20032024_compressed.nc",
    tdiff_nc = "../Allan_Yixuan/TDiffsum_AMJJA_8days_Bu_bfv20032024_compressed.nc",
    mask     = "results/Species_Maps/Beech_mask.tif"
  ),
  Oak   = list(
    psi_nc   = "../Allan_Yixuan/PSImean_AMJJA_8days_Ei_bfv_20032024_compressed.nc",
    tdiff_nc = "../Allan_Yixuan/TDiffsum_AMJJA_8days_Ei_bfv20032024_compressed.nc",
    mask     = "results/Species_Maps/Oak_mask.tif"
  ),
  Spruce= list(
    psi_nc   = "../Allan_Yixuan/PSImean_AMJJA_8days_Fi_bfv_20032024_compressed.nc",
    tdiff_nc = "../Allan_Yixuan/TDiffsum_AMJJA_8days_Fi_bfv20032024_compressed.nc",
    mask     = "results/Species_Maps/Spruce_mask.tif"
  ),
  Pine  = list(
    psi_nc   = "../Allan_Yixuan/PSImean_AMJJA_8days_Ki_bfv_20032024_compressed.nc",
    tdiff_nc = "../Allan_Yixuan/TDiffsum_AMJJA_8days_Ki_bfv20032024_compressed.nc",
    mask     = "results/Species_Maps/Pine_mask.tif"
  )
)

start_date <- "2003-01-01"
years <- 2003:2024

transfer_psi_to_df <- function(nc_file, start_date) {
  
  # Step 1: Open the NetCDF file
  nc <- nc_open(nc_file)
  
  # Step 2: Extract variables
  time <- ncvar_get(nc, "time")
  depth <- ncvar_get(nc, "depth")
  x <- ncvar_get(nc, "x")
  y <- ncvar_get(nc, "y")
  psi <- ncvar_get(nc, "psi")
  
  # Convert time from "days since" format to actual dates
  dates <- as.Date(start_date) + time
  
  # Close the NetCDF file
  nc_close(nc)
  
  # Step 3: Assign meaningful dimension names
  dimnames(psi) <- list(
    x = x,
    y = y,
    depth = depth,
    time = as.character(dates)
  )
  
  # Step 4: Transform the 4D array into a data frame
  psi_melted <- melt(psi, varnames = c("x", "y", "depth", "time"), value.name = "soil_water_potential")
  
  psi_melted <- na.omit(psi_melted)
  return(psi_melted)
}
save_psi_melted <- function(psi_melted, file_path) {
  write.csv(psi_melted, file_path, row.names = FALSE)
}

filter_psi <- function(psi_melted, month_day, depth) {
  psi_melted <- na.omit(psi_melted)
  psi_filter <- subset(psi_melted, 
                       format(as.Date(time), "%m-%d") == month_day & depth == depth)
  return(psi_filter)
}

save_psi_raster <- function(psi_melted, month_day, depth, output_dir) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  psi <- filter_psi(psi_melted, month_day, depth)
  psi$time <- as.Date(psi$time, format = "%Y-%m-%d")
  years <- 2003:2024
  for (year in years) {
    psi_year <- psi %>% filter(format(time, "%Y") == as.character(year))
    psi_selected <- psi_year[, c("x", "y", "soil_water_potential")]
    psi_raster <- rast(psi_selected, type = "xyz")
    crs(psi_raster) <- "epsg:31467"
    output_path <- file.path(output_dir, paste0("psi_", year, ".tif"))
    writeRaster(psi_raster, output_path, overwrite = TRUE)
    cat("Raster saved for year:", year, "at", output_path, "\n")
  }
}

depth_val = 100
month_name = "August"
month_day = "08-29"
species_name = "Beech"
Beech.df <- transfer_psi_to_df(species_config$Beech$psi_nc, start_date)
Beech.df.filter <- filter_psi(Beech.df, month_day, depth)
  
  
output_dir = file.path(sprintf("results_monthly_%d", depth_val), month_name, species_name)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

save_psi_raster(psi_melted = df, month_day = "08-29", depth = depth_val, output_dir)



# Check PSI at same pixels for different depth at same time

Beech.day <- subset(Beech.df, time == "2003-04-15")
set.seed(123)  
unique_pixels_day <- unique(Beech.day[, c("x", "y")])
sampled_pixels_day <- unique_pixels_day[sample(1:nrow(unique_pixels_day), 10), ]

Beech.sampled.day <- Beech.day[Beech.day$x %in% sampled_pixels_day$x & Beech.day$y %in% sampled_pixels_day$y, ]

library(ggplot2)

ggplot(Beech.sampled.day, aes(x = depth, y = soil_water_potential, color = factor(paste(x, y)))) +
  geom_line() +
  geom_point() +
  scale_x_reverse() +  
  labs(x = "Depth (cm)", y = "Soil Water Potential", color = "Pixel (x, y)") +
  theme_minimal()

### test results ###
library(terra)
library(dplyr)
library(ggplot2)
library(viridis)

# 1) Load your two rasters
r50  <- rast("results_monthly_50/April/Beech/psi_2024.tif")
r100 <- rast("results_monthly_100/April/Beech/psi_2024.tif")

par(mfrow = c(1, 2), 
    mar    = c(5, 4, 4, 2))
plot(r50,
     main = "Depth 50")
plot(r100,
     main = "Depth 100")
par(mfrow = c(1, 1))

# 2) Convert each to a data.frame with x, y, and the psi value
df50  <- as.data.frame(r50,  xy = TRUE) %>%
  rename(psi = 3) %>%
  mutate(depth = "Depth = 50 cm")

df100 <- as.data.frame(r100, xy = TRUE) %>%
  rename(psi = 3) %>%
  mutate(depth = "Depth = 100 cm")

# 3) Combine and compute global range
df  <- bind_rows(df50, df100)
rng <- range(df$psi, na.rm = TRUE)

# 4) Plot with viridis, white background, and panel borders
ggplot(df, aes(x = x, y = y, fill = psi)) +
  geom_raster() +
  facet_wrap(~ depth, ncol = 2) +
  coord_equal() +
  scale_fill_viridis_c(
    option = "viridis",
    limits = rng,
    name   = expression(psi~"(2024)")
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    # white plot background
    plot.background    = element_rect(fill = "white", color = NA),
    panel.background   = element_rect(fill = "white", color = NA),
    # black border around each facet
    panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.5),
    # remove grid lines
    panel.grid         = element_blank(),
    # adjust legend bar thickness
    legend.key.height  = unit(2, "cm"),
    legend.title       = element_text(size = 11),
    legend.text        = element_text(size = 9),
    strip.text         = element_text(size = 12)
  )

