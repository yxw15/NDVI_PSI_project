# -----------------------------------------------
# R script: Generate soil‐under‐species rasters in MODIS space
#            & compute percentage of each soil type under each species
# -----------------------------------------------
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

# 1. Load required libraries
library(terra)
library(dplyr)
library(ggplot2)

# 2. Load soil stack and TreeSpecies MODIS stack
soil_code_stack <- rast("soil_map/soil_reclass_layers.tif")
species.MODIS    <- rast("../TreeSpeciesGermany/TreeSpeciesMODIS.tif")

# 3. Extract and threshold each species band (threshold > 0.7 → presence = 1; else NA)
Beech.MODIS  <- species.MODIS[[2]]
Oak.MODIS    <- species.MODIS[[4]]
Spruce.MODIS <- species.MODIS[[6]]
Pine.MODIS   <- species.MODIS[[7]]

species_list <- list(
  Beech  = Beech.MODIS,
  Oak    = Oak.MODIS,
  Spruce = Spruce.MODIS,
  Pine   = Pine.MODIS
)

species_bin_MODIS <- list()
for (sp in names(species_list)) {
  r <- species_list[[sp]]
  r[r > 0.7]  <- 1
  r[r <= 0.7] <- NA
  names(r) <- sp
  species_bin_MODIS[[sp]] <- r
}

# 4. Define the soil “template” to get CRS & grid
#    Use the entire soil stack (all layers) as the template
soil_orig <- soil_code_stack

# 5. Project each species binary mask onto the soil grid (all layers share same CRS/grid)
species_bin_soil <- list()
for (sp in names(species_bin_MODIS)) {
  sp_modis    <- species_bin_MODIS[[sp]]
  sp_on_soil  <- project(sp_modis, soil_orig, method = "near")
  names(sp_on_soil) <- sp
  species_bin_soil[[sp]] <- sp_on_soil
}

# 6. Loop over each species to:
#    (a) generate "soil_under_<species>_MODIS.tif"
#    (b) compute percentage of each soil layer under that species
all_species_stats <- list()

for (sp in names(species_bin_soil)) {
  message("Processing species: ", sp, " …")
  
  # (a) Species presence mask in soil CRS (TRUE where 1, NA elsewhere)
  sp_soil_mask <- species_bin_soil[[sp]] == 1
  
  # (b) Mask the entire soil stack by that species presence
  soil_masked <- mask(soil_code_stack, sp_soil_mask)
  
  # (c) Reproject masked soil stack back to MODIS grid
  soil_masked_MODIS <- project(
    x      = soil_masked,
    y      = species_bin_MODIS[[sp]],
    method = "near"
  )
  
  # (d) Write output raster: soil_map/<Species>_soil_MODIS.tif
  out_raster_fname <- file.path("soil_map", paste0(sp, "_soil_MODIS.tif"))
  writeRaster(
    soil_masked_MODIS,
    filename  = out_raster_fname,
    overwrite = TRUE
  )
  message(" → Saved raster: ", out_raster_fname)
  
  # (e) Compute statistics in soil CRS
  total_pixels <- global(!is.na(species_bin_soil[[sp]]), "sum", na.rm = TRUE)[1,1]
  
  df <- data.frame(
    layer   = names(soil_code_stack),
    percent = numeric(nlyr(soil_code_stack)),
    stringsAsFactors = FALSE
  )
  for (i in seq_along(df$layer)) {
    valid_pixels <- global(!is.na(soil_masked[[i]]), "sum", na.rm = TRUE)[1,1]
    df$percent[i] <- (valid_pixels / total_pixels) * 100
  }
  df$species <- sp
  all_species_stats[[sp]] <- df
}

# 7. Combine all species statistics into one data frame
results_all <- bind_rows(all_species_stats)

# 8. Read soil descriptions and join
soil_descrip <- read.csv("soil_map/feinbod_lookup_with_english.csv",
                         stringsAsFactors = FALSE)

results_all <- results_all %>%
  left_join(soil_descrip, by = c("layer" = "feinbod")) %>%
  mutate(percent = round(percent, 2)) %>%
  select(species, feinbod_descr_en, percent)

# 9. Set plotting order and color palette
species_order <- c("Oak", "Beech", "Spruce", "Pine")
cb_palette <- c(
  Oak    = "#E69F00",  # Orange
  Beech  = "#0072B2",  # Deep blue
  Spruce = "#009E73",  # Bluish-green
  Pine   = "#F0E442"   # Yellow
)
results_all$species <- factor(results_all$species, levels = species_order)


# 11. Loop over species to generate and save bar plots
for (sp in levels(results_all$species)) {
  species_df <- results_all %>% filter(species == sp)
  
  p <- ggplot(species_df, aes(x = reorder(feinbod_descr_en, percent), y = percent)) +
    geom_bar(stat = "identity", fill = cb_palette[sp]) +
    coord_flip() +
    labs(
      x     = "soil type",
      y     = "percent of pixels (%)",
      title = ""
    ) +
    theme_minimal() +
    theme(
      axis.text.x       = element_text(angle = 0, hjust = 0.5),
      plot.background   = element_rect(fill = "white", color = "white"),
      panel.background  = element_rect(fill = "white"),
      legend.position   = "none",
      plot.title        = element_text(hjust = 0.5, size = 18, face = "bold"),
      axis.title        = element_text(face = "bold", size = 14),
      axis.text         = element_text(color = "black", size = 12),
      panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank()
    )
  
  print(p)
  
  ggsave(
    filename = file.path("soil_map", paste0("soil_statistics_", sp, ".png")),
    plot     = p,
    width    = 9,
    height   = 6,
    dpi      = 300
  )
  message(" → Saved plot: soil_map/soil_statistics_", sp, ".png\n")
}

# 12. (Optional) Example plot for Beech – first soil layer under Beech in MODIS space
plot(
  rast("soil_map/Beech_soil_MODIS.tif")[[1]],
  main = "Soil Layer 1 under Beech (in MODIS grid)"
)


# -----------------------------------------------
# R script: Extract soil_code values under each species
#            (i.e., create a MODIS‐resolution raster of soil codes for each species)
# -----------------------------------------------

# 0. Set working directory (adjust as needed)
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

# 1. Load required libraries
library(terra)
library(dplyr)

# 2. Load the single‐layer soil_code raster (EPSG:31467, 1 km resolution)
soil_code <- rast("soil_map/soil_code.tif")
#   Values range from 1 to 34; each integer represents a soil type.

# 3. Load the TreeSpecies MODIS stack (EPSG:4326, ~0.0023° resolution)
species.MODIS <- rast("../TreeSpeciesGermany/TreeSpeciesMODIS.tif")

# 4. Extract and threshold each species band (threshold > 0.7 → presence = 1; else NA)
Beech.MODIS  <- species.MODIS[[2]]
Oak.MODIS    <- species.MODIS[[4]]
Spruce.MODIS <- species.MODIS[[6]]
Pine.MODIS   <- species.MODIS[[7]]

species_list <- list(
  Beech  = Beech.MODIS,
  Oak    = Oak.MODIS,
  Spruce = Spruce.MODIS,
  Pine   = Pine.MODIS
)

species_bin_MODIS <- list()
for (sp in names(species_list)) {
  r <- species_list[[sp]]
  r[r > 0.7]  <- 1
  r[r <= 0.7] <- NA
  names(r) <- sp
  species_bin_MODIS[[sp]] <- r
}

# 5. Project each species binary mask onto the soil_code grid (EPSG:31467, 1 km)
#    So that we can mask soil_code in its native CRS/resolution.
species_bin_soil <- list()
for (sp in names(species_bin_MODIS)) {
  sp_modis   <- species_bin_MODIS[[sp]]
  sp_on_soil <- project(sp_modis, soil_code, method = "near")
  names(sp_on_soil) <- sp
  species_bin_soil[[sp]] <- sp_on_soil
}

# 6. Loop over each species to:
#    (a) mask soil_code by the species presence (in soil CRS)
#    (b) reproject that masked soil_code back to MODIS grid
#    (c) save the resulting raster of soil codes under the species
soil_code_under_species_MODIS <- list()

for (sp in names(species_bin_soil)) {
  message("Processing species: ", sp, " …")
  
  # (a) Species presence mask in soil CRS (TRUE where 1, NA elsewhere)
  sp_soil_mask <- species_bin_soil[[sp]] == 1
  
  # (b) Mask soil_code by that presence mask
  soil_masked <- mask(soil_code, sp_soil_mask)
  #    → soil_masked has integer soil_code values where species is present, NA elsewhere
  
  # (c) Reproject masked soil_code back to MODIS grid (EPSG:4326, ~0.0023°)
  soil_masked_MODIS <- project(
    x      = soil_masked,
    y      = species_bin_MODIS[[sp]],
    method = "near"
  )
  names(soil_masked_MODIS) <- sp
  
  # (d) Save output raster: soil_map/<Species>_soilCode_MODIS.tif
  out_fname <- file.path("soil_map", paste0(sp, "_soilCode_MODIS.tif"))
  writeRaster(
    soil_masked_MODIS,
    filename  = out_fname,
    overwrite = TRUE
  )
  message(" → Saved raster: ", out_fname)
  
  # (e) Store in list for any further use
  soil_code_under_species_MODIS[[sp]] <- soil_masked_MODIS
}

# 8. (Optional) Summarize soil_code frequency under each species, in MODIS CRS
#    This counts how many MODIS pixels carry each soil_code per species.
soil_frequencies <- list()

for (sp in names(soil_code_under_species_MODIS)) {
  r <- soil_code_under_species_MODIS[[sp]]
  freq_table <- freq(r, digits = 0, value = TRUE)
  # 'freq' returns a two‐column matrix: value (soil_code) and count (number of pixels)
  freq_df <- as.data.frame(freq_table)
  colnames(freq_df) <- c("soil_code", "count")
  freq_df$species <- sp
  soil_frequencies[[sp]] <- freq_df
}

# Combine into one data frame
soil_freq_all <- bind_rows(soil_frequencies) %>%
  arrange(species, soil_code)

# Print summary of counts
print(soil_freq_all)

# -----------------------------------------------
# At this point:
# • You have four GeoTIFFs in "soil_map/":
#     - Beech_soilCode_MODIS.tif
#     - Oak_soilCode_MODIS.tif
#     - Spruce_soilCode_MODIS.tif
#     - Pine_soilCode_MODIS.tif
#   Each is a single‐layer SpatRaster in MODIS CRS/resolution whose values
#   are the integer soil_code where the species was present, and NA elsewhere.
#
# • `soil_freq_all` is a data.frame showing, for each species, how many MODIS pixels
#   carry each soil_code.  If a particular soil_code does not appear under a species,
#   it will simply be absent from that species’s frequency table.
# -----------------------------------------------

# -----------------------------------------------
# R code: Plot Clay, Silt, and Sand Percentages for Beech over Germany using ggplot2
# -----------------------------------------------

# 1. Load required libraries
library(terra)
library(dplyr)
library(ggplot2)
library(gridExtra)   # for arranging multiple ggplots

# 2. Read in Beech‐specific soil‐code raster (MODIS CRS, 0.0023° resolution)
Beech_soil_code <- rast("soil_map/Beech_soilCode_MODIS.tif")

# 3. Read in soil descriptions and soil texture classes
soil_descrip <- read.csv(
  "soil_map/feinbod_lookup_with_english.csv",
  stringsAsFactors = FALSE
)

soil_texture <- read.csv(
  "soil_map/soil_texture_classes.csv",
  stringsAsFactors = FALSE
)

# 4. Build lookup table: feinbod_code → soil‐code string → midpoint clay/silt/sand
lut <- soil_descrip %>%
  select(feinbod_code, feinbod) %>%
  inner_join(
    soil_texture,
    by = c("feinbod" = "code")
  ) %>%
  mutate(
    clay_mid = (clay_lower + clay_upper) / 2,
    silt_mid = (silt_lower + silt_upper) / 2,
    sand_mid = (sand_lower + sand_upper) / 2
  ) %>%
  select(feinbod_code, clay_mid, silt_mid, sand_mid)

# 5. Use terra::subs() to replace each feinbod_code with the corresponding midpoint
#    This avoids potential issues with classify() and works directly on the integer codes.

# (a) Clay midpoint raster
clay_raster <- subst(
  Beech_soil_code,
  from = lut$feinbod_code,
  to   = lut$clay_mid
)
names(clay_raster) <- "clay"

# (b) Silt midpoint raster
silt_raster <- subst(
  Beech_soil_code,
  from = lut$feinbod_code,
  to   = lut$silt_mid
)
names(silt_raster) <- "silt"

# (c) Sand midpoint raster
sand_raster <- subst(
  Beech_soil_code,
  from = lut$feinbod_code,
  to   = lut$sand_mid
)
names(sand_raster) <- "sand"

# 6. Convert each raster to a data frame for ggplot
#    We include only non‐NA pixels (i.e., where Beech is present)
clay_df <- as.data.frame(clay_raster, xy = TRUE) %>%
  rename(value = clay) %>%
  filter(!is.na(value))

silt_df <- as.data.frame(silt_raster, xy = TRUE) %>%
  rename(value = silt) %>%
  filter(!is.na(value))

sand_df <- as.data.frame(sand_raster, xy = TRUE) %>%
  rename(value = sand) %>%
  filter(!is.na(value))

# 7. Create ggplot objects

library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)
library(gridExtra)  # For arranging plots

boundary_germany <- ne_countries(scale = "medium", 
                                 country = "Germany", 
                                 returnclass = "sf")

## 1. CLAY: white → blue
p_clay <- ggplot(clay_df, aes(x = x, y = y, color = value)) +
  geom_point(size = 0.5) +
  geom_sf(data = boundary_germany, fill = NA, color = "black", inherit.aes = FALSE) +
  coord_sf() +
  scale_color_gradientn(
    colours = c("white", "lightblue", "dodgerblue", "#0072B2"),
    name    = "clay (%)"
  ) +
  labs(title = "", x = NULL, y = NULL) +
  theme(
    legend.key.width = unit(1.3, "cm"),  # or try 2.5–3 for even wider
    legend.margin = margin(t = 0, r = 20, b = 0, l = 20),  # adds left/right space
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.background = element_rect(fill = "white", color = "white"),
    panel.background = element_rect(fill = "white"),
    legend.background = element_rect(fill = "white", color = "white"),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_text(color = "black", size = 14),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    strip.text = element_text(face = "bold", size = 12)
  )
print(p_clay)

## 2. SILT: white → green
p_silt <- ggplot(silt_df, aes(x = x, y = y, color = value)) +
  geom_point(size = 0.5) +
  geom_sf(data = boundary_germany, fill = NA, color = "black", inherit.aes = FALSE) +
  coord_sf() +
  scale_color_gradientn(
    colours = c("white", "#b9f6ca", "#00bfae", "#00675b"),
    name    = "silt (%)"
  ) +
  labs(title = "", x = NULL, y = NULL) +
  theme(
    legend.key.width = unit(1.3, "cm"),  # or try 2.5–3 for even wider
    legend.margin = margin(t = 0, r = 20, b = 0, l = 20),  # adds left/right space
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.background = element_rect(fill = "white", color = "white"),
    panel.background = element_rect(fill = "white"),
    legend.background = element_rect(fill = "white", color = "white"),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_text(color = "black", size = 14),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    strip.text = element_text(face = "bold", size = 12)
  )

## 3. SAND: white → orange
p_sand <- ggplot(sand_df, aes(x = x, y = y, color = value)) +
  geom_point(size = 0.5) +
  geom_sf(data = boundary_germany, fill = NA, color = "black", inherit.aes = FALSE) +
  coord_sf() +
  scale_color_gradientn(
    colours = c("white", "#ffe0b2", "orange", "#ff6600"),
    name    = "sand (%)"
  ) +
  labs(title = "", x = NULL, y = NULL) +
  theme(
    legend.key.width = unit(1.3, "cm"),  # or try 2.5–3 for even wider
    legend.margin = margin(t = 0, r = 20, b = 0, l = 20),  # adds left/right space
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.background = element_rect(fill = "white", color = "white"),
    panel.background = element_rect(fill = "white"),
    legend.background = element_rect(fill = "white", color = "white"),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_text(color = "black", size = 14),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    strip.text = element_text(face = "bold", size = 12)
  )

## Arrange all three plots side by side
# grid.arrange(p_clay, p_silt, p_sand, ncol = 3)
library(patchwork)
soil_grid <- p_clay + p_silt + p_sand + plot_layout(ncol = 3)
ggsave("soil_map/Beech_soil_fractions_grid.png", plot = soil_grid, width = 12, height = 5, dpi = 300)

## (Optional) Save each plot
ggsave("soil_map/Beech_clay_percentage.png", plot = p_clay, width = 5, height = 5, dpi = 300)
ggsave("soil_map/Beech_silt_percentage.png", plot = p_silt, width = 5, height = 5, dpi = 300)
ggsave("soil_map/Beech_sand_percentage.png", plot = p_sand, width = 5, height = 5, dpi = 300)


library(terra)
library(dplyr)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)

# Boundary once (doesn't need to be in function)
boundary_germany <- ne_countries(scale = "medium", 
                                 country = "Germany", 
                                 returnclass = "sf")

# The plotting function
make_soil_fraction_plots <- function(
    species, 
    raster_path, 
    lut, 
    out_dir = "soil_map"
) {
  # 1. Read raster
  soil_code <- rast(raster_path)
  
  # 2. Map codes to fractions
  clay_raster <- subst(soil_code, from = lut$feinbod_code, to = lut$clay_mid)
  names(clay_raster) <- "clay"
  silt_raster <- subst(soil_code, from = lut$feinbod_code, to = lut$silt_mid)
  names(silt_raster) <- "silt"
  sand_raster <- subst(soil_code, from = lut$feinbod_code, to = lut$sand_mid)
  names(sand_raster) <- "sand"
  
  # 3. To data.frame
  clay_df <- as.data.frame(clay_raster, xy = TRUE) %>% rename(value = clay) %>% filter(!is.na(value))
  silt_df <- as.data.frame(silt_raster, xy = TRUE) %>% rename(value = silt) %>% filter(!is.na(value))
  sand_df <- as.data.frame(sand_raster, xy = TRUE) %>% rename(value = sand) %>% filter(!is.na(value))
  
  # 4. Plot objects (style is identical to Beech above)
  p_clay <- ggplot(clay_df, aes(x = x, y = y, color = value)) +
    geom_point(size = 0.5) +
    geom_sf(data = boundary_germany, fill = NA, color = "black", inherit.aes = FALSE) +
    coord_sf() +
    scale_color_gradientn(
      colours = c("white", "lightblue", "dodgerblue", "#0072B2"),
      name    = "clay (%)"
    ) +
    labs(title = "", x = NULL, y = NULL) +
    theme(
      legend.key.width = unit(1.3, "cm"),
      legend.margin = margin(t = 0, r = 20, b = 0, l = 20),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(color = "black", size = 14),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    )
  
  p_silt <- ggplot(silt_df, aes(x = x, y = y, color = value)) +
    geom_point(size = 0.5) +
    geom_sf(data = boundary_germany, fill = NA, color = "black", inherit.aes = FALSE) +
    coord_sf() +
    scale_color_gradientn(
      colours = c("white", "#b9f6ca", "#00bfae", "#00675b"),
      name    = "silt (%)"
    ) +
    labs(title = "", x = NULL, y = NULL) +
    theme(
      legend.key.width = unit(1.3, "cm"),
      legend.margin = margin(t = 0, r = 20, b = 0, l = 20),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(color = "black", size = 14),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    )
  
  p_sand <- ggplot(sand_df, aes(x = x, y = y, color = value)) +
    geom_point(size = 0.5) +
    geom_sf(data = boundary_germany, fill = NA, color = "black", inherit.aes = FALSE) +
    coord_sf() +
    scale_color_gradientn(
      colours = c("white", "#ffe0b2", "orange", "#ff6600"),
      name    = "sand (%)"
    ) +
    labs(title = "", x = NULL, y = NULL) +
    theme(
      legend.key.width = unit(1.3, "cm"),
      legend.margin = margin(t = 0, r = 20, b = 0, l = 20),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(color = "black", size = 14),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    )
  
  # Combined grid
  soil_grid <- p_clay + p_silt + p_sand + plot_layout(ncol = 3, guides = "collect") &
    theme(legend.position = "top")
  
  # 5. Save
  ggsave(file.path(out_dir, paste0(species, "_soil_fractions_grid.png")), plot = soil_grid, width = 12, height = 5, dpi = 300)
  ggsave(file.path(out_dir, paste0(species, "_clay_percentage.png")), plot = p_clay, width = 5, height = 5, dpi = 300)
  ggsave(file.path(out_dir, paste0(species, "_silt_percentage.png")), plot = p_silt, width = 5, height = 5, dpi = 300)
  ggsave(file.path(out_dir, paste0(species, "_sand_percentage.png")), plot = p_sand, width = 5, height = 5, dpi = 300)
}

# 6. Run for all species
species_list <- c("Beech", "Oak", "Pine", "Spruce")
for (species in species_list) {
  raster_path <- file.path("soil_map", paste0(species, "_soilCode_MODIS.tif"))
  make_soil_fraction_plots(
    species = species,
    raster_path = raster_path,
    lut = lut,
    out_dir = "soil_map"
  )
}


# -----------------------------------------------
# R code: Composite figure with columns = species (Oak, Beech, Spruce, Pine),
#          rows = fractions (Clay, Silt, Sand), species names atop each column,
#          single shared legend at the bottom, and no grey border around panels
# -----------------------------------------------

# 0. Set working directory
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

# 1. Load required libraries
library(terra)
library(dplyr)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)

# 2. Read Germany boundary once
boundary_germany <- ne_countries(
  scale = "medium",
  country = "Germany",
  returnclass = "sf"
)

# 3. Read LUT (feinbod_code → clay_mid, silt_mid, sand_mid)
soil_descrip <- read.csv(
  "soil_map/feinbod_lookup_with_english.csv",
  stringsAsFactors = FALSE
)
soil_texture <- read.csv(
  "soil_map/soil_texture_classes.csv",
  stringsAsFactors = FALSE
)

lut <- soil_descrip %>%
  select(feinbod_code, feinbod) %>%
  inner_join(
    soil_texture,
    by = c("feinbod" = "code")
  ) %>%
  mutate(
    clay_mid = (clay_lower + clay_upper) / 2,
    silt_mid = (silt_lower + silt_upper) / 2,
    sand_mid = (sand_lower + sand_upper) / 2
  ) %>%
  select(feinbod_code, clay_mid, silt_mid, sand_mid)

# 4. Function to create clay/silt/sand plots (no individual legends, white backgrounds)
make_fraction_plots <- function(raster_path) {
  # (a) Load the soil_code raster
  soil_code <- rast(raster_path)
  
  # (b) Reclassify into clay/silt/sand midpoints
  clay_r <- subst(soil_code, from = lut$feinbod_code, to = lut$clay_mid)
  names(clay_r) <- "clay"
  silt_r <- subst(soil_code, from = lut$feinbod_code, to = lut$silt_mid)
  names(silt_r) <- "silt"
  sand_r <- subst(soil_code, from = lut$feinbod_code, to = lut$sand_mid)
  names(sand_r) <- "sand"
  
  # (c) Convert each to data.frame & drop NA
  clay_df <- as.data.frame(clay_r, xy = TRUE) %>%
    rename(value = clay) %>%
    filter(!is.na(value))
  silt_df <- as.data.frame(silt_r, xy = TRUE) %>%
    rename(value = silt) %>%
    filter(!is.na(value))
  sand_df <- as.data.frame(sand_r, xy = TRUE) %>%
    rename(value = sand) %>%
    filter(!is.na(value))
  
  # (d) Base theme for all fraction plots with white background and no panel border color
  base_theme <- theme_minimal() +
    theme(
      axis.text.x        = element_text(angle = 0, hjust = 0.5),
      plot.background    = element_rect(fill = "white", color = "white"),
      panel.background   = element_rect(fill = "white"),
      legend.background  = element_rect(fill = "white", color = "white"),
      plot.title         = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      plot.subtitle      = element_text(hjust = 0.5, size = 14),
      axis.title         = element_blank(),
      axis.text          = element_text(color = "black", size = 14),
      panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major   = element_blank(),
      panel.grid.minor   = element_blank(),
      legend.position    = "bottom",
      legend.key.width   = unit(1.3, "cm"),    # ← make each legend key wider
      legend.key.height  = unit(0.5, "cm"),  # ← optionally make it taller as well
      legend.text        = element_text(size = 14),
      legend.title        = element_text(size = 14),
      strip.background   = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text         = element_text(face = "bold", size = 12)
    )
  
  # (e) Create clay plot
  p_clay <- ggplot(clay_df, aes(x = x, y = y, color = value)) +
    geom_point(size = 0.5) +
    geom_sf(data = boundary_germany, fill = NA, color = "black", inherit.aes = FALSE) +
    coord_sf(expand = FALSE) +
    scale_color_gradientn(
      colours = c("white", "lightblue", "dodgerblue", "#0072B2"),
      name    = "clay (%)"
    ) +
    base_theme
  
  # (f) Create silt plot
  p_silt <- ggplot(silt_df, aes(x = x, y = y, color = value)) +
    geom_point(size = 0.5) +
    geom_sf(data = boundary_germany, fill = NA, color = "black", inherit.aes = FALSE) +
    coord_sf(expand = FALSE) +
    scale_color_gradientn(
      colours = c("white", "#b9f6ca", "#00bfae", "#00675b"),
      name    = "silt (%)"
    ) +
    base_theme
  
  # (g) Create sand plot
  p_sand <- ggplot(sand_df, aes(x = x, y = y, color = value)) +
    geom_point(size = 0.5) +
    geom_sf(data = boundary_germany, fill = NA, color = "black", inherit.aes = FALSE) +
    coord_sf(expand = FALSE) +
    scale_color_gradientn(
      colours = c("white", "#ffe0b2", "orange", "#ff6600"),
      name    = "sand (%)"
    ) +
    base_theme
  
  return(list(clay = p_clay, silt = p_silt, sand = p_sand))
}

# 5. Generate plots for each species in order: Oak, Beech, Spruce, Pine
species_order <- c("Oak", "Beech", "Spruce", "Pine")
clay_plots <- list()
silt_plots <- list()
sand_plots <- list()

for (sp in species_order) {
  path <- file.path("soil_map", paste0(sp, "_soilCode_MODIS.tif"))
  p_list <- make_fraction_plots(path)
  
  # Add species name as column title (only on clay plots, top row)
  clay_plots[[sp]] <- p_list$clay +
    labs(title = sp) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.margin = margin(t = 5, r = 5, b = 2, l = 5)
    )
  
  silt_plots[[sp]] <- p_list$silt + theme(plot.margin = margin(t = 2, r = 5, b = 2, l = 5))
  sand_plots[[sp]] <- p_list$sand + theme(plot.margin = margin(t = 2, r = 5, b = 5, l = 5))
}

# 6. Combine into a single list in row-major order:
#    Row 1: clay for Oak, Beech, Spruce, Pine
#    Row 2: silt for Oak, Beech, Spruce, Pine
#    Row 3: sand for Oak, Beech, Spruce, Pine
all_plots <- c(
  clay_plots[species_order],
  silt_plots[species_order],
  sand_plots[species_order]
)

# 7. Arrange with patchwork: 3 rows × 4 columns, collect a shared legend at bottom,
#    and set overall plot background to white (removing any grey border)
composite <- wrap_plots(all_plots, ncol = 4, nrow = 3, guides = "collect") &
  theme(
    legend.position  = "bottom",
    plot.background   = element_rect(fill = "white", color = NA),
    plot.margin       = margin(5, 5, 5, 5),      # small margin on all sides
    panel.spacing     = unit(0, "cm")
  )

# 8. Display the composite
print(composite)

# 9. Save the final composite figure
ggsave(
  filename = "soil_map/all_species_soil_fractions_composite.png",
  plot     = composite,
  width    = 12,   # 4 columns × 3 inches each ≈ 12″ wide
  height   = 12,    # 3 rows × 3 inches each  ≈ 9″ tall
  dpi      = 300,
  bg       = "white"
)

# -----------------------------------------------
# Composite + Standardized Bar Plot of Soil Texture Fractions
# -----------------------------------------------

# 0. Set working directory
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

# 1. Load libraries
library(terra)
library(dplyr)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)

# 2. Read Germany boundary
boundary_germany <- ne_countries(
  scale = "medium",
  country = "Germany",
  returnclass = "sf"
)

# 3. Read LUT
soil_descrip <- read.csv("soil_map/feinbod_lookup_with_english.csv", stringsAsFactors = FALSE)
soil_texture <- read.csv("soil_map/soil_texture_classes.csv", stringsAsFactors = FALSE)

lut <- soil_descrip %>%
  select(feinbod_code, feinbod) %>%
  inner_join(soil_texture, by = c("feinbod" = "code")) %>%
  mutate(
    clay_mid = (clay_lower + clay_upper) / 2,
    silt_mid = (silt_lower + silt_upper) / 2,
    sand_mid = (sand_lower + sand_upper) / 2
  ) %>%
  select(feinbod_code, clay_mid, silt_mid, sand_mid)

# 4. Plot + data extraction function
make_fraction_plots <- function(raster_path) {
  soil_code <- rast(raster_path)
  
  clay_r <- subst(soil_code, from = lut$feinbod_code, to = lut$clay_mid); names(clay_r) <- "clay"
  silt_r <- subst(soil_code, from = lut$feinbod_code, to = lut$silt_mid); names(silt_r) <- "silt"
  sand_r <- subst(soil_code, from = lut$feinbod_code, to = lut$sand_mid); names(sand_r) <- "sand"
  
  clay_df <- as.data.frame(clay_r, xy = TRUE) %>% rename(value = clay) %>% filter(!is.na(value))
  silt_df <- as.data.frame(silt_r, xy = TRUE) %>% rename(value = silt) %>% filter(!is.na(value))
  sand_df <- as.data.frame(sand_r, xy = TRUE) %>% rename(value = sand) %>% filter(!is.na(value))
  
  base_theme <- theme_minimal() +
    theme(
      axis.text.x        = element_text(angle = 0, hjust = 0.5),
      plot.background    = element_rect(fill = "white", color = "white"),
      panel.background   = element_rect(fill = "white"),
      legend.background  = element_rect(fill = "white", color = "white"),
      plot.title         = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title         = element_blank(),
      axis.text          = element_text(color = "black", size = 14),
      panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major   = element_blank(),
      panel.grid.minor   = element_blank(),
      legend.position    = "bottom",
      legend.key.width   = unit(1.3, "cm"),
      legend.key.height  = unit(0.5, "cm"),
      legend.text        = element_text(size = 14),
      legend.title       = element_text(size = 14),
      strip.background   = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text         = element_text(face = "bold", size = 12)
    )
  
  p_clay <- ggplot(clay_df, aes(x = x, y = y, color = value)) +
    geom_point(size = 0.5) +
    geom_sf(data = boundary_germany, fill = NA, color = "black", inherit.aes = FALSE) +
    coord_sf(expand = FALSE) +
    scale_color_gradientn(
      colours = c("white", "lightblue", "dodgerblue", "#0072B2"),
      name = "clay (%)"
    ) + base_theme
  
  p_silt <- ggplot(silt_df, aes(x = x, y = y, color = value)) +
    geom_point(size = 0.5) +
    geom_sf(data = boundary_germany, fill = NA, color = "black", inherit.aes = FALSE) +
    coord_sf(expand = FALSE) +
    scale_color_gradientn(
      colours = c("white", "#b9f6ca", "#00bfae", "#00675b"),
      name = "silt (%)"
    ) + base_theme
  
  p_sand <- ggplot(sand_df, aes(x = x, y = y, color = value)) +
    geom_point(size = 0.5) +
    geom_sf(data = boundary_germany, fill = NA, color = "black", inherit.aes = FALSE) +
    coord_sf(expand = FALSE) +
    scale_color_gradientn(
      colours = c("white", "#ffe0b2", "orange", "#ff6600"),
      name = "sand (%)"
    ) + base_theme
  
  return(list(
    clay = p_clay, silt = p_silt, sand = p_sand,
    clay_df = clay_df, silt_df = silt_df, sand_df = sand_df
  ))
}

# 5. Generate plots and collect data
species_order <- c("Oak", "Beech", "Spruce", "Pine")
clay_plots <- list(); silt_plots <- list(); sand_plots <- list()
soil_fraction_df <- data.frame()

for (sp in species_order) {
  path <- file.path("soil_map", paste0(sp, "_soilCode_MODIS.tif"))
  p_list <- make_fraction_plots(path)
  
  clay_plots[[sp]] <- p_list$clay + labs(title = sp) +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          plot.margin = margin(t = 5, r = 5, b = 2, l = 5))
  
  silt_plots[[sp]] <- p_list$silt + theme(plot.margin = margin(t = 2, r = 5, b = 2, l = 5))
  sand_plots[[sp]] <- p_list$sand + theme(plot.margin = margin(t = 2, r = 5, b = 5, l = 5))
  
  soil_fraction_df <- bind_rows(
    soil_fraction_df,
    p_list$clay_df %>% mutate(species = sp, fraction = "clay") %>% select(species, fraction, value),
    p_list$silt_df %>% mutate(species = sp, fraction = "silt") %>% select(species, fraction, value),
    p_list$sand_df %>% mutate(species = sp, fraction = "sand") %>% select(species, fraction, value)
  )
}

# 6. Combine into composite
all_plots <- c(clay_plots[species_order], silt_plots[species_order], sand_plots[species_order])

composite <- wrap_plots(all_plots, ncol = 4, nrow = 3, guides = "collect") &
  theme(
    legend.position = "bottom",
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin     = margin(5, 5, 5, 5),
    panel.spacing   = unit(0, "cm")
  )

# 7. Display & save composite
print(composite)
ggsave("soil_map/all_species_soil_fractions_composite.png", composite, width = 12, height = 12, dpi = 300, bg = "white")

# 8. Save raw data
write.csv(soil_fraction_df, "soil_map/soil_fraction_values.csv", row.names = FALSE)

# 9. Standardized stacked bar plot
mean_fraction_df <- soil_fraction_df %>%
  group_by(species, fraction) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  group_by(species) %>%
  mutate(proportion = 100 * mean_value / sum(mean_value)) %>%
  ungroup() %>%
  mutate(
    species  = factor(species, levels = c("Oak", "Beech", "Spruce", "Pine")),
    fraction = factor(fraction, levels = c("clay", "silt", "sand"))
  )

# 10. Theme for bar plot (with Y-axis title)
bar_theme <- theme_minimal() +
  theme(
    axis.text.x        = element_text(angle = 0, hjust = 0.5, size = 14, color = "black"),
    axis.text.y        = element_text(color = "black", size = 14),
    axis.title.x       = element_blank(),
    axis.title.y       = element_text(color = "black", size = 14),
    plot.background    = element_rect(fill = "white", color = "white"),
    panel.background   = element_rect(fill = "white"),
    panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.title         = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
    legend.position    = "top",
    legend.key.width   = unit(1.3, "cm"),
    legend.key.height  = unit(0.5, "cm"),
    legend.text        = element_text(size = 14),
    legend.title       = element_text(size = 14),
    legend.background  = element_rect(fill = "white", color = "white"),
    strip.background   = element_rect(fill = "white", color = "black", linewidth = 0.5),
    strip.text         = element_text(face = "bold", size = 12)
  )

bar_plot <- ggplot(mean_fraction_df, aes(x = species, y = proportion, fill = fraction)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("clay" = "#0072B2", "silt" = "#00bfae", "sand" = "#ff6600")) +
  labs(title = "", y = "percentage (%)") +
  bar_theme

# 11. Display & save bar plot
print(bar_plot)
ggsave("soil_map/soil_texture_barplot_standardized.png", bar_plot, width = 8, height = 6, dpi = 300, bg = "white")
