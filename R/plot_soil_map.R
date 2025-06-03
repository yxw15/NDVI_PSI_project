setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

library(ncdf4)
soil_map <- nc_open("../ALLAN_PIA_SoilMoisture/BUEK1000N/BUEK1000Nv232_BOFORM_GK3_1km.nc")
# crs(soil_map) <- "epsg:31467"
nc_close(soil_map)

# adjust the path if needed
soil_info <- read.csv(
  file        = "../ALLAN_PIA_SoilMoisture/BUEK1000N/bodenprofile_bk1000n_bfviewer.csv",
  header      = TRUE,              # first row contains column names
  stringsAsFactors = FALSE,        # keep character columns as strings
  na.strings  = c("", "NA")        # treat empty fields as NA
)

# 1. Extract all unique, non‐NA feinbod values in sorted order:
fein_levels <- sort(
  unique( soil_info$feinbod[ !is.na(soil_info$feinbod) ] )
)

# 2. For each feinbod code, pull the first matching description from soil_info:
#    (If a given feinbod appears more than once, we take the first non‐NA description.)
fein_descrs <- sapply(fein_levels, function(x) {
  # all descriptions where feinbod == x
  des <- unique( soil_info$feinbod_descr[ soil_info$feinbod == x ] )
  # drop any NA, then take the first
  des <- des[ !is.na(des) ]
  if (length(des) == 0) return( NA_character_ )
  return(des[1])
}, USE.NAMES = FALSE)

# 3. Build a data.frame with the integer code, feinbod, and its description:
fein_lookup <- data.frame(
  feinbod_code   = seq_along(fein_levels),
  feinbod        = fein_levels,
  feinbod_descr  = fein_descrs,
  stringsAsFactors = FALSE
)

# 4. Write out as CSV (no row‐names):
write.csv(
  fein_lookup,
  file      = "soil_map/feinbod_lookup.csv",
  row.names = FALSE,
  quote     = TRUE
)

fein_descr_en <- c(
  # 1
  "Hh (humic topsoil)",                   
  # 2
  "Hhe–Hhs (transitional humic top-/subsoil)",  
  # 3
  "Hhs (humic subsoil)",                        
  # 4
  "Hn (humic peat soil)",                       
  # 5
  "Ls2 (slightly sandy loam)",                   
  # 6
  "Ls3 (moderately sandy loam)",                 
  # 7
  "Ls4 (strongly sandy loam)",                   
  # 8
  "Lt2 (slightly clayey loam)",                  
  # 9
  "Lt3 (moderately clayey loam)",                
  # 10
  "Lts (sandy-clayey loam)",                     
  # 11
  "Lu (silty loam)",                             
  # 12
  "Sl2 (slightly loamy sand)",                   
  # 13
  "Sl3 (moderately loamy sand)",                 
  # 14
  "Sl4 (strongly loamy sand)",                   
  # 15
  "Slu (silty-loamy sand)",                      
  # 16
  "Ss (pure sand)",                              
  # 17
  "St2 (slightly clayey sand)",                  
  # 18
  "St3 (moderately clayey sand)",                
  # 19
  "Su2 (slightly silty sand)",                   
  # 20
  "Su3 (moderately silty sand)",                 
  # 21
  "Su4 (strongly silty sand)",                    
  # 22
  "Tl (loamy clay)",                             
  # 23
  "Ts2 (slightly sandy clay)",                    
  # 24
  "Ts3 (moderately sandy clay)",                  
  # 25
  "Ts4 (strongly sandy clay)",                    
  # 26
  "Tt (pure clay)",                              
  # 27
  "Tu2 (slightly silty clay)",                    
  # 28
  "Tu3 (moderately silty clay)",                  
  # 29
  "Tu4 (strongly silty clay)",                    
  # 30
  "Uls (sandy-loamy silt)",                       
  # 31
  "Us (sandy silt)",                              
  # 32
  "Ut2 (slightly clayey silt)",                    
  # 33
  "Ut3 (moderately clayey silt)",                  
  # 34
  "Ut4 (strongly clayey silt)",                    
  # 35
  "Uu (pure silt)"                                
)

# (Make sure length(fein_descr_en) == nrow(fein_lookup).)

# 3. Add this as a new column to your data.frame:
fein_lookup$feinbod_descr_en <- fein_descr_en

# 4. Write out the new CSV with four columns:
#      feinbod_code, feinbod, feinbod_descr, feinbod_descr_en
write.csv(
  fein_lookup,
  file      = "soil_map/feinbod_lookup_with_english.csv",
  row.names = FALSE,
  quote     = TRUE
)

soil_info2 <- merge(
  x      = soil_info,
  y      = fein_lookup,
  by     = "feinbod",
  all.x  = TRUE
)


library(terra)
soil_map <- rast("soil_map/boform_wald.tif")

lut_all <- unique(
  soil_info2[, c("boform_id", "feinbod_code")]
)

lut_mat <- cbind(
  as.numeric(lut_all$boform_id),
  as.numeric(lut_all$feinbod_code)
)

soil_reclass <- classify(soil_map, lut_mat)
writeRaster(soil_reclass, 
            filename = "soil_map/soil_code.tif", 
            overwrite= TRUE)

vals_df <- terra::unique(soil_reclass)    # returns a data.frame with one column
classes  <- vals_df[[1]]                  # pull out the numeric codes column
# If there might be NA in vals_df, drop it:
classes  <- classes[!is.na(classes)]

# (C) For each unique class code, create a one‐layer raster where only pixels == class stay;
#     everything else becomes NA.  Then stack them.
layer_list <- lapply(classes, function(cl) {
  lyr <- soil_reclass
  lyr[soil_reclass != cl] <- NA
  return(lyr)
})

# (D) Combine into a multi‐layer SpatRaster
soil_reclass_layers <- rast(layer_list)

# (E) Optionally, name each band by its class code:
fein_names <- fein_lookup$feinbod[match(classes, fein_lookup$feinbod_code)]
names(soil_reclass_layers) <- fein_names

# (F) Write out the 34‐band (or however many classes there actually were) TIFF:
writeRaster(
  soil_reclass_layers,
  "soil_map/soil_reclass_layers.tif",
  overwrite = TRUE
)

library(ggplot2)
library(dplyr)

# -------------------------------------------------------------------
# 1. Convert the single‐band SpatRaster to a data.frame
# -------------------------------------------------------------------

# Suppose `soil_reclass` is your SpatRaster with values 1:34.
# We use terra::as.data.frame() with xy=TRUE to get columns "x", "y", and "boform_wald".
df_rast <- as.data.frame(soil_reclass, xy = TRUE, na.rm = FALSE)

# By default, terra will name the value column after your raster’s varname (here "boform_wald").
# Rename it to "feinbod_code" so it matches your lookup table column:
df_rast <- df_rast %>%
  rename(feinbod_code = boform_wald)

# You’ll now have a data.frame with columns:
#   x, y, feinbod_code
# where `feinbod_code` is an integer 1:34 (or NA).

# -------------------------------------------------------------------
# 2. Join with the lookup table to get the “feinbod” string for each code
# -------------------------------------------------------------------

# We only need “feinbod_code” ↔ “feinbod” in this join.  Pixels where feinbod_code is NA will remain NA.
df_joined <- df_rast %>%
  left_join(
    fein_lookup %>% select(feinbod_code, feinbod_descr_en),
    by = "feinbod_code"
  )

# At this point, `df_joined` has columns:
#   x, y, feinbod_code, feinbod
# For any pixel where feinbod_code was NA (or didn’t match), `feinbod` will be NA.

# If you want to drop every pixel whose code is NA (i.e. mask background), do:
df_plot <- df_joined %>% filter(!is.na(feinbod_descr_en))

# -------------------------------------------------------------------
# 3. Plot with ggplot2, mapping fill to the factor “feinbod”
# -------------------------------------------------------------------
soil_code_map <- ggplot(df_plot, aes(x = x, y = y, fill = feinbod_descr_en)) +
  geom_raster(na.rm = TRUE) +
  scale_fill_viridis_d(
    name = "soil type",     
    guide = guide_legend(
      ncol = 1,             
      keyheight = unit(0.8, "lines")
    )
  ) +
  coord_equal() +
  labs(
    x = "easting (m)",
    y = "northing (m)",
    title = "DHDN / 3-degree Gauss-Kruger zone 3 (EPSG:31467)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x       = element_text(angle = 0, hjust = 0.5, size = 14, color = "black"),
    axis.text.y       = element_text(size = 14, color = "black"),
    axis.title        = element_text(face = "bold", size = 16),
    plot.title        = element_text(hjust = 0.5, size = 18, face = "bold"),
    plot.background   = element_rect(fill = "white", color = "white"),
    panel.background  = element_rect(fill = "white"),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    # legend.position   = "top",
    legend.title = element_text(face = "bold",size = 14),
    legend.background = element_rect(fill = "white", color = "white"),
    legend.text       = element_text(size = 14),
    strip.background  = element_rect(fill = "white", color = "black", linewidth = 0.5),
    strip.text        = element_text(face = "bold", size = 12)
  )
output_path_soil_code_map <- "soil_map/soil_code_map.png"
ggsave(filename = output_path_soil_code_map, plot = soil_code_map, width = 9, height = 6, dpi = 300)


# inspect structure and a few rows
str(soil)
head(soil)
summary(soil)

unique(soil$feinbod)
unique(soil$feinbod_descr)

library(terra)
nc_file <- "../ALLAN_PIA_SoilMoisture/BUEK1000N/BUEK1000Nv232_BOFORM_GK3_1km.nc"
wald_rast <- rast(nc_file, subds = "boform_wald")
wald_rast
plot(wald_rast, main = "Soil Profile ID of Forest (boform_wald)")
writeRaster(wald_rast, "soil_map/boform_wald.tif", overwrite = TRUE)

f <- freq(wald_rast)               # a data.frame with columns 'value' and 'count'
n_unique <- nrow(f)                # one row per distinct value
cat("Number of unique boform_wald values:", n_unique, "\n")

species.map <- rast("../TreeSpeciesGermany/species_class_sum_INT1U.tif")
# Beech
message("⏳ Creating Beech layer…")
Beech <- ifel(species.map == 2, species.map, NA)
message("   • Writing ‘Beech.tif’ to disk…")
writeRaster(
  Beech, 
  filename  = "species_map_original/Beech.tif", 
  overwrite = TRUE
)
message("✅ Beech.tif written!\n")

# Oak
message("⏳ Creating Oak layer…")
Oak <- ifel(species.map == 5, species.map, NA)
message("   • Writing ‘Oak.tif’ to disk…")
writeRaster(
  Oak, 
  filename  = "species_map_original/Oak.tif", 
  overwrite = TRUE
)
message("✅ Oak.tif written!\n")

# Spruce
message("⏳ Creating Spruce layer…")
Spruce <- ifel(species.map == 8, species.map, NA)
message("   • Writing ‘Spruce.tif’ to disk…")
writeRaster(
  Spruce, 
  filename  = "species_map_original/Spruce.tif", 
  overwrite = TRUE
)
message("✅ Spruce.tif written!\n")

# Pine
message("⏳ Creating Pine layer…")
Pine <- ifel(species.map == 9, species.map, NA)
message("   • Writing ‘Pine.tif’ to disk…")
writeRaster(
  Pine, 
  filename  = "species_map_original/Pine.tif", 
  overwrite = TRUE
)
message("✅ Pine.tif written!\n")


### Plot the soil percentage of each species ###
soil_code <- rast("soil_map/soil_code.tif")
soil_code_layer <- rast("soil_map/soil_reclass_layers.tif")
species.MODIS <- rast("../TreeSpeciesGermany/TreeSpeciesMODIS.tif")
# Beech.orig <- rast("species_map_original/Beech.tif")
# Beech.prj.soil <- project(Beech.orig, soil_code)
Beech.MODIS <- species.MODIS[[2]]
Beech.MODIS[ Beech.MODIS > 0.7 ] <- 3
Beech.MODIS[ Beech.MODIS <= 0.7 ] <- NA

library(terra)

# Step 1: Create a mask where Beech is present
soil_masked <- mask(soil_code_layer, Beech.prj.soil)

# Step 3: Count total beech pixels (in the mask)
total_beech_pixels <- global(beech_mask, "sum", na.rm = TRUE)[1, 1]

# Step 4: For each soil class layer, count how many beech pixels have that class
results <- data.frame(
  layer = names(soil_code_layer),
  percent_of_beech = numeric(nlyr(soil_code_layer))
)

for (i in seq_along(results$layer)) {
  # Count how many Beech pixels are classified as this soil class
  class_count <- global(!is.na(soil_masked[[i]]), "sum", na.rm = TRUE)[1, 1]
  
  # Percentage of beech area this layer covers
  results$percent_of_beech[i] <- (class_count / total_beech_pixels) * 100
}

# Optional sanity check
results <- results %>%
  mutate(percent_of_beech = round(percent_of_beech, 2))

# View result
print(results)

soil_descrip <- read.csv("soil_map/feinbod_lookup_with_english.csv")
library(dplyr)

# Join the English descriptions
results_joined <- left_join(results, soil_descrip, by = c("layer" = "feinbod"))

# Keep only relevant columns
results_final <- results_joined %>%
  select(layer, feinbod_descr_en, percent_of_beech) %>%
  arrange(desc(percent_of_beech))

# View final result
print(results_final)

library(ggplot2)

cb_palette <- c(
  "Oak"    = "#E69F00",
  "Beech"  = "#0072B2",
  "Spruce" = "#009E73",
  "Pine"   = "#F0E442"
)

# Convert percent_of_beech to numeric (if it's character)
results_final$percent_of_beech <- as.numeric(results_final$percent_of_beech)

# Plot
ggplot(results_final, aes(x = reorder(feinbod_descr_en, percent_of_beech), 
                          y = percent_of_beech)) +
  geom_bar(stat = "identity", fill = "forestgreen") +
  coord_flip() +
  labs(
    x = "soil type",
    y = "percent of beech pixels (%)",
    title = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.x       = element_text(angle = 0, hjust = 0.5),
    plot.background   = element_rect(fill = "white", color = "white"),
    panel.background  = element_rect(fill = "white"),
    legend.background = element_rect(fill = "white", color = "white"),
    plot.title        = element_text(hjust = 0.5, size = 18, face = "bold"),
    plot.subtitle     = element_text(hjust = 0.5, size = 14),
    axis.title        = element_text(face = "bold", size = 16),
    axis.text         = element_text(color = "black", size = 14),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    legend.position   = "top",
    legend.text       = element_text(size = 14),
    strip.background  = element_rect(fill = "white", color = "black", linewidth = 0.5),
    strip.text        = element_text(face = "bold", size = 12)
  )


### Percentage of four species ###
library(terra)
library(dplyr)
library(ggplot2)
library(readr)

# Load soil data
soil_code_layer <- rast("soil_map/soil_reclass_layers.tif")
soil_descrip <- read.csv("soil_map/feinbod_lookup_with_english.csv")

# Load individual species rasters
species_rasters <- list(
  Oak    = rast("species_map_original/Oak.tif"),
  Beech  = rast("species_map_original/Beech.tif"),
  Spruce = rast("species_map_original/Spruce.tif"),
  Pine   = rast("species_map_original/Pine.tif")
)

# Project species rasters to match soil
species_rasters <- lapply(species_rasters, function(r) project(r, soil_code_layer))

# Container for results
all_species_results <- list()

# Loop through each species
for (sp in names(species_rasters)) {
  sp_raster <- species_rasters[[sp]]
  sp_mask <- ifel(!is.na(sp_raster), 1, NA)
  masked_soil <- mask(soil_code_layer, sp_mask)
  
  total_pixels <- global(!is.na(sp_raster), "sum", na.rm = TRUE)[1, 1]
  
  df <- data.frame(
    layer = names(soil_code_layer),
    percent = numeric(nlyr(soil_code_layer))
  )
  
  for (i in seq_along(df$layer)) {
    valid_pixels <- global(!is.na(masked_soil[[i]]), "sum", na.rm = TRUE)[1, 1]
    df$percent[i] <- (valid_pixels / total_pixels) * 100
  }
  
  df$species <- sp
  all_species_results[[sp]] <- df
}

# Combine and join descriptions
results_all <- bind_rows(all_species_results) %>%
  left_join(soil_descrip, by = c("layer" = "feinbod")) %>%
  mutate(percent = round(percent, 2)) %>%
  select(species, feinbod_descr_en, percent)

# Color palette
species_order <- c("Oak", "Beech", "Spruce", "Pine")
cb_palette <- c("Oak"   = "#E69F00",  # Orange
                "Beech" = "#0072B2",  # Deep blue
                "Spruce"= "#009E73",  # Bluish-green
                "Pine"  = "#F0E442")  # Yellow
results_all$species <- factor(results_all$species, levels = species_order)

# Plot
soil_statistics <- ggplot(results_all, aes(x = reorder(feinbod_descr_en, percent), y = percent, fill = species)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  scale_fill_manual(values = cb_palette) +
  labs(
    x = "soil type",
    y = "percent of species pixels (%)",
    fill = "",
    title = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.x       = element_text(angle = 0, hjust = 0.5),
    plot.background   = element_rect(fill = "white", color = "white"),
    panel.background  = element_rect(fill = "white"),
    legend.background = element_rect(fill = "white", color = "white"),
    plot.title        = element_text(hjust = 0.5, size = 18, face = "bold"),
    plot.subtitle     = element_text(hjust = 0.5, size = 14),
    axis.title        = element_text(face = "bold", size = 16),
    axis.text         = element_text(color = "black", size = 14),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    legend.position   = "top",
    legend.text       = element_text(size = 14),
    strip.background  = element_rect(fill = "white", color = "black", linewidth = 0.5),
    strip.text        = element_text(face = "bold", size = 12)
  )

print(soil_statistics)
output_path_soil_statistics <- "soil_map/soil_statistics.png"
ggsave(filename = output_path_soil_statistics, plot = soil_statistics, width = 9, height = 6, dpi = 300)


# Loop over species
for (sp in unique(results_all$species)) {
  species_df <- results_all %>% filter(species == sp)
  
  p <- ggplot(species_df, aes(x = reorder(feinbod_descr_en, percent), y = percent)) +
    geom_bar(stat = "identity", fill = cb_palette[sp]) +
    coord_flip() +
    labs(
      x = "soil type",
      y = paste0("percent of pixels (%)"),
      title = ""
    ) +
    theme_minimal() +
    theme(
      axis.text.x       = element_text(angle = 0, hjust = 0.5),
      plot.background   = element_rect(fill = "white", color = "white"),
      panel.background  = element_rect(fill = "white"),
      legend.position   = "none",
      plot.title        = element_text(hjust = 0.5, size = 18, face = "bold"),
      axis.title        = element_text(face = "bold", size = 16),
      axis.text         = element_text(color = "black", size = 14),
      panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank()
    )
  
  print(p)
  # Save the plot
  ggsave(
    filename = paste0("soil_map/soil_statistics_", sp, ".png"),
    plot = p,
    width = 9,
    height = 6,
    dpi = 300
  )
}
