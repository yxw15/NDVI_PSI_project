# =============================================================================
# Plot 4 species rasters (Oak, Beech, Spruce, Pine) with Germany boundary
# - Facet order: Oak, Beech, Spruce, Pine
# - Color by species using sp_cols
# - Force lon/lat axes (EPSG:4326)
# =============================================================================

rm(list = ls())

setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

library(terra)
library(tidyverse)
library(sf)
library(rnaturalearth)
library(grid)   # unit()

# ----------------------------
# Theme
# ----------------------------
base_theme <- theme_minimal() +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.text = element_text(color = "black", size = 14),
    legend.position = "bottom",
    plot.title  = element_text(hjust = 0.5, size = 18, color = "black"),
    axis.title  = element_text(size = 16),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
    axis.text.y = element_text(angle = 0, hjust = 0.5, size = 12),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey92", linewidth = 0.25),
    panel.border = element_blank(),
    strip.text = element_text(size = 14)
  )

# ----------------------------
# Species colors + facet order
# ----------------------------
sp_cols   <- c(Oak="#E69F00", Beech="#0072B2", Spruce="#009E73", Pine="#F0E442")
sp_levels <- c("Oak","Beech","Spruce","Pine")

# ----------------------------
# Germany boundary (lon/lat)
# ----------------------------
boundary_germany <- ne_countries(
  scale = "medium",
  country = "Germany",
  returnclass = "sf"
) |>
  st_transform(4326)

# ----------------------------
# Read rasters
# ----------------------------
r_oak    <- rast("species_map_MODIS/Oak.tif")
r_beech  <- rast("species_map_MODIS/Beech.tif")
r_spruce <- rast("species_map_MODIS/Spruce.tif")
r_pine   <- rast("species_map_MODIS/Pine.tif")

r_to_df <- function(r, species_name) {
  as.data.frame(r, xy = TRUE, na.rm = TRUE) |>
    rename(lon = x, lat = y, value = 3) |>
    mutate(species = species_name)
}

df_species <- bind_rows(
  r_to_df(r_oak,    "Oak"),
  r_to_df(r_beech,  "Beech"),
  r_to_df(r_spruce, "Spruce"),
  r_to_df(r_pine,   "Pine")
) |>
  mutate(
    species = factor(species, levels = sp_levels)
  )

p <- ggplot() +
  # species locations (POINTS ONLY)
  geom_point(
    data = df_species,
    aes(x = lon, y = lat, color = species),
    size = 0.5,
    alpha = 0.7
  ) +
  # Germany boundary
  geom_sf(
    data = boundary_germany,
    fill = NA,
    color = "black",
    linewidth = 0.6
  ) +
  facet_wrap(~ species, nrow = 1) +
  scale_color_manual(values = sp_cols, guide = "none") +
  coord_sf(
    crs = 4326,
    expand = FALSE
  ) +
  labs(
    x = "longitude",
    y = "latitude",
    title = ""
  ) +
  base_theme

print(p)

ggsave("results_rootzone/Figures_till2022/SI_PSI/SI_species_map.png", plot = p, width = 13, height = 5.5, dpi = 300)
