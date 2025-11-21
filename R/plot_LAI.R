setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

library(terra)
library(dplyr)
library(ggplot2)

library(rnaturalearth)
library(sf)

# Download Germany boundary as sf object
boundary_germany <- ne_countries(country = "Germany", scale = "medium", returnclass = "sf")
boundary_germany <- st_transform(boundary_germany, crs = 4326)

# 1. Find all LAI files for DOY 209 or 241, years 2003–2022
path <- "/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/LAI"
pattern <- "Lai_500m.*doy(200[3-9]|201[0-9]|202[0-2])(209|241)_.*\\.tif$"
tif_files <- list.files(
  path = path,
  pattern = pattern,
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

# 2. Compute mean LAI across files
lai.files <- rast(tif_files)
lai.files.mean <- mean(lai.files)

# 3. Read species masks
Oak    <- rast("species_map_MODIS/Oak.tif")
Beech  <- rast("species_map_MODIS/Beech.tif")
Spruce <- rast("species_map_MODIS/Spruce.tif")
Pine   <- rast("species_map_MODIS/Pine.tif")

# 4. Project mean LAI to each mask, mask to species area only
lai.oak    <- mask(project(lai.files.mean, Oak), Oak)
lai.beech  <- mask(project(lai.files.mean, Beech), Beech)
lai.spruce <- mask(project(lai.files.mean, Spruce), Spruce)
lai.pine   <- mask(project(lai.files.mean, Pine), Pine)

# 5. Convert rasters to data frames, add species label
df.oak    <- as.data.frame(lai.oak, xy=TRUE)    %>% mutate(species="Oak")
df.beech  <- as.data.frame(lai.beech, xy=TRUE)  %>% mutate(species="Beech")
df.spruce <- as.data.frame(lai.spruce, xy=TRUE) %>% mutate(species="Spruce")
df.pine   <- as.data.frame(lai.pine, xy=TRUE)   %>% mutate(species="Pine")

# 6. Combine data frames, standardize column name for LAI
lai_df <- bind_rows(df.oak, df.beech, df.spruce, df.pine)
names(lai_df)[3] <- "LAI"
lai_df$species <- factor(lai_df$species, levels = c("Oak", "Beech", "Spruce", "Pine"))

# 7. Custom plot theme
custom_theme <- theme(
  axis.text.x       = element_text(angle=0, hjust=0.5),
  plot.background   = element_rect(fill="white", color="white"),
  panel.background  = element_rect(fill="white"),
  legend.background = element_rect(fill="white", color="white"),
  plot.title        = element_text(hjust=0.5, size=18, face="bold"),
  plot.subtitle     = element_text(hjust=0.5, size=14),
  axis.title        = element_text(face="bold", size=16),
  axis.text         = element_text(color="black", size=14),
  panel.border      = element_rect(color="black", fill=NA, linewidth=0.5),
  panel.grid.major  = element_blank(),
  panel.grid.minor  = element_blank(),
  legend.position   = "top",
  legend.text       = element_text(size=14),
  strip.background  = element_rect(fill="white", color="black", linewidth=0.5),
  strip.text        = element_text(face="bold", size=12)
)

cb_palette <- c(
  "Oak"    = "#E69F00",
  "Beech"  = "#0072B2",
  "Spruce" = "#009E73",
  "Pine"   = "#F0E442"
)

plot_lai_density <- ggplot(lai_df, aes(x = LAI, fill = species, color = species)) +
  geom_density(alpha = 0.8) +
  facet_wrap(~ species, scales = "free_y") +
  labs(title = "",
       x = "LAI",
       y = "density",
       color = "",  # Remove color legend title
       fill = "") +
  scale_color_manual(values = cb_palette) +
  scale_fill_manual(values = cb_palette) +
  theme_minimal() +
  custom_theme

print(plot_lai_density)  
ggsave("results_lai/LAI_density_by_species.png", plot_lai_density, width = 10, height = 8)

# 8. Plot LAI map and bar functions

plot_lai_map <- function(df, var, title,
                     boundary = boundary_germany,
                     palette = c("white", "#e8fce8", "#bafcb6", "lightgreen", "#59d269", "#228b22", "darkgreen"),
                     limits = c(0, 5)) {
  ggplot(df) +
    geom_point(aes(x = x, y = y, color = .data[[var]]), size = 0.5) +
    geom_sf(data = boundary, fill = NA, color = "black", inherit.aes = FALSE, linewidth = 0.8) +
    scale_color_gradientn(
      colours = palette,
      values = scales::rescale(seq(limits[1], limits[2], length.out = length(palette)), from = limits),
      limits = limits,
      oob = scales::squish,
      name = "leaf area index (LAI)"
    ) +
    facet_wrap(~species, nrow = 1) +
    coord_sf(crs = 4326, expand = FALSE) +
    labs(x = "longitude", y = "latitude") +
    guides(color = guide_colorbar(title.position = "left", title.hjust = 0.5)) +
    custom_theme
}

p_map <- plot_lai_map(lai_df, "LAI", "Mean LAI by Species (DOY 209 & 241)")
ggsave("results_rootzone/Figures_till2022/supplementary/mean_LAI_map_till2022.png", p_map, width = 14, height = 6)

plot_mean_lai_bar <- function(df,
                              var = "LAI",
                              species_levels = c("Oak", "Beech", "Spruce", "Pine"),
                              title = "",
                              palette = cb_palette) {
  summary_df <- df %>%
    group_by(species) %>%
    summarize(mean_value = mean(.data[[var]], na.rm = TRUE)) %>%
    mutate(species = factor(species, levels = species_levels))
  
  ggplot(summary_df, aes(x = species, y = mean_value, fill = species)) +
    geom_col() +
    scale_fill_manual(values = palette, name = "") +
    labs(x = "", y = paste("mean leaf area index (LAI)"), title = title) +
    custom_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
}

p_bar <- plot_mean_lai_bar(lai_df)
ggsave("results_rootzone/Figures_till2022/supplementary/mean_LAI_bar_till2022.png", p_bar, width = 8, height = 6)
