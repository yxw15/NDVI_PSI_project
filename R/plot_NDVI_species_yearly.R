# ---- Packages ----
library(terra)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

# ---- Your theme ----
theme_consistent <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.background   = element_rect(fill = "white", color = "white"),
      panel.background  = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title        = element_text(hjust = 0.5, size = 18, face = "bold"),
      plot.subtitle     = element_text(hjust = 0.5, size = 14),
      axis.title        = element_text(face = "bold", size = 14),
      axis.text.y       = element_text(color = "black", size = 12),
      axis.text.x       = element_text(color = "black", size = 12),
      panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      legend.position   = "top",
      legend.text       = element_text(size = 12),
      strip.background  = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text        = element_text(face = "bold", size = 12)
    )
}

# ---- Plot function ----
# Assumes species rasters already have same CRS/resolution/extent as the NDVI quantiles.
plot_ndvi_quantiles_species_2x4 <- function(
    target_year,
    q209_path = "../WZMAllDOYs/Quantiles_209.nc",
    q241_path = "../WZMAllDOYs/Quantiles_241.nc",
    oak_path   = "species_map_MODIS/Oak.tif",
    beech_path = "species_map_MODIS/Beech.tif",
    spruce_path= "species_map_MODIS/Spruce.tif",
    pine_path  = "species_map_MODIS/Pine.tif",
    out_width  = 12, out_height = 8.5, out_dpi = 300
){
  # Load NDVI quantiles for day 209 (Jul) and day 241 (Aug)
  r209 <- rast(q209_path)
  r241 <- rast(q241_path)
  
  year_str <- as.character(target_year)
  i209 <- which(format(as.Date(time(r209)), "%Y") == year_str)
  i241 <- which(format(as.Date(time(r241)), "%Y") == year_str)
  if (length(i209) != 1) stop("Could not uniquely find year ", target_year, " in ", q209_path)
  if (length(i241) != 1) stop("Could not uniquely find year ", target_year, " in ", q241_path)
  
  q209 <- r209[[i209]]; names(q209) <- "NDVIq"
  q241 <- r241[[i241]]; names(q241) <- "NDVIq"
  
  # Species masks (same grid)
  oak    <- rast(oak_path)
  beech  <- rast(beech_path)
  spruce <- rast(spruce_path)
  pine   <- rast(pine_path)
  
  ref <- q209
  ex  <- ext(ref)
  xlim <- c(ex$xmin, ex$xmax)
  ylim <- c(ex$ymin, ex$ymax)
  
  # Germany border
  germany_sf <- rnaturalearth::ne_countries(country = "Germany", scale = "large", returnclass = "sf")
  germany_sf <- st_transform(germany_sf, crs(ref))
  
  # Mask helper: keep >0 presence
  mask_by_species <- function(ndvi, species_rast) {
    mask(ndvi, species_rast)
  }
  
  make_df <- function(ndvi, species_rast, species_name, md_label){
    m <- mask_by_species(ndvi, species_rast)
    as.data.frame(m, xy = TRUE, na.rm = TRUE) |>
      transform(species = species_name, monthday = md_label)
  }
  
  species_order <- c("Oak","Beech","Spruce","Pine")
  
  df209 <- rbind(
    make_df(q209, oak,    "Oak",    "DOY 209 (July)"),
    make_df(q209, beech,  "Beech",  "DOY 209 (July)"),
    make_df(q209, spruce, "Spruce", "DOY 209 (July)"),
    make_df(q209, pine,   "Pine",   "DOY 209 (July)")
  )
  df241 <- rbind(
    make_df(q241, oak,    "Oak",    "DOY 241 (August)"),
    make_df(q241, beech,  "Beech",  "DOY 241 (August)"),
    make_df(q241, spruce, "Spruce", "DOY 241 (August)"),
    make_df(q241, pine,   "Pine",   "DOY 241 (August)")
  )
  
  alldf <- rbind(df209, df241)
  alldf$species  <- factor(alldf$species, levels = species_order)
  alldf$monthday <- factor(alldf$monthday, levels = c("DOY 209 (July)", "DOY 241 (August)"))
  
  # Palette
  ndvi_colors <- rev(c("blue", "dodgerblue", "cyan", "yellow", "orange", "red"))
  
  # Plot
  p <- ggplot(alldf, aes(x = x, y = y)) +
    geom_point(aes(fill = NDVIq), shape = 21, size = 1.2, stroke = 0, alpha = 0.8) + # filled squares
    geom_sf(data = germany_sf, fill = NA, color = "black", linewidth = 0.25, inherit.aes = FALSE) +
    facet_grid(monthday ~ species) +
    scale_fill_gradientn(colours = ndvi_colors, name = "NDVI quantiles (rank)") +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE, crs = crs(ref)) +
    labs(x = "longitude", y = "latitude") +
    theme_consistent(base_size = 12) +
    guides(fill = guide_colorbar(barwidth = grid::unit(6, "cm"),
                                 barheight = grid::unit(0.5, "cm")))
  
  # Save to results_rootzone/Figures/<YEAR>/JulAug/NDVI_all_species_2x4_JulAug_<YEAR>.png
  out_dir <- file.path("results_rootzone", "Figures", target_year, "JulAug")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  outfile <- file.path(out_dir, sprintf("NDVI_all_species_2x4_JulAug_%d.png", target_year))
  
  ggsave(outfile, p, width = out_width, height = out_height, dpi = out_dpi)
  message("Saved: ", outfile)
  invisible(outfile)
}

# ---- Run for your years ----
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")
years_to_make <- c(2003, 2005, 2018, 2022)
invisible(lapply(years_to_make, plot_ndvi_quantiles_species_2x4))
