plot_soil_fraction_composite <- function(
    species_order = c("Oak", "Beech", "Spruce", "Pine"),
    soil_map_dir = "soil_map",
    result_data_dir = "results_rootzone/Data",
    result_fig_dir = "results_rootzone/Figures",
    lut_path = file.path(soil_map_dir, "feinbod_lookup_with_english.csv"),
    texture_path = file.path(soil_map_dir, "soil_texture_classes.csv"),
    output_plot = file.path(result_fig_dir, "all_species_soil_fractions_composite.png"),
    output_csv  = file.path(result_data_dir, "soil_fraction_values.csv"),
    plot_width  = 12,
    plot_height = 12,
    plot_dpi    = 300
) {
  # Load required libraries
  library(terra)
  library(dplyr)
  library(ggplot2)
  library(sf)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(patchwork)
  
  # Read Germany boundary
  boundary_germany <- ne_countries(
    scale = "medium",
    country = "Germany",
    returnclass = "sf"
  )
  
  # Read lookup tables and build LUT
  soil_descrip <- read.csv(lut_path, stringsAsFactors = FALSE)
  soil_texture <- read.csv(texture_path, stringsAsFactors = FALSE)
  lut <- soil_descrip %>%
    select(feinbod_code, feinbod) %>%
    inner_join(soil_texture, by = c("feinbod" = "code")) %>%
    mutate(
      clay_mid = (clay_lower + clay_upper) / 2,
      silt_mid = (silt_lower + silt_upper) / 2,
      sand_mid = (sand_lower + sand_upper) / 2
    ) %>%
    select(feinbod_code, clay_mid, silt_mid, sand_mid)
  
  # Function to generate fraction plots and data for one raster
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
        axis.text.x      = element_text(angle = 0, hjust = 0.5),
        plot.background  = element_rect(fill = "white", color = "white"),
        panel.background = element_rect(fill = "white"),
        legend.background= element_rect(fill = "white", color = "white"),
        plot.title       = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
        axis.title       = element_blank(),
        axis.text        = element_text(color = "black", size = 14),
        panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position  = "bottom",
        legend.key.width = unit(1.3, "cm"),
        legend.key.height= unit(0.5, "cm"),
        legend.text      = element_text(size = 14),
        legend.title     = element_text(size = 14),
        strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
        strip.text       = element_text(face = "bold", size = 12)
      )
    
    p_clay <- ggplot(clay_df, aes(x = x, y = y, color = value)) +
      geom_point(size = 0.5) +
      geom_sf(data = boundary_germany, fill = NA, color = "black", inherit.aes = FALSE) +
      coord_sf(expand = FALSE) +
      scale_color_gradientn(colours = c("white", "lightblue", "dodgerblue", "#0072B2"), name = "clay (%)") +
      base_theme
    
    p_silt <- ggplot(silt_df, aes(x = x, y = y, color = value)) +
      geom_point(size = 0.5) +
      geom_sf(data = boundary_germany, fill = NA, color = "black", inherit.aes = FALSE) +
      coord_sf(expand = FALSE) +
      scale_color_gradientn(colours = c("white", "#b9f6ca", "#00bfae", "#00675b"), name = "silt (%)") +
      base_theme
    
    p_sand <- ggplot(sand_df, aes(x = x, y = y, color = value)) +
      geom_point(size = 0.5) +
      geom_sf(data = boundary_germany, fill = NA, color = "black", inherit.aes = FALSE) +
      coord_sf(expand = FALSE) +
      scale_color_gradientn(colours = c("white", "#ffe0b2", "orange", "#ff6600"), name = "sand (%)") +
      base_theme
    
    list(
      clay = p_clay,
      silt = p_silt,
      sand = p_sand,
      data = bind_rows(
        clay_df %>% mutate(fraction = "clay"),
        silt_df %>% mutate(fraction = "silt"),
        sand_df %>% mutate(fraction = "sand")
      )
    )
  }
  
  # Prepare storage for plots and data
  plot_list <- list()
  soil_fraction_df <- tibble()
  
  # Generate plots per species
  for (sp in species_order) {
    raster_file <- file.path(soil_map_dir, paste0(sp, "_soilCode_MODIS.tif"))
    res <- make_fraction_plots(raster_file)
    plot_list[[sp]] <- res[c("clay", "silt", "sand")]
    soil_fraction_df <- bind_rows(soil_fraction_df, res$data %>% mutate(species = sp))
  }
  
  # Define desired row (fraction) order
  fraction_order <- c("silt", "clay", "sand")
  
  # Build composite list: row-by-row
  all_plots_list <- list()
  for (fr in fraction_order) {
    for (sp in species_order) {
      p <- plot_list[[sp]][[fr]] +
        labs(title = sp) +
        theme(plot.margin = margin(5, 5, 2, 5))
      all_plots_list <- c(all_plots_list, list(p))
    }
  }
  
  # Assemble and save composite
  composite <- wrap_plots(all_plots_list,
                          ncol = length(species_order),
                          nrow = length(fraction_order),
                          guides = "collect") &
    theme(
      legend.position = "bottom",
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin     = margin(5, 5, 5, 5),
      panel.spacing   = unit(0, "cm")
    )
  
  ggsave(output_plot, composite,
         width = plot_width,
         height = plot_height,
         dpi = plot_dpi,
         bg = "white")
  
  # Save raw data
  write.csv(soil_fraction_df, output_csv, row.names = FALSE)
  
  message("Composite plot saved to: ", output_plot)
  message("Soil fraction data saved to: ", output_csv)
}


plot_soil_fraction_composite()
