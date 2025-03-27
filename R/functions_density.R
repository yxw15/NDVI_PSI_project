plot_density_Quantiles_PSI_linear <- function(data, swp_bin_width = 50, output_path) {
  
  # Define species order
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  
  # Filter data for the four species and set factor levels
  species_df <- data %>% 
    filter(species %in% species_order) %>%
    mutate(species = factor(species, levels = species_order))
  
  # Bin Quantiles (1-22) and Soil Water Potential (using swp_bin_width)
  species_binned <- species_df %>%
    mutate(
      quantile_bin = cut(Quantiles, breaks = seq(0, 22, 1), include.lowest = TRUE, right = FALSE),
      swp_bin = cut(soil_water_potential, 
                    breaks = seq(floor(min(soil_water_potential) / swp_bin_width) * swp_bin_width,
                                 ceiling(max(soil_water_potential) / swp_bin_width) * swp_bin_width,
                                 swp_bin_width),
                    include.lowest = TRUE, right = FALSE)
    ) %>%
    drop_na(quantile_bin, swp_bin)
  
  # Calculate density percentages for each bin combination, normalizing by each swp_bin (column)
  density_df <- species_binned %>%
    group_by(species, quantile_bin, swp_bin) %>%
    summarize(count = n(), .groups = "drop") %>%
    group_by(species, swp_bin) %>%  # Group by species and soil water potential bin
    mutate(density_percent = (count / sum(count)) * 100) %>%
    ungroup()
  
  # Join densities back to binned data
  species_density <- species_binned %>%
    left_join(density_df, by = c("species", "quantile_bin", "swp_bin"))
  
  # Create the plot with fitted linear regression lines and facets for each species
  p <- ggplot(species_density, aes(x = soil_water_potential, y = Quantiles, color = density_percent)) +
    geom_point(size = 2.5, alpha = 0.9) +
    geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 1) +
    scale_color_viridis_c(name = "Density (%)", option = "plasma", direction = -1,
                          breaks = seq(0, ceiling(max(species_density$density_percent, na.rm = TRUE)), by = 1)) +
    labs(
      title = "Quantiles vs. Soil Water Potential with Density by Species",
      x = "Soil Water Potential",
      y = "Quantiles"
    ) +
    facet_wrap(~species, ncol = 2) +
    theme_minimal()
  
  # Print and save the plot to the specified output path
  print(p)
  ggsave(output_path, plot = p, width = 12, height = 8, dpi = 300)
}

plot_density_Quantiles_TDiff_linear <- function(data, tdiff_bin_width = 3, output_path) {
  
  # Define species order
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  
  # Filter data for the four species and set factor levels
  species_df <- data %>% 
    filter(species %in% species_order) %>%
    mutate(species = factor(species, levels = species_order))
  
  # Bin Quantiles (1-22) and Transpiration Deficit (using tdiff_bin_width)
  species_binned <- species_df %>%
    mutate(
      quantile_bin = cut(Quantiles, breaks = seq(0, 22, 1), include.lowest = TRUE, right = FALSE),
      tdiff_bin = cut(transpiration_deficit, 
                      breaks = seq(floor(min(transpiration_deficit) / tdiff_bin_width) * tdiff_bin_width,
                                   ceiling(max(transpiration_deficit) / tdiff_bin_width) * tdiff_bin_width,
                                   tdiff_bin_width),
                      include.lowest = TRUE, right = FALSE)
    ) %>%
    drop_na(quantile_bin, tdiff_bin)
  
  # Calculate density percentages for each bin combination,
  # normalizing by each transpiration deficit bin (tdiff_bin) per species
  density_df <- species_binned %>%
    group_by(species, quantile_bin, tdiff_bin) %>%
    summarize(count = n(), .groups = "drop") %>%
    group_by(species, tdiff_bin) %>%  # Group by species and tdiff_bin
    mutate(density_percent = (count / sum(count)) * 100) %>%
    ungroup()
  
  # Join densities back to binned data
  species_density <- species_binned %>%
    left_join(density_df, by = c("species", "quantile_bin", "tdiff_bin"))
  
  # Create the plot with fitted linear regression lines and facets for each species
  p <- ggplot(species_density, aes(x = transpiration_deficit, y = Quantiles, color = density_percent)) +
    geom_point(size = 2.5, alpha = 0.9) +
    geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 1) +
    scale_color_viridis_c(name = "Density (%)", option = "plasma", direction = -1,
                          breaks = seq(0, ceiling(max(species_density$density_percent, na.rm = TRUE)), by = 1)) +
    labs(
      title = "Quantiles vs. Transpiration Deficit with Density by Species",
      x = "Transpiration Deficit",
      y = "Quantiles"
    ) +
    facet_wrap(~species, ncol = 2) +
    theme_minimal()
  
  # Print and save the plot to the specified output path
  print(p)
  ggsave(output_path, plot = p, width = 12, height = 8, dpi = 300)
}

plot_density_TDiff_PSI_linear <- function(data, swp_bin_width = 50, tdiff_bin_width = 3, output_path) {
  
  # Define species order
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  
  # Filter data for the four species and set factor levels
  species_df <- data %>% 
    filter(species %in% species_order) %>%
    mutate(species = factor(species, levels = species_order))
  
  # Bin Soil Water Potential and Transpiration Deficit using provided bin widths
  species_binned <- species_df %>%
    mutate(
      swp_bin = cut(soil_water_potential, 
                    breaks = seq(floor(min(soil_water_potential) / swp_bin_width) * swp_bin_width,
                                 ceiling(max(soil_water_potential) / swp_bin_width) * swp_bin_width,
                                 swp_bin_width),
                    include.lowest = TRUE, right = FALSE),
      tdiff_bin = cut(transpiration_deficit, 
                      breaks = seq(floor(min(transpiration_deficit) / tdiff_bin_width) * tdiff_bin_width,
                                   ceiling(max(transpiration_deficit) / tdiff_bin_width) * tdiff_bin_width,
                                   tdiff_bin_width),
                      include.lowest = TRUE, right = FALSE)
    ) %>%
    drop_na(swp_bin, tdiff_bin)
  
  # Calculate density percentages for each bin combination,
  # normalizing by each soil water potential bin (swp_bin) within each species
  density_df <- species_binned %>%
    group_by(species, swp_bin, tdiff_bin) %>%
    summarize(count = n(), .groups = "drop") %>%
    group_by(species, swp_bin) %>%  # Group by species and each soil water potential bin (column)
    mutate(density_percent = (count / sum(count)) * 100) %>%
    ungroup()
  
  # Join densities back to binned data
  species_density <- species_binned %>%
    left_join(density_df, by = c("species", "swp_bin", "tdiff_bin"))
  
  # Create the plot with fitted regression lines and facets for each species
  p <- ggplot(species_density, aes(x = soil_water_potential, y = transpiration_deficit, color = density_percent)) +
    geom_point(size = 2.5, alpha = 0.9) +
    geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 1) +
    scale_color_viridis_c(name = "Density (%)", option = "plasma", direction = -1,
                          breaks = seq(0, ceiling(max(species_density$density_percent, na.rm = TRUE)), by = 1)) +
    labs(
      title = "Transpiration Deficit vs. Soil Water Potential with Density by Species",
      x = "Soil Water Potential",
      y = "Transpiration Deficit"
    ) +
    facet_wrap(~species, ncol = 2) +
    theme_minimal()
  
  # Print and save the plot to the specified output path
  print(p)
  ggsave(output_path, plot = p, width = 12, height = 8, dpi = 300)
}

plot_density_Quantiles_PSI_raster <- function(data, swp_bin_width = 50, output_path) {
  # Load required packages
  library(dplyr)
  library(terra)
  
  # Define species order
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  
  # Filter data for the four species and set factor levels
  species_df <- data %>% 
    filter(species %in% species_order) %>%
    mutate(species = factor(species, levels = species_order))
  
  # Bin Quantiles (1-22) and Soil Water Potential using swp_bin_width
  species_binned <- species_df %>%
    mutate(
      quantile_bin = cut(Quantiles, 
                         breaks = seq(0, 22, 1), 
                         include.lowest = TRUE, 
                         right = FALSE),
      swp_bin = cut(soil_water_potential, 
                    breaks = seq(floor(min(soil_water_potential) / swp_bin_width) * swp_bin_width,
                                 ceiling(max(soil_water_potential) / swp_bin_width) * swp_bin_width,
                                 swp_bin_width),
                    include.lowest = TRUE, 
                    right = FALSE)
    ) %>%
    drop_na(quantile_bin, swp_bin)
  
  # Get the factor levels (bins) for quantiles and soil water potential
  quant_levels <- levels(species_binned$quantile_bin)  # should be 22 bins
  swp_levels   <- levels(species_binned$swp_bin)         # number of PSI classes based on swp_bin_width
  n_quant <- length(quant_levels)
  n_swp   <- length(swp_levels)
  
  # Compute the soil water potential breaks (used for binning) so we can set the x extent later.
  swp_breaks <- seq(floor(min(species_df$soil_water_potential) / swp_bin_width) * swp_bin_width,
                    ceiling(max(species_df$soil_water_potential) / swp_bin_width) * swp_bin_width,
                    swp_bin_width)
  
  # Create a list to store the raster for each species
  raster_list <- list()
  
  # Loop over each species
  for (sp in species_order) {
    # Subset the binned data for the current species
    sp_data <- species_binned %>% filter(species == sp)
    
    # Initialize an empty matrix for counts: rows = quantile bins (1:22), columns = PSI classes (n_swp)
    bin_rast <- matrix(0, nrow = n_quant, ncol = n_swp)
    
    # Loop over each PSI class (soil water potential bin) and each quantile bin.
    # We reverse the row order so that lower quantile values appear at the bottom.
    for (i in 1:n_swp) {
      for (j in 1:n_quant) {
        count_val <- sp_data %>% 
          filter(swp_bin == swp_levels[i], quantile_bin == quant_levels[j]) %>%
          nrow()
        bin_rast[n_quant + 1 - j, i] <- count_val
      }
      print(paste("Species", sp, "- swp bin", i))
    }
    
    # Normalize each column (i.e. each PSI class) so that the column sums to 1.
    bin_rast_norm <- apply(bin_rast, 2, function(x) {
      if (sum(x) == 0) rep(0, length(x)) else x / sum(x)
    })
    
    # Convert the normalized matrix to a raster.
    # The number of rows corresponds to the quantile bins and the columns to the swp bins.
    rast_obj <- rast(bin_rast_norm)
    
    # Set spatial extent: x extent is based on the soil water potential breaks,
    # y extent is from 1 to n_quant.
    ext(rast_obj) <- c(swp_breaks[1], swp_breaks[length(swp_breaks)], 1, n_quant)
    
    # Store the raster in the list with the species name
    raster_list[[sp]] <- rast_obj
  }
  
  # Create a common color palette
  col_pal <- colorRampPalette(c("grey", "yellow", "red"))
  
  # Open a graphics device to save the multi-panel plot (PNG example)
  png(filename = output_path, width = 12, height = 8, units = "in", res = 300)
  
  # Set up a 2 x 2 plotting layout
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  
  # Loop over species and plot each raster
  for (sp in species_order) {
    plot(raster_list[[sp]], 
         col = col_pal(100),
         main = paste("Species:", sp),
         xlab = "Soil Water Potential", 
         ylab = "Quantiles")
  }
  
  dev.off()
}

plot_density_Quantiles_TDiff_raster <- function(data, tdiff_bin_width = 3, output_path) {
  # Load required packages
  library(dplyr)
  library(terra)
  
  # Define species order
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  
  # Filter data for the four species and set factor levels
  species_df <- data %>% 
    filter(species %in% species_order) %>%
    mutate(species = factor(species, levels = species_order))
  
  # Bin Quantiles (0-22) and Transpiration Deficit (using tdiff_bin_width)
  species_binned <- species_df %>%
    mutate(
      quantile_bin = cut(Quantiles, 
                         breaks = seq(0, 22, 1), 
                         include.lowest = TRUE, 
                         right = FALSE),
      tdiff_bin = cut(transpiration_deficit, 
                      breaks = seq(floor(min(transpiration_deficit) / tdiff_bin_width) * tdiff_bin_width,
                                   ceiling(max(transpiration_deficit) / tdiff_bin_width) * tdiff_bin_width,
                                   tdiff_bin_width),
                      include.lowest = TRUE, 
                      right = FALSE)
    ) %>%
    drop_na(quantile_bin, tdiff_bin)
  
  # Determine the bin levels and counts
  quant_levels <- levels(species_binned$quantile_bin)  # expected to be 22 bins
  tdiff_levels <- levels(species_binned$tdiff_bin)        # number of transpiration deficit bins
  n_quant <- length(quant_levels)
  n_tdiff <- length(tdiff_levels)
  
  # Compute the transpiration deficit breaks (for setting the x extent)
  tdiff_breaks <- seq(floor(min(species_df$transpiration_deficit) / tdiff_bin_width) * tdiff_bin_width,
                      ceiling(max(species_df$transpiration_deficit) / tdiff_bin_width) * tdiff_bin_width,
                      tdiff_bin_width)
  
  # Create a list to store the raster for each species
  raster_list <- list()
  
  # Loop over each species
  for (sp in species_order) {
    # Subset the binned data for the current species
    sp_data <- species_binned %>% filter(species == sp)
    
    # Initialize an empty matrix for counts: rows = quantile bins, columns = transpiration deficit bins
    bin_rast <- matrix(0, nrow = n_quant, ncol = n_tdiff)
    
    # Loop over each transpiration deficit bin (columns) and each quantile bin (rows)
    # Reverse row order so that lower Quantiles appear at the bottom.
    for (i in 1:n_tdiff) {
      for (j in 1:n_quant) {
        count_val <- sp_data %>% 
          filter(tdiff_bin == tdiff_levels[i], quantile_bin == quant_levels[j]) %>%
          nrow()
        bin_rast[n_quant + 1 - j, i] <- count_val
      }
      print(paste("Species", sp, "- tdiff bin", i))
    }
    
    # Normalize each column (each transpiration deficit bin) so the column sums to 1
    bin_rast_norm <- apply(bin_rast, 2, function(x) {
      if (sum(x) == 0) rep(0, length(x)) else x / sum(x)
    })
    
    # Convert the normalized matrix to a raster using terra::rast
    rast_obj <- rast(bin_rast_norm)
    
    # Set spatial extent:
    # - x-axis from the first to the last transpiration deficit break,
    # - y-axis from 1 to the number of quantile bins.
    ext(rast_obj) <- c(tdiff_breaks[1], tdiff_breaks[length(tdiff_breaks)], 1, n_quant)
    
    # Store the raster in the list with the species name
    raster_list[[sp]] <- rast_obj
  }
  
  # Create a common color palette (grey → yellow → red)
  col_pal <- colorRampPalette(c("grey", "yellow", "red"))
  
  # Open a graphics device to save the multi-panel plot (PNG example)
  png(filename = output_path, width = 12, height = 8, units = "in", res = 300)
  
  # Set up a 2 x 2 plotting layout (one panel per species)
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  
  # Loop over species and plot each raster
  for (sp in species_order) {
    plot(raster_list[[sp]], 
         col = col_pal(100),
         main = paste("Species:", sp),
         xlab = "Transpiration Deficit", 
         ylab = "Quantiles")
  }
  
  dev.off()
}

plot_density_TDiff_PSI_raster <- function(data, swp_bin_width = 50, tdiff_bin_width = 3, output_path) {
  # Load required packages
  library(dplyr)
  library(terra)
  
  # Define species order
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  
  # Filter data for the four species and set factor levels
  species_df <- data %>% 
    filter(species %in% species_order) %>%
    mutate(species = factor(species, levels = species_order))
  
  # Bin Soil Water Potential and Transpiration Deficit using provided bin widths
  species_binned <- species_df %>%
    mutate(
      swp_bin = cut(soil_water_potential, 
                    breaks = seq(floor(min(soil_water_potential) / swp_bin_width) * swp_bin_width,
                                 ceiling(max(soil_water_potential) / swp_bin_width) * swp_bin_width,
                                 swp_bin_width),
                    include.lowest = TRUE, right = FALSE),
      tdiff_bin = cut(transpiration_deficit, 
                      breaks = seq(floor(min(transpiration_deficit) / tdiff_bin_width) * tdiff_bin_width,
                                   ceiling(max(transpiration_deficit) / tdiff_bin_width) * tdiff_bin_width,
                                   tdiff_bin_width),
                      include.lowest = TRUE, right = FALSE)
    ) %>%
    drop_na(swp_bin, tdiff_bin)
  
  # Get the factor levels for swp and transpiration deficit bins
  swp_levels <- levels(species_binned$swp_bin)   # PSI classes based on soil water potential
  tdiff_levels <- levels(species_binned$tdiff_bin) # transpiration deficit bins
  n_swp <- length(swp_levels)
  n_tdiff <- length(tdiff_levels)
  
  # Compute the breaks for swp and tdiff to define raster extents
  swp_breaks <- seq(floor(min(species_df$soil_water_potential) / swp_bin_width) * swp_bin_width,
                    ceiling(max(species_df$soil_water_potential) / swp_bin_width) * swp_bin_width,
                    swp_bin_width)
  tdiff_breaks <- seq(floor(min(species_df$transpiration_deficit) / tdiff_bin_width) * tdiff_bin_width,
                      ceiling(max(species_df$transpiration_deficit) / tdiff_bin_width) * tdiff_bin_width,
                      tdiff_bin_width)
  
  # Create a list to store the raster for each species
  raster_list <- list()
  
  # Loop over each species to build the count matrix, normalize, and convert to a raster
  for (sp in species_order) {
    # Subset the binned data for the current species
    sp_data <- species_binned %>% filter(species == sp)
    
    # Initialize an empty matrix with rows = tdiff bins and columns = swp bins
    bin_rast <- matrix(0, nrow = n_tdiff, ncol = n_swp)
    
    # Loop over each soil water potential bin (columns) and transpiration deficit bin (rows)
    # Reverse the row order so that lower tdiff values appear at the bottom.
    for (i in 1:n_swp) {
      for (j in 1:n_tdiff) {
        count_val <- sp_data %>% 
          filter(swp_bin == swp_levels[i], tdiff_bin == tdiff_levels[j]) %>%
          nrow()
        bin_rast[n_tdiff + 1 - j, i] <- count_val
      }
      print(paste("Species", sp, "- swp bin", i))
    }
    
    # Normalize each column so that each soil water potential bin sums to 1.
    bin_rast_norm <- apply(bin_rast, 2, function(x) {
      if (sum(x) == 0) rep(0, length(x)) else x / sum(x)
    })
    
    # Convert the normalized matrix to a raster using terra::rast
    rast_obj <- rast(bin_rast_norm)
    
    # Set spatial extent:
    #   x-axis: from the first to the last soil water potential break,
    #   y-axis: from the first to the last transpiration deficit break.
    ext(rast_obj) <- c(swp_breaks[1], swp_breaks[length(swp_breaks)],
                       tdiff_breaks[1], tdiff_breaks[length(tdiff_breaks)])
    
    # Store the raster with the species name as key
    raster_list[[sp]] <- rast_obj
  }
  
  # Create a common color palette: grey → yellow → red
  col_pal <- colorRampPalette(c("grey", "yellow", "red"))
  
  # Open a graphics device (here, a PNG file) to save the multi-panel plot
  png(filename = output_path, width = 12, height = 8, units = "in", res = 300)
  
  # Set up a 2 x 2 plotting layout (one panel per species)
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  
  # Loop over each species and plot its raster
  for (sp in species_order) {
    plot(raster_list[[sp]], 
         col = col_pal(100),
         main = paste("Species:", sp),
         xlab = "Soil Water Potential", 
         ylab = "Transpiration Deficit")
  }
  
  dev.off()
}

setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")
load("results/Data/All_Species_Quantiles_PSI_TDiff.RData")

plot_density_Quantiles_PSI_linear(data = all_results_df,
                                  output_path = "results/density/Quantiles_PSI_density_linear.png")

plot_density_Quantiles_TDiff_linear(data = all_results_df,
                                    output_path = "results/density/Quantiles_TDiff_density_linear.png")

plot_density_TDiff_PSI_linear(data = all_results_df,
                              output_path = "results/density/TDiff_PSI_density_linear.png")

plot_density_Quantiles_PSI_raster(data = all_results_df,
                                  output_path = "results/density/Quantiles_PSI_density_raster.png")

plot_density_Quantiles_TDiff_raster(data = all_results_df,
                                    output_path = "results/density/Quantiles_TDiff_density_raster.png")

plot_density_TDiff_PSI_raster(data = all_results_df,
                              output_path = "results/density/TDiff_PSI_density_raster.png")

