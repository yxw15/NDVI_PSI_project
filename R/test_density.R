library(dplyr)
library(terra)
library(tidyr)
library(viridis)

swp_bin_width <- 50
output_path <- "results/density/Quantiles_PSI_density_raster.png"

# Define species order
species_order <- c("Oak", "Beech", "Spruce", "Pine")

# Filter all_results_df for the four species and set factor levels
species_df <- all_results_df %>% 
  filter(species %in% species_order) %>%
  mutate(species = factor(species, levels = species_order))

# If soil water potential values are negative and you want a positive scale,
# uncomment the following line:
# species_df <- species_df %>% mutate(soil_water_potential = abs(soil_water_potential))

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
  rast_obj <- rast(bin_rast_norm)
  
  # Calculate the cell size (width) based on the x extent.
  cell_size <- (swp_breaks[length(swp_breaks)] - swp_breaks[1]) / n_swp
  
  # Set spatial extent:
  # - x extent is based on the soil water potential breaks
  # - y extent is adjusted so that each cell is square (cell height = cell_size)
  ext(rast_obj) <- c(swp_breaks[1], swp_breaks[length(swp_breaks)], 1, 1 + cell_size * n_quant)
  
  # Store the raster in the list with the species name
  raster_list[[sp]] <- rast_obj
}

# Create a common color palette using the inferno palette from viridis for higher contrast
col_pal <- inferno(100)

# Open a graphics device to save the multi-panel plot (PNG example)
png(filename = output_path, width = 12, height = 12, units = "in", res = 300)

# Set up a 2 x 2 plotting layout with adjusted margins
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

# Enforce a square (1:1) aspect ratio for each panel so that the cell proportions remain square
par(asp = 1)

# Loop over species and plot each raster
for (sp in species_order) {
  plot(raster_list[[sp]], 
       col = col_pal,
       main = paste("Species:", sp),
       xlab = "Soil Water Potential", 
       ylab = "Quantiles")
}

dev.off()
