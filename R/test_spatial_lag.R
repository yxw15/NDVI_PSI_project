setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

load("results/Data/All_Species_Quantiles_PSI_TDiff.RData")

# Load necessary libraries
library(spdep)      # For spatial neighbor and weights creation
library(spatialreg) # For spatial regression (lagsarlm)
library(ggplot2)    # For plotting
library(dplyr)      # For data manipulation

# Define species order and color palette
species_order <- c("Oak", "Beech", "Spruce", "Pine")
cb_palette <- c("Oak"    = "#E69F00",   # Orange
                "Beech"  = "#0072B2",   # Deep blue
                "Spruce" = "#009E73",   # Bluish-green
                "Pine"   = "#F0E442")   # Yellow

# Set the species factor levels in your data
all_results_df$species <- factor(all_results_df$species, levels = species_order)

# Filter data to only include the four species of interest
data_filtered <- subset(all_results_df, species %in% species_order)

# Define the results folder (no per-species subfolders)
results_folder <- "results/spatial_lag"
if(!dir.exists(results_folder)){
  dir.create(results_folder, recursive = TRUE)
}

# Initialize list to store fitted data for combined plotting
fitted_list <- list()

# Loop over each species
for(sp in species_order){
  # Subset data for the current species
  data_sp <- subset(data_filtered, species == sp)
  
  # Create coordinates matrix from x and y columns
  coords <- as.matrix(data_sp[, c("x", "y")])
  
  # Create neighbors list based on a distance threshold of 0.2 degree.
  nb <- dnearneigh(coords, 0, 0.2)
  
  # If no neighbors are found (could happen with sparse data), skip this species
  if(length(nb) == 0 || all(sapply(nb, length) == 0)){
    message(paste("No neighbors found for species", sp, "- skipping spatial analysis."))
    next
  }
  
  # Convert the neighbor list to a spatial weights list using row-standardized weights
  listw <- nb2listw(nb, style = "W")
  
  # Fit the spatial lag model with 'Quantiles' as the response and 'soil_water_potential' as predictor
  model_lag <- lagsarlm(Quantiles ~ soil_water_potential, data = data_sp, listw = listw)
  
  # Save the spatial lag model object as an RDS file in the results folder
  saveRDS(model_lag, file = file.path(results_folder, paste0("model_lag_", sp, ".rds")))
  
  # Extract fitted values and add them to the data
  data_sp$fitted <- fitted(model_lag)
  
  # Save the data with fitted values to a CSV file in the results folder
  write.csv(data_sp, file = file.path(results_folder, paste0("data_", sp, "_with_fitted.csv")), row.names = FALSE)
  
  # Store the fitted data for combined plotting
  fitted_list[[sp]] <- data_sp
}

# Combine the fitted data for all species
combined_fitted <- do.call(rbind, fitted_list) %>%
  arrange(species, soil_water_potential)

# Create a combined plot: original points and fitted lines for each species, colored by the palette
p_combined <- ggplot(combined_fitted, aes(x = soil_water_potential, y = Quantiles, color = species)) +
  geom_point(alpha = 0.6) +
  geom_line(aes(y = fitted), size = 1) +
  scale_color_manual(values = cb_palette) +
  labs(title = "Spatial Lag Model: Quantiles vs. Soil Water Potential",
       x = "Soil Water Potential",
       y = "Quantiles") +
  theme_minimal()

# Save the combined plot as a PNG file in the results folder
ggsave(filename = file.path(results_folder, "combined_fitted_line.png"),
       plot = p_combined, width = 8, height = 6)

# Optionally, display the plot
print(p_combined)
