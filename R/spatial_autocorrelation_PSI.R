# ===============================================================
# Set Working Directory & Load Data
# ===============================================================
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

# Load the combined data (for further analysis and ACF plotting)
load("results/Data/All_Species_Quantiles_PSI_TDiff.RData")

# ===============================================================
# Load Required Libraries
# ===============================================================
library(sp)
library(dplyr)
library(spdep)
library(ggplot2)

# ===============================================================
# Filter for Beech species in 2018 and Sample Data
# ===============================================================
data <- all_results_df %>% filter(species == "Beech" & year == "2018")

# Create a coordinates matrix from the x and y columns of the sample
coords <- cbind(data$x, data$y)
PSI <- data$soil_water_potential

# ===============================================================
# Calculate Moran's I for Multiple Distance Thresholds
# ===============================================================
# Define distance thresholds every 0.2 degree
dist_seq <- seq(0.6, 1.0, by = 0.2)

# Initialize a list to store individual data frames for each distance
df_list <- list()

# Loop over distance thresholds
for (i in seq_along(dist_seq)) {
  d <- dist_seq[i]
  
  # Define neighbors within distance 'd' (from 0 up to d)
  nb <- dnearneigh(coords, 0, d)
  
  # Create spatial weights (row-standardized, handling isolated points)
  lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
  
  # Check if any neighbors exist; if none, assign NA values
  if (sum(unlist(nb)) == 0) {
    temp_df <- data.frame(distance = d, Moran_I = NA, p_value = NA)
  } else {
    # Run the Moran's I test
    test <- moran.test(PSI, lw, zero.policy = TRUE)
    
    # Extract Moran's I statistic and p-value
    moran_I <- test$estimate[1]
    p_val <- test$p.value
    
    # Save the results for this distance as a one-row data frame
    temp_df <- data.frame(distance = d, Moran_I = moran_I, p_value = p_val)
  }
  
  # Save the individual result as a CSV file with the naming format "PSI_spatial_[distance].csv"
  file_name <- sprintf("PSI_spatial_%.1f.csv", d)
  write.csv(temp_df, file = file_name, row.names = FALSE)
  
  # Store the result in the list for combining later
  df_list[[i]] <- temp_df
}

# Combine all individual data frames into one
results_spatial <- do.call(rbind, df_list)

# ===============================================================
# Plot Moran's I vs. Distance
# ===============================================================
p <- ggplot(results_spatial, aes(x = distance, y = Moran_I)) +
  geom_line() +
  geom_point() +
  labs(x = "Distance Threshold",
       y = "Moran's I",
       title = "Moran's I vs. Distance for PSI") +
  theme_minimal()

# Display the plot
print(p)

# ===============================================================
# Save the Combined Results and Plot
# ===============================================================
# Save the combined results data frame as a CSV file named "PSI_spatial_distance.csv"
write.csv(results_spatial, file = "results_spatial/PSI_spatial_distance.csv", row.names = FALSE)

# Save the plot as a PNG file
ggsave("results_spatial/results_PSI_spatial_plot.png", plot = p, width = 8, height = 6)
