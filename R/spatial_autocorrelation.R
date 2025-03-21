# ===============================================================
# Set Working Directory & Load Data
# ===============================================================
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

# Load the combined data (for further analysis and ACF plotting)
load("results/Data/All_Species_Quantiles_PSI_TDiff.RData")

# ===============================================================
# Load Required Libraries
# ===============================================================
library(dplyr)
library(spdep)
library(ggplot2)

# ===============================================================
# Data Preparation: Filter for Beech and Remove Missing Values
# ===============================================================
Beech_df <- all_results_df %>% 
  filter(species == "Beech") %>% 
  filter(!is.na(x), !is.na(y), !is.na(year), !is.na(Quantiles))

# Compute an upper bound for maximum distance using the bounding box diagonal
coords_all <- cbind(Beech_df$x, Beech_df$y)
max_x <- max(coords_all[, 1])
min_x <- min(coords_all[, 1])
max_y <- max(coords_all[, 2])
min_y <- min(coords_all[, 2])
max_dist <- sqrt((max_x - min_x)^2 + (max_y - min_y)^2)

# Create distance bins (0.1 degree increments)
bins <- seq(0.1, max_dist, by = 0.1)

# ===============================================================
# Moran's I Calculation by Year and Distance Bin
# ===============================================================
# Initialize an empty data frame to store results
results <- data.frame(year = character(),
                      distance = numeric(),
                      moransI = numeric(),
                      stringsAsFactors = FALSE)

# Loop over each year in the Beech dataset
for (yr in unique(Beech_df$year)) {
  yr_data <- Beech_df %>% filter(year == yr)
  
  # Skip years with too few points to compute spatial autocorrelation
  if(nrow(yr_data) < 3) next
  
  coords <- cbind(yr_data$x, yr_data$y)
  
  # Loop over each distance bin
  for (bin in bins) {
    lower <- bin - 0.1
    upper <- bin
    
    # Identify neighbors within the current distance interval (lower, upper]
    nb <- dnearneigh(coords, lower, upper)
    
    # If no neighbors are found for any point, record NA for this bin
    if (sum(sapply(nb, length)) == 0) {
      results <- rbind(results, data.frame(year = yr, distance = bin, moransI = NA))
    } else {
      # Create spatial weights list with zero.policy = TRUE to handle empty neighbor sets
      lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
      # Compute Moran's I using zero.policy = TRUE in the test as well
      mi <- moran.test(yr_data$Quantiles, lw, zero.policy = TRUE)$estimate["Moran I statistic"]
      results <- rbind(results, data.frame(year = yr, distance = bin, moransI = mi))
    }
  }
}

# Remove NA values from the results
results_clean <- results %>% filter(!is.na(moransI))

# For each distance bin, calculate the mean Moran's I over all years
avg_results <- results_clean %>% 
  group_by(distance) %>% 
  summarise(mean_moransI = mean(moransI, na.rm = TRUE), .groups = "drop")

# ===============================================================
# Plotting the Results
# ===============================================================
# Plot Moran's I vs. distance for each year (faceted by year)
p1 <- ggplot(results_clean, aes(x = distance, y = moransI)) +
  geom_line() +
  facet_wrap(~ year, scales = "free_y") +
  labs(x = "Distance (degrees)",
       y = "Moran's I",
       title = "Spatial Autocorrelation (Quantiles) for Beech by Year") +
  theme_minimal()

# Plot the mean Moran's I vs. distance averaged over all years
p2 <- ggplot(avg_results, aes(x = distance, y = mean_moransI)) +
  geom_line() +
  labs(x = "Distance (degrees)",
       y = "Mean Moran's I",
       title = "Mean Spatial Autocorrelation (Quantiles) for Beech Averaged Over Years") +
  theme_minimal()

# Display the plots
print(p1)
print(p2)

# ===============================================================
# Save the Results: Data Frames and Figures
# ===============================================================
# Save the results in a list object
results_spatial <- list(yearly = results_clean, averaged = avg_results)

# Save the list as an RData file
save(results_spatial, file = "results_spatial/Beech_results_spatial.RData")

# Also, write the data frames to CSV files
write.csv(results_clean, "results_spatial/Beech_yearly_moransI_results.csv", row.names = FALSE)
write.csv(avg_results, "results_spatial/Beech_averaged_moransI_results.csv", row.names = FALSE)

# Save the figures as PNG files (adjust width, height, and dpi as needed)
ggsave("results_spatial/Beech_moran_by_year.png", plot = p1, width = 8, height = 6, dpi = 300)
ggsave("results_spatial/Beech_mean_moran.png", plot = p2, width = 8, height = 6, dpi = 300)
