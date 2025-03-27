setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")
# Load necessary libraries
library(tidyverse)

# Define species order
species_order <- c("Oak", "Beech", "Spruce", "Pine")

# Load dataset
load("results/Data/All_Species_Quantiles_PSI_TDiff.RData")

# Filter data for the four species
species_df <- all_results_df %>% 
  filter(species %in% species_order) %>%
  mutate(species = factor(species, levels = species_order))

# Bin Quantiles (1-22) and Soil Water Potential (each 50 units)
species_binned <- species_df %>%
  mutate(
    quantile_bin = cut(Quantiles, breaks = seq(0, 22, 1), include.lowest = TRUE, right = FALSE),
    swp_bin = cut(soil_water_potential, breaks = seq(floor(min(soil_water_potential)/50)*50, ceiling(max(soil_water_potential)/50)*50, 50), include.lowest = TRUE, right = FALSE)
  ) %>%
  drop_na(quantile_bin, swp_bin)

# Calculate density percentages for each bin combination within each species
density_df <- species_binned %>%
  group_by(species, quantile_bin, swp_bin) %>%
  summarize(count = n(), .groups = "drop") %>%
  group_by(species) %>%
  mutate(density_percent = (count / sum(count)) * 100) %>%
  ungroup()

# Join densities back to binned data
species_density <- species_binned %>%
  left_join(density_df, by = c("species", "quantile_bin", "swp_bin"))

# Plot the data points colored by density percentage with fitted regression lines and facets for each species
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

# Print and save the plot
print(p)
ggsave("results/density/Quantiles_PSI_linear_density.png", plot = p, width = 12, height = 8, dpi = 300)