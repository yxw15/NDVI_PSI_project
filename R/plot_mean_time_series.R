setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")
load("results/Data/All_species_month_year_Quantiles_PSI_TDiff.RData")

library(dplyr)
library(ggplot2)
library(tidyr)

# Define color palette for species
cb_palette <- c(
  "oak"   = "#E69F00",  # Orange
  "beech" = "#0072B2",  # Deep blue
  "spruce"= "#009E73",  # Bluish-green
  "pine"  = "#F0E442"   # Yellow
)

# Filter for July and August
summer_df <- final_df %>%
  filter(month %in% c("July", "August"))

# Calculate mean values by species
summary_df <- summer_df %>%
  group_by(species) %>%
  summarise(
    Mean_Quantiles = mean(Quantiles, na.rm = TRUE),
    Mean_Soil_Water_Potential = mean(soil_water_potential, na.rm = TRUE),
    Mean_Transpiration_Deficit = mean(transpiration_deficit, na.rm = TRUE)
  ) %>%
  mutate(species = factor(tolower(species), levels = c("oak", "beech", "spruce", "pine")))

# Reshape for plotting
summary_long <- summary_df %>%
  pivot_longer(
    cols = starts_with("Mean_"),
    names_to = "Parameter",
    values_to = "Mean_Value"
  ) %>%
  mutate(
    Parameter = case_when(
      Parameter == "Mean_Quantiles" ~ "NDVI quantiles",
      Parameter == "Mean_Soil_Water_Potential" ~ "soil water potential",
      Parameter == "Mean_Transpiration_Deficit" ~ "transpiration deficit",
      TRUE ~ Parameter
    )
  )

# Plot
p <- ggplot(summary_long, aes(x = species, y = Mean_Value, fill = species)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Parameter, scales = "free_y") +
  scale_fill_manual(values = cb_palette) +
  labs(
    # title = "",
    subtitle = "",
    x = "",
    y = "mean value",
    fill = ""
  ) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.background = element_rect(fill = "white", color = "white"),
    panel.background = element_rect(fill = "white"),
    legend.background = element_rect(fill = "white", color = "white"),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
    plot.subtitle = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_text(color = "black", size = 14),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 14),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    strip.text = element_text(face = "bold", size = 12)
  )

ggsave(filename = "results/key_displays_July_August/mean_time_series.png", plot = p, device = "png", width = 10, height = 8, dpi = 300)
