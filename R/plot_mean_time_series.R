library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)

setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")
load("results/Data/All_species_month_year_Quantiles_PSI_TDiff.RData")

# Define color palette for species (Title Case)
cb_palette <- c(
  "Oak"    = "#E69F00",
  "Beech"  = "#0072B2",
  "Spruce" = "#009E73",
  "Pine"   = "#F0E442"
)

# Filter for July and August
summer_df <- final_df %>%
  filter(month %in% c("July", "August"))

# Calculate mean values by species, convert to Title Case, set factor order
summary_df <- summer_df %>%
  group_by(species) %>%
  summarise(
    # Mean_Quantiles             = mean(Quantiles, na.rm = TRUE),
    Mean_Soil_Water_Potential  = mean(soil_water_potential, na.rm = TRUE),
    Mean_Transpiration_Deficit = mean(transpiration_deficit, na.rm = TRUE)
  ) %>%
  mutate(
    species = str_to_title(species),
    species = factor(species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  )

# Reshape for plotting and relabel with case_when()
summary_long <- summary_df %>%
  pivot_longer(
    cols      = starts_with("Mean_"),
    names_to  = "Parameter",
    values_to = "Mean_Value"
  ) %>%
  mutate(
    Parameter = case_when(
      # Parameter == "Mean_Quantiles"              ~ "NDVI quantiles",
      Parameter == "Mean_Soil_Water_Potential"   ~ "Soil water potential",
      Parameter == "Mean_Transpiration_Deficit"  ~ "Transpiration deficit",
      TRUE                                       ~ Parameter
    )
  )

# Plot
p <- ggplot(summary_long, aes(x = species, y = Mean_Value, fill = species)) +
  geom_col(position = "dodge") +
  facet_wrap(~ Parameter, scales = "free_y") +
  scale_fill_manual(values = cb_palette) +
  labs(
    x    = NULL,
    y    = "mean value",
    fill = NULL
  ) +
  theme(
    axis.text.x        = element_text(angle = 0, hjust = 0.5),
    plot.background    = element_rect(fill = "white", color = "white"),
    panel.background   = element_rect(fill = "white"),
    legend.background  = element_rect(fill = "white", color = "white"),
    plot.title         = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
    plot.subtitle      = element_text(hjust = 0.5, size = 14),
    axis.title         = element_text(face = "bold", size = 16),
    axis.text          = element_text(color = "black", size = 14),
    panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "top",
    legend.text        = element_text(size = 14),
    strip.background   = element_rect(fill = "white", color = "black", linewidth = 0.5),
    strip.text         = element_text(face = "bold", size = 12)
  )

# Save
ggsave(
  filename = "results/key_displays_July_August/mean_time_series.png",
  plot     = p,
  device   = "png",
  width    = 10,
  height   = 8,
  dpi      = 300
)
