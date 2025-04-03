setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")
load("results/Data/All_Species_Quantiles_PSI_TDiff.RData")

custom_theme <- theme(
  axis.text.x = element_text(hjust = 1, size = 12),
  axis.text = element_text(color = "black", size = 12),
  plot.background = element_rect(fill = "white", color = "white"),
  panel.background = element_rect(fill = "white"),
  legend.background = element_rect(fill = "white", color = "white"),
  plot.title = element_text(hjust = 0.5, size = 20, face = "bold", color = "black"),
  axis.title = element_text(face = "bold", size = 14),
  panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  legend.position = "top",
  strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
  strip.text = element_text(face = "bold", size = 14)
)

# =====================================================
# Relationship: Oak Quantiles vs. Pine Transpiration Deficit
# =====================================================

# Summarize data by year for Oak Quantiles and Pine Transpiration Deficit
Quantiles.Oak.byYear <- all_results_df %>%
  filter(species == "Oak") %>%
  group_by(year) %>%
  summarise(mean_quantile = mean(Quantiles, na.rm = TRUE))

TranspirationDeficit.Pine.byYear <- all_results_df %>%
  filter(species == "Pine") %>%
  group_by(year) %>%
  summarise(mean_transpiration_deficit = mean(transpiration_deficit, na.rm = TRUE))

# Merge the yearly summaries
combined_df_TDiff <- inner_join(Quantiles.Oak.byYear, TranspirationDeficit.Pine.byYear, by = "year")
print(combined_df_TDiff)

# Reshape for time series plotting
combined_long_TDiff <- combined_df_TDiff %>%
  pivot_longer(cols = c(mean_quantile, mean_transpiration_deficit),
               names_to = "measure", values_to = "value")

# -------------------------------
# 1. Time Series Panels Plot
# -------------------------------
p1_TDiff <- ggplot(combined_long_TDiff, aes(x = year, y = value, group = 1)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_wrap(~ measure, scales = "free_y", ncol = 1,
             labeller = labeller(measure = c("mean_quantile" = "Panel A: Oak Quantiles",
                                             "mean_transpiration_deficit" = "Panel B: Pine Transpiration Deficit"))) +
  labs(title = "Oak Quantiles & Pine Transpiration Deficit", 
       x = "Year", y = "Value") +
  theme_minimal() +
  custom_theme
print(p1_TDiff)
ggsave("results/key_displays/Time_Series_Quantiles_oak_TDiff_pine.png", 
       plot = p1_TDiff, width = 12, height = 8, dpi = 300)

# -------------------------------
# 2. Scatter Plot with Regression (Combined Data)
# -------------------------------
correlation_TDiff <- cor(combined_df_TDiff$mean_quantile, combined_df_TDiff$mean_transpiration_deficit, 
                         use = "complete.obs")
print(paste("Correlation between Oak Quantiles and Pine Transpiration Deficit:", 
            round(correlation_TDiff, 2)))

p2_TDiff <- ggplot(combined_df_TDiff, aes(x = mean_transpiration_deficit, y = mean_quantile)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(title = paste("Correlation (r =", round(correlation_TDiff, 2), ")"),
       x = "Mean Pine Transpiration Deficit", y = "Mean Oak Quantiles") +
  theme_minimal() +
  custom_theme
print(p2_TDiff)
ggsave("results/key_displays/Scatter_Correlation_Quantiles_oak_TDiff_pine.png", 
       plot = p2_TDiff, width = 12, height = 8, dpi = 300)

# Calculate mean soil water potential by year for Oak and Pine
PSI.Oak.byYear <- all_results_df %>%
  filter(species == "Oak") %>%
  group_by(year) %>%
  summarise(mean_soil_water_potential_oak = mean(soil_water_potential, na.rm = TRUE))

PSI.Pine.byYear <- all_results_df %>%
  filter(species == "Pine") %>%
  group_by(year) %>%
  summarise(mean_soil_water_potential_Pine = mean(soil_water_potential, na.rm = TRUE))

# Merge the two summaries by year
combined_df <- inner_join(PSI.Oak.byYear, PSI.Pine.byYear, by = "year")
print(combined_df)

# Calculate correlation between Oak and Pine soil water potential
correlation_value <- cor(combined_df$mean_soil_water_potential_oak, 
                         combined_df$mean_soil_water_potential_Pine, 
                         use = "complete.obs")
print(paste("Correlation between Oak and Pine Soil Water Potential:", round(correlation_value, 2)))

# Create scatter plot with regression line using custom_theme
p <- ggplot(combined_df, aes(x = mean_soil_water_potential_oak, y = mean_soil_water_potential_Pine)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(title = paste("Correlation (r =", round(correlation_value, 2), ")"),
       x = "Mean Oak Soil Water Potential",
       y = "Mean Pine Soil Water Potential") +
  theme_minimal() +
  custom_theme
print(p)
