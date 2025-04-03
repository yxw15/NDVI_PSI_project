# =====================================================
# Set Working Directory, Load Data, and Required Libraries
# =====================================================
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")
load("results/Data/All_Species_Quantiles_PSI_TDiff.RData")
library(dplyr)
library(ggplot2)
library(tidyr)

# A sample of the original data:
# > head(all_results_df)
# # A tibble: 6 × 7
#       x     y year  Quantiles soil_water_potential transpiration_deficit species
#   <dbl> <dbl> <chr>     <dbl>                <dbl>                 <dbl> <chr>  
# 1  9.62  54.8 2003          9               -33.7                 0.159  Beech  
# 2  9.62  54.8 2004         16               -11.8                 0.0879 Beech  
# 3  9.62  54.8 2005          7               -15.7                 0.0642 Beech  
# 4  9.62  54.8 2006          3               -71.5                 0.138  Beech  
# 5  9.62  54.8 2007         18                -6.97                0.126  Beech  
# 6  9.62  54.8 2008         21               -61.7                 0.0602 Beech

# =====================================================
# Define a Custom Theme for All Plots
# =====================================================
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
# Section 1: Oak Quantiles vs. Pine Soil Water Potential
# =====================================================

# --- 1a. Yearly Summaries & Time Series Panels ---
Quantiles.Oak.byYear <- all_results_df %>%
  filter(species == "Oak") %>%
  group_by(year) %>%
  summarise(mean_quantile = mean(Quantiles, na.rm = TRUE))

PSI.Pine.byYear <- all_results_df %>%
  filter(species == "Pine") %>%
  group_by(year) %>%
  summarise(mean_soil_water_potential = mean(soil_water_potential, na.rm = TRUE))

combined_df_SWP <- inner_join(Quantiles.Oak.byYear, PSI.Pine.byYear, by = "year")
print(combined_df_SWP)

combined_long_SWP <- combined_df_SWP %>%
  pivot_longer(cols = c(mean_quantile, mean_soil_water_potential),
               names_to = "measure", values_to = "value")

p1_SWP <- ggplot(combined_long_SWP, aes(x = year, y = value, group = 1)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_wrap(~ measure, scales = "free_y", ncol = 1,
             labeller = labeller(measure = c("mean_quantile" = "Panel A: Oak Quantiles",
                                             "mean_soil_water_potential" = "Panel B: Pine Soil Water Potential"))) +
  labs(title = "Time Series Panels: Oak Quantiles & Pine Soil Water Potential", 
       x = "Year", y = "Value") +
  theme_minimal() +
  custom_theme
print(p1_SWP)
ggsave("results/key_displays/Time_Series_Quantiles_oak_PSI_pine.png", 
       plot = p1_SWP, width = 12, height = 8, dpi = 300)

# --- 1b. Scatter Plot with Regression (Yearly Data) ---
correlation_SWP <- cor(combined_df_SWP$mean_quantile, combined_df_SWP$mean_soil_water_potential, use = "complete.obs")
print(paste("Correlation between Oak Quantiles and Pine Soil Water Potential:", round(correlation_SWP, 2)))

p2_SWP <- ggplot(combined_df_SWP, aes(x = mean_quantile, y = mean_soil_water_potential)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(title = paste("Correlation (r =", round(correlation_SWP, 2), ")"),
       x = "Mean Oak Quantiles", y = "Mean Pine Soil Water Potential") +
  theme_minimal() +
  custom_theme
print(p2_SWP)
ggsave("results/key_displays/Scatter_Correlation_Quantiles_oak_PSI_pine.png", 
       plot = p2_SWP, width = 12, height = 8, dpi = 300)

# --- 1c. Regression Plot Using Random Sampling ---
set.seed(123)
oak_sample_SWP <- all_results_df %>%
  filter(species == "Oak", !is.na(Quantiles)) %>%
  sample_n(600000)

pine_sample_SWP <- all_results_df %>%
  filter(species == "Pine", !is.na(soil_water_potential)) %>%
  sample_n(600000)

reg_df_SWP <- data.frame(Quantiles.Oak = oak_sample_SWP$Quantiles,
                         PSI.Pine = pine_sample_SWP$soil_water_potential)

model_SWP <- lm(Quantiles.Oak ~ PSI.Pine, data = reg_df_SWP)
summary(model_SWP)

p3_SWP <- ggplot(reg_df_SWP, aes(x = PSI.Pine, y = Quantiles.Oak)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(title = "Linear Regression: Oak Quantiles vs. Pine Soil Water Potential",
       x = "Soil Water Potential for Pine", y = "Oak Quantiles") +
  theme_minimal() +
  custom_theme
print(p3_SWP)
ggsave("results/key_displays/Regression_Quantiles_oak_PSI_pine.png", 
       plot = p3_SWP, width = 12, height = 8, dpi = 300)

# --- 1d. Binned Mean Points Plot for Soil Water Potential ---
bin_width_SWP <- 50
total_pixels_SWP <- nrow(reg_df_SWP)
binned_points_SWP <- reg_df_SWP %>%
  mutate(swp_bin = cut(PSI.Pine,
                       breaks = seq(floor(min(PSI.Pine) / bin_width_SWP) * bin_width_SWP,
                                    ceiling(max(PSI.Pine) / bin_width_SWP) * bin_width_SWP,
                                    bin_width_SWP),
                       include.lowest = TRUE,
                       right = FALSE)) %>%
  group_by(swp_bin) %>%
  summarize(count = n(),
            mean_swp = mean(PSI.Pine, na.rm = TRUE),
            mean_quantile = mean(Quantiles.Oak, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(count >= total_pixels_SWP * 0.001)  # Remove bins with <0.1% of total pixels

p4_SWP <- ggplot(binned_points_SWP, aes(x = mean_swp, y = mean_quantile)) +
  geom_point(size = 3, color = "red") +
  labs(title = "Binned Mean: Oak Quantiles vs. Pine Soil Water Potential",
       x = "Mean Soil Water Potential (50 unit bins)",
       y = "Mean Oak Quantiles") +
  theme_minimal() +
  custom_theme
print(p4_SWP)
ggsave("results/key_displays/Binned_Mean_Quantiles_oak_PSI_pine.png", 
       plot = p4_SWP, width = 12, height = 8, dpi = 300)

# =====================================================
# Section 2: Oak Quantiles vs. Pine Transpiration Deficit
# =====================================================

# --- 2a. Yearly Summaries & Time Series Panels ---
Quantiles.Oak.byYear <- all_results_df %>%
  filter(species == "Oak") %>%
  group_by(year) %>%
  summarise(mean_quantile = mean(Quantiles, na.rm = TRUE))

TranspirationDeficit.Pine.byYear <- all_results_df %>%
  filter(species == "Pine") %>%
  group_by(year) %>%
  summarise(mean_transpiration_deficit = mean(transpiration_deficit, na.rm = TRUE))

combined_df_TDiff <- inner_join(Quantiles.Oak.byYear, TranspirationDeficit.Pine.byYear, by = "year")
print(combined_df_TDiff)

combined_long_TDiff <- combined_df_TDiff %>%
  pivot_longer(cols = c(mean_quantile, mean_transpiration_deficit),
               names_to = "measure", values_to = "value")

p1_TDiff <- ggplot(combined_long_TDiff, aes(x = year, y = value, group = 1)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_wrap(~ measure, scales = "free_y", ncol = 1,
             labeller = labeller(measure = c("mean_quantile" = "Panel A: Oak Quantiles",
                                             "mean_transpiration_deficit" = "Panel B: Pine Transpiration Deficit"))) +
  labs(title = "Time Series Panels: Oak Quantiles & Pine Transpiration Deficit", 
       x = "Year", y = "Value") +
  theme_minimal() +
  custom_theme
print(p1_TDiff)
ggsave("results/key_displays/Time_Series_Quantiles_oak_TDiff_pine.png", 
       plot = p1_TDiff, width = 12, height = 8, dpi = 300)

# --- 2b. Scatter Plot with Regression (Yearly Data) ---
correlation_TDiff <- cor(combined_df_TDiff$mean_quantile, combined_df_TDiff$mean_transpiration_deficit, 
                         use = "complete.obs")
print(paste("Correlation between Oak Quantiles and Pine Transpiration Deficit:", 
            round(correlation_TDiff, 2)))

p2_TDiff <- ggplot(combined_df_TDiff, aes(x = mean_quantile, y = mean_transpiration_deficit)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(title = paste("Correlation (r =", round(correlation_TDiff, 2), ")"),
       x = "Mean Oak Quantiles", y = "Mean Pine Transpiration Deficit") +
  theme_minimal() +
  custom_theme
print(p2_TDiff)
ggsave("results/key_displays/Scatter_Correlation_Quantiles_oak_TDiff_pine.png", 
       plot = p2_TDiff, width = 12, height = 8, dpi = 300)

# --- 2c. Regression Plot Using Random Sampling ---
set.seed(123)
oak_sample_TDiff <- all_results_df %>%
  filter(species == "Oak", !is.na(Quantiles)) %>%
  sample_n(600000)

pine_sample_TDiff <- all_results_df %>%
  filter(species == "Pine", !is.na(transpiration_deficit)) %>%
  sample_n(600000)

reg_df_TDiff <- data.frame(Quantiles.Oak = oak_sample_TDiff$Quantiles,
                           TDiff.Pine = pine_sample_TDiff$transpiration_deficit)

model_TDiff <- lm(Quantiles.Oak ~ TDiff.Pine, data = reg_df_TDiff)
summary(model_TDiff)

p3_TDiff <- ggplot(reg_df_TDiff, aes(x = TDiff.Pine, y = Quantiles.Oak)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(title = "Linear Regression: Oak Quantiles vs. Pine Transpiration Deficit",
       x = "Pine Transpiration Deficit", y = "Oak Quantiles") +
  theme_minimal() +
  custom_theme
print(p3_TDiff)
ggsave("results/key_displays/Regression_Quantiles_oak_TDiff_pine.png", 
       plot = p3_TDiff, width = 12, height = 8, dpi = 300)

# --- 2d. Density Plot with Regression for Transpiration Deficit ---
plot_density_reg_TDiff <- function(data, tdiff_bin_width = 3, output_path) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  binned_df <- data %>%
    mutate(
      quantile_bin = cut(Quantiles.Oak,
                         breaks = seq(0, 22, 1),
                         include.lowest = TRUE,
                         right = FALSE),
      tdiff_bin = cut(TDiff.Pine,
                      breaks = seq(floor(min(TDiff.Pine) / tdiff_bin_width) * tdiff_bin_width,
                                   ceiling(max(TDiff.Pine) / tdiff_bin_width) * tdiff_bin_width,
                                   tdiff_bin_width),
                      include.lowest = TRUE,
                      right = FALSE)
    ) %>%
    drop_na(quantile_bin, tdiff_bin)
  
  density_df <- binned_df %>%
    group_by(quantile_bin, tdiff_bin) %>%
    summarize(count = n(), .groups = "drop") %>%
    group_by(tdiff_bin) %>%
    mutate(density_percent = (count / sum(count)) * 100) %>%
    ungroup()
  
  density_data <- left_join(binned_df, density_df, by = c("quantile_bin", "tdiff_bin"))
  
  max_density <- ceiling(max(density_data$density_percent, na.rm = TRUE))
  density_breaks <- seq(0, max_density, length.out = 5)
  
  p <- ggplot(density_data, aes(x = TDiff.Pine, y = Quantiles.Oak, color = density_percent)) +
    geom_point(size = 2.5, alpha = 0.9) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1) +
    scale_color_viridis_c(
      name = "Density (%)",
      option = "inferno",
      direction = -1,
      trans = "sqrt",
      breaks = density_breaks
    ) +
    labs(
      title = "Quantiles (Oak) vs. Transpiration Deficit (Pine) with Density",
      x = "Pine Transpiration Deficit",
      y = "Oak Quantiles"
    ) +
    theme_minimal() +
    custom_theme
  
  print(p)
  ggsave(output_path, plot = p, width = 12, height = 8, dpi = 300)
}

plot_density_reg_TDiff(reg_df_TDiff, tdiff_bin_width = 3, 
                       output_path = "results/key_displays/Quantiles_oak_TDiff_pine_density.png")

# --- 2e. Binned Mean Points Plot for Transpiration Deficit ---
bin_width_TDiff <- 3
total_pixels_TDiff <- nrow(reg_df_TDiff)
binned_points_TDiff <- reg_df_TDiff %>%
  mutate(tdiff_bin = cut(TDiff.Pine,
                         breaks = seq(floor(min(TDiff.Pine) / bin_width_TDiff) * bin_width_TDiff,
                                      ceiling(max(TDiff.Pine) / bin_width_TDiff) * bin_width_TDiff,
                                      bin_width_TDiff),
                         include.lowest = TRUE,
                         right = FALSE)) %>%
  group_by(tdiff_bin) %>%
  summarize(count = n(),
            mean_tdiff = mean(TDiff.Pine, na.rm = TRUE),
            mean_quantile = mean(Quantiles.Oak, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(count >= total_pixels_TDiff * 0.0001)  # Remove bins with <0.01% of total pixels

p4_TDiff <- ggplot(binned_points_TDiff, aes(x = mean_tdiff, y = mean_quantile)) +
  geom_point(size = 3, color = "red") +
  labs(title = "Binned Mean: Oak Quantiles vs. Pine Transpiration Deficit",
       x = "Mean Pine Transpiration Deficit (3 unit bins)",
       y = "Mean Oak Quantiles") +
  theme_minimal() +
  custom_theme
print(p4_TDiff)
ggsave("results/key_displays/Binned_Mean_Quantiles_oak_TDiff_pine.png", 
       plot = p4_TDiff, width = 12, height = 8, dpi = 300)
