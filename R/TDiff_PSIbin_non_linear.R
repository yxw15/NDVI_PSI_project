# Load required libraries
library(lme4)      # For mixed-effects modeling
library(dplyr)     # For data manipulation
library(ggplot2)   # For plotting
library(gridExtra) # For arranging multiple plots
library(patchwork) # For combining ggplots

setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")
load("results/Data/NDVI_PSI_TDiff_species_df.RData")

# Process the data and remove any rows with missing values
TDiff_PSIbin_df <- TDiff_PSIbin(NDVI_PSI_TDiff_species)
TDiff_PSIbin_df <- na.omit(TDiff_PSIbin_df)

# Define a custom color palette and name the colors by species.
cb_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442")
names(cb_palette) <- c("Oak", "Beech", "Spruce", "Pine")

# Set species as a factor in the desired order.
species_levels <- c("Oak", "Beech", "Spruce", "Pine")
TDiff_PSIbin_df$species <- factor(TDiff_PSIbin_df$species, levels = species_levels)

## Panel A: Mixed-Effects Model Plot

# Fit the mixed-effects model with a third-order polynomial for bin_median.
model <- lmer(avg_transpiration_deficit ~ poly(bin_median, 3) + 
                (poly(bin_median, 3) | species),
              data = TDiff_PSIbin_df)

# Create prediction data for each species.
pred_data <- TDiff_PSIbin_df %>%
  group_by(species) %>%
  summarise(min_bin = min(bin_median),
            max_bin = max(bin_median)) %>%
  group_by(species) %>%
  do(data.frame(bin_median = seq(.$min_bin, .$max_bin, length.out = 100))) %>%
  ungroup()

# Ensure the species factor is in the same order.
pred_data$species <- factor(pred_data$species, levels = species_levels)

# Predict the fitted values using the model (including random effects).
pred_data$predicted <- predict(model, newdata = pred_data, re.form = NULL)

# Calculate 30% of the global maximum average transpiration deficit (for visual reference)
max_deficit <- max(TDiff_PSIbin_df$avg_transpiration_deficit, na.rm = TRUE)
line_val <- 0.3 * max_deficit

# Create the mixed-effects model plot with the horizontal line added.
plot_mixed <- ggplot(TDiff_PSIbin_df, aes(x = bin_median, y = avg_transpiration_deficit, color = species)) +
  geom_point() +
  geom_line(data = pred_data, aes(x = bin_median, y = predicted, color = species), size = 1) +
  geom_hline(yintercept = line_val, linetype = "dashed", color = "black", size = 1) +
  scale_color_manual(values = cb_palette) +
  labs(x = "Soil Water Potential (bin_median)",
       y = "Average Transpiration Deficit") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.background = element_rect(fill = "white", color = "white"),
    panel.background = element_rect(fill = "white"),
    legend.background = element_rect(fill = "white", color = "white"),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    panel.border = element_rect(color = NA, fill = NA, linewidth = 0),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    strip.text = element_text(face = "bold", size = 12)
  )

## Panel B: Bar Plot of Local Slopes at 30% of Maximum Transpiration Deficit

# For each species, compute the local slope at the point where the predicted value equals 30% of its maximum.
# We define a window of ±5 points around the index of the threshold.
local_slope_data <- pred_data %>%
  group_by(species) %>%
  do({
    df <- .
    max_pred <- max(df$predicted, na.rm = TRUE)
    thresh_val <- 0.3 * max_pred
    idx <- which.min(abs(df$predicted - thresh_val))
    # Define a window around the threshold index (adjust window size as needed)
    win_start <- max(1, idx - 5)
    win_end <- min(nrow(df), idx + 5)
    df_window <- df[win_start:win_end, ]
    # Fit a simple linear model on the local window.
    lm_local <- lm(predicted ~ bin_median, data = df_window)
    summ <- summary(lm_local)
    slope_val <- coef(lm_local)["bin_median"]
    slope_se <- summ$coefficients["bin_median", "Std. Error"]
    p_val <- summ$coefficients["bin_median", "Pr(>|t|)"]
    r2_val <- summ$r.squared
    data.frame(threshold_value = thresh_val,
               slope = slope_val,
               slope_se = slope_se,
               p_value = p_val,
               r2 = r2_val)
  }) %>%
  ungroup()

# Prepare the data for the bar plot.
local_slope_data <- local_slope_data %>%
  mutate(slope_abs = abs(slope),
         label_text = ifelse(p_value < 0.05,
                             paste0("p<0.05\nR²: ", round(r2, 3)),
                             paste0("Slope: ", round(slope, 4), "\n",
                                    "p: ", signif(p_value, 3), "\n",
                                    "R²: ", round(r2, 3))))

# Create the bar plot.
# The bars are ordered (from smallest to largest absolute slope) using reorder().
p_bar <- ggplot(local_slope_data, aes(x = reorder(species, slope_abs), y = slope_abs, fill = species)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_errorbar(aes(ymin = slope_abs - slope_se, ymax = slope_abs + slope_se),
                width = 0.2, color = "black") +
  geom_text(aes(label = label_text),
            position = position_stack(vjust = 0.5),
            size = 4, color = "black") +
  scale_fill_manual(values = cb_palette) +
  labs(x = "Species",
       y = "Absolute Slope at 30% Transpiration Deficit") +
  coord_flip() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_blank(),
    axis.text.y = element_text(color = "black"),
    plot.background = element_rect(fill = "white", color = "white"),
    panel.background = element_rect(fill = "white"),
    legend.background = element_rect(fill = "white", color = "white"),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    panel.border = element_rect(color = NA, fill = NA, linewidth = 0),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    strip.text = element_text(face = "bold", size = 12)
  )

## Combine Panel A and Panel B Side by Side

final_plot <- plot_mixed + p_bar + plot_layout(widths = c(2, 1))
print(final_plot)

# Save the final plot to a file.
ggsave("results/Figures/TDiff_PSIbin_non_linear.png", plot = final_plot, width = 12, height = 6, dpi = 300)
