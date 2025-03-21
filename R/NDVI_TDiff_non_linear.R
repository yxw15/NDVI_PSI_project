setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")
load("results/Data/NDVI_PSI_TDiff_species_df.RData")

# Load required libraries
library(lme4)      # For mixed-effects modeling
library(dplyr)     # For data manipulation
library(ggplot2)   # For plotting
library(gridExtra) # For arranging multiple plots
library(patchwork) # For combining ggplots

NDVI_TDiffbin_df <- NDVI_TDiffbin(NDVI_PSI_TDiff_species)
NDVI_TDiffbin_df <- na.omit(NDVI_TDiffbin_df)

# Define a custom color palette and set the species order.
cb_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442")
names(cb_palette) <- c("Oak", "Beech", "Spruce", "Pine")
species_levels <- c("Oak", "Beech", "Spruce", "Pine")
NDVI_TDiffbin_df$species <- factor(NDVI_TDiffbin_df$species, levels = species_levels)

## Panel A: Mixed-Effects Model Plot for avg_quantiles vs bin_median

# Fit a mixed-effects model using a third-order polynomial for bin_median as a fixed effect,
# with a random third-order polynomial by species.
model_quant <- lmer(avg_quantiles ~ poly(bin_median, 3) + (poly(bin_median, 3) | species),
                    data = NDVI_TDiffbin_df)

# Create prediction data for each species.
pred_data_quant <- NDVI_TDiffbin_df %>%
  group_by(species) %>%
  summarise(min_bin = min(bin_median),
            max_bin = max(bin_median),
            .groups = "drop") %>%
  group_by(species) %>%
  do(data.frame(bin_median = seq(.$min_bin, .$max_bin, length.out = 100))) %>%
  ungroup()

# Ensure the species factor is in the desired order.
pred_data_quant$species <- factor(pred_data_quant$species, levels = species_levels)

# Generate predicted values (including random effects).
pred_data_quant$predicted <- predict(model_quant, newdata = pred_data_quant, re.form = NULL)

# Set the reference line at quantile = 10.
line_val_quant <- 10

# Create Panel A plot.
plot_mixed_quant <- ggplot(NDVI_TDiffbin_df, aes(x = bin_median, y = avg_quantiles, color = species)) +
  geom_point() +
  geom_line(data = pred_data_quant, aes(x = bin_median, y = predicted, color = species), size = 1) +
  geom_hline(yintercept = line_val_quant, linetype = "dashed", color = "black", size = 1) +
  scale_color_manual(values = cb_palette) +
  labs(x = "Transpiration Deficit (bin_median)",
       y = "Average Quantiles") +
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

## Panel B: Bar Plot of Local Slopes at Predicted Quantiles = 10

# For each species, compute the local slope at the point where the predicted quantile is closest to 10.
local_slope_data_quant <- pred_data_quant %>%
  group_by(species) %>%
  do({
    df <- .
    thresh_val <- 10  # fixed threshold of 10 for predicted quantiles
    # Find the index where the predicted value is closest to the threshold.
    idx <- which.min(abs(df$predicted - thresh_val))
    # Define a window (±5 points) around the threshold index.
    win_start <- max(1, idx - 5)
    win_end <- min(nrow(df), idx + 5)
    df_window <- df[win_start:win_end, ]
    # Fit a local linear model.
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
local_slope_data_quant <- local_slope_data_quant %>%
  mutate(slope_abs = abs(slope),
         label_text = ifelse(p_value < 0.05,
                             paste0("p<0.05\nR²: ", round(r2, 3)),
                             paste0("Slope: ", round(slope, 4), "\n",
                                    "p: ", signif(p_value, 3), "\n",
                                    "R²: ", round(r2, 3))))

# Create the bar plot.
p_bar_quant <- ggplot(local_slope_data_quant, aes(x = reorder(species, slope_abs), y = slope_abs, fill = species)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_errorbar(aes(ymin = slope_abs - slope_se, ymax = slope_abs + slope_se), width = 0.2, color = "black") +
  geom_text(aes(label = label_text), position = position_stack(vjust = 0.5), size = 4, color = "black") +
  scale_fill_manual(values = cb_palette) +
  labs(x = "Species",
       y = "Absolute Slope at Predicted Quantiles = 10") +
  coord_flip() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
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

## Combine Panels A and B Side by Side

final_plot_quant <- plot_mixed_quant + p_bar_quant + plot_layout(widths = c(2, 1))
print(final_plot_quant)

# Save the final combined plot.
ggsave("results/Figures/NDVI_TDiffbin_non_linear.png", plot = final_plot_quant, width = 12, height = 6, dpi = 300)
