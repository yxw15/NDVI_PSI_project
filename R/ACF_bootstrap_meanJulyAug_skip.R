# =============================================================================
# Script Name: Multi-Species Pixel-wise ACF Analysis for Soil Water Potential 
#              (Dataframe Version) with Skip Logic and Depth-specific CSV Folders
#
# Script Description:
#   This script analyzes time series data for several parameters across four species 
#   (Oak, Beech, Spruce, Pine) using an aggregated dataframe. Each pixel (defined by 
#   its x and y coordinates) has yearly measurements that are first averaged and then 
#   used to compute the autocorrelation function (ACF) for that pixel's time series.
#   This version skips re-processing if the corresponding output PNG already exists,
#   and writes CSV outputs into depth-specific subfolders under "results_acf".
#
# Analysis steps:
#   1. Data Preparation:
#      - Filter data for July and August and pivot it so that each pixel forms a 
#        complete time series (columns = years).
#   2. Computing the ACF:
#      - Compute the ACF for each pixel’s time series up to a specified maximum lag.
#   3. Bootstrapping Confidence Intervals:
#      - Use bootstrap resampling on each pixel’s time series to generate confidence 
#        intervals for the ACF.
#   4. Reshaping Data & Saving CSVs:
#      - Reshape the ACF values and bootstrap confidence intervals into a long format.
#      - Save ACF and CI summary CSVs in depth-specific folders under "results_acf/depth_<cm>".
#   5. Plotting:
#      - Generate a boxplot of pixel-wise ACF values for each lag with a ribbon showing 
#        the 25th and 75th percentiles.
#   6. Saving the Output:
#      - Save each plot in "results/Figures/depth<cm>/" unless it already exists.
#
# Required Libraries:
#   - dplyr, ggplot2, tidyr
#
# Usage:
#   - Adjust working directory and file paths as necessary.
#   - Ensure the aggregated dataframe (combined) contains columns:
#         x, y, year, Quantiles, soil_water_potential, transpiration_deficit, species, month
# =============================================================================

# Clear workspace (optional)
rm(list = ls())

# Set working directory (adjust if needed)
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

# Load required libraries
library(dplyr)
library(ggplot2)
library(tidyr)

# Define color palette for species
cb_palette <- c(
  "Oak"   = "#E69F00",  # Orange 🌳
  "Beech" = "#0072B2",  # Deep blue 🌿
  "Spruce"= "#009E73",  # Bluish-green 🌲
  "Pine"  = "#F0E442"   # Yellow 🌲
)

# =============================================================================
# Functions for Dataframe Analysis
# =============================================================================

pivotTimeSeries <- function(df, parameter_col) {
  message("🚀 Pivoting dataframe to time series matrix...")
  ts_df <- df %>% 
    select(x, y, year, !!sym(parameter_col)) %>%
    arrange(x, y, year)
  ts_wide <- pivot_wider(ts_df, names_from = year, values_from = !!sym(parameter_col))
  ts_matrix <- as.matrix(ts_wide[, -c(1,2)])
  rownames(ts_matrix) <- paste(ts_wide$x, ts_wide$y, sep = "_")
  valid_ts_matrix <- ts_matrix[apply(ts_matrix, 1, function(x) all(!is.na(x))), ]
  message("✅ Found ", nrow(valid_ts_matrix), " valid pixel series out of ", nrow(ts_matrix), " 🐾")
  return(valid_ts_matrix)
}

computeACFMatrix <- function(valid_vals, lag.max) {
  message("🔄 Computing ACF for each pixel time series...")
  acf_matrix <- t(apply(valid_vals, 1, function(ts) {
    acf(ts, plot = FALSE, lag.max = lag.max)$acf
  }))
  message("✅ ACF matrix computed! 💪")
  return(acf_matrix)
}

bootstrapCI <- function(valid_vals, lag.max, n_boot = 1000, alpha = 0.05) {
  message("💫 Bootstrapping confidence intervals...")
  set.seed(123)
  n_pix <- nrow(valid_vals)
  boot_lower <- matrix(NA, nrow = n_pix, ncol = lag.max + 1)
  boot_upper <- matrix(NA, nrow = n_pix, ncol = lag.max + 1)
  for (i in seq_len(n_pix)) {
    message(sprintf("🔢 Bootstrapping pixel %d of %d... 🧩", i, n_pix))
    ts <- valid_vals[i, ]
    acf_boot <- replicate(n_boot, {
      ts_perm <- sample(ts)
      as.numeric(acf(ts_perm, plot = FALSE, lag.max = lag.max)$acf)
    })
    boot_lower[i, ] <- apply(acf_boot, 1, quantile, probs = alpha/2)
    boot_upper[i, ] <- apply(acf_boot, 1, quantile, probs = 1-alpha/2)
  }
  message("✅ Bootstrapping complete! 🎉")
  return(list(lower = boot_lower, upper = boot_upper))
}

reshapeACFData <- function(acf_matrix, boot_ci_lower, boot_ci_upper,
                           species, var_short, out_dir_base = "results_acf",
                           depth = NULL) {
  depth_dir <- if (!is.null(depth)) file.path(out_dir_base, paste0("depth_", depth)) else out_dir_base
  if (!dir.exists(depth_dir)) dir.create(depth_dir, recursive = TRUE)
  message("📂 Saving CSVs to ", depth_dir)
  
  acf_df   <- as.data.frame(acf_matrix)
  lower_df <- as.data.frame(boot_ci_lower)
  upper_df <- as.data.frame(boot_ci_upper)
  acf_df$Pixel   <- rownames(acf_matrix)
  lower_df$Pixel <- rownames(acf_matrix)
  upper_df$Pixel <- rownames(acf_matrix)
  acf_long   <- pivot_longer(acf_df,   -Pixel, names_to = "Lag", values_to = "ACF")
  lower_long <- pivot_longer(lower_df, -Pixel, names_to = "Lag", values_to = "Lower_CI")
  upper_long <- pivot_longer(upper_df, -Pixel, names_to = "Lag", values_to = "Upper_CI")
  acf_combined <- left_join(acf_long, lower_long, by = c("Pixel","Lag")) %>%
    left_join(upper_long, by = c("Pixel","Lag"))
  acf_combined$LagNum <- as.integer(sub("V", "", acf_combined$Lag)) - 1
  
  ci_summary <- acf_combined %>%
    group_by(LagNum) %>%
    summarise(
      ci_lower = quantile(Lower_CI, 0.25, na.rm = TRUE),
      ci_upper = quantile(Upper_CI, 0.75, na.rm = TRUE)
    )
  
  write.csv(acf_combined, file = file.path(depth_dir, paste0(species, "_", var_short, "_acf_combined.csv")), row.names = FALSE)
  write.csv(ci_summary,   file = file.path(depth_dir, paste0(species, "_", var_short, "_ci_summary.csv")), row.names = FALSE)
  message("✅ CSVs saved for ", species, "_", var_short, "! 🌟")
  
  return(list(acf_combined = acf_combined, ci_summary = ci_summary))
}

generateACFPlot <- function(acf_combined, ci_summary, lag.max, species, var_short, cb_palette) {
  message("📊 Generating plot for ", species, " - ", var_short, "...")
  acf_combined_filtered <- acf_combined %>% filter(LagNum > 0)
  ci_summary_filtered    <- ci_summary    %>% filter(LagNum > 0)
  p <- ggplot(acf_combined_filtered, aes(x = LagNum, y = ACF)) +
    geom_ribbon(data = ci_summary_filtered,
                aes(x = LagNum, ymin = ci_lower, ymax = ci_upper),
                fill = cb_palette[species], alpha = 0.3, inherit.aes = FALSE) +
    geom_boxplot(aes(group = LagNum), fill = cb_palette[species], alpha = 0.6, outlier.size = 0.5) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    scale_x_continuous(breaks = sort(unique(acf_combined_filtered$LagNum)),
                       labels = sort(unique(acf_combined_filtered$LagNum))) +
    coord_cartesian(ylim = c(-1, 1)) +
    labs(title = "",
         x = "Lag (year)", y = "ACF") +
    theme_minimal(base_size = 14) +
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
  message("✅ Plot ready! 🖼️")
  return(p)
}

saveACFPlot <- function(plot, species, var_short, out_dir, width = 12, height = 6) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  outfile <- file.path(out_dir, paste0(species, "_", var_short, "_ACF_CI_Boxplot.png"))
  ggsave(outfile, plot, width = width, height = height)
  message("💾 Saved plot to ", outfile)
}

processDataFrame <- function(df, parameter_col, lag.max = 21, n_boot = 1000, alpha = 0.05,
                             var_short, species, cb_palette,
                             fig_out_dir = "results/Figures",
                             acf_out_dir = "results_acf",
                             depth = NULL) {
  message("🌟 Starting processing for ", species, ", variable=", var_short, ", depth=", depth, "cm")
  ts_mat   <- pivotTimeSeries(df, parameter_col)
  lag_act  <- min(lag.max, ncol(ts_mat) - 1)
  acf_mat  <- computeACFMatrix(ts_mat, lag_act)
  ci       <- bootstrapCI(ts_mat, lag_act, n_boot, alpha)
  reshaped <- reshapeACFData(acf_mat, ci$lower, ci$upper,
                             species, var_short,
                             out_dir_base = acf_out_dir,
                             depth = depth)
  p <- generateACFPlot(reshaped$acf_combined, reshaped$ci_summary, lag_act, species, var_short, cb_palette)
  saveACFPlot(p, species, var_short, file.path(fig_out_dir, paste0("depth", depth)))
  message("🎯 Completed processing for ", species, " variable=", var_short, " at depth=", depth)
}

# =============================================================================
# Execution Loop: Species, Parameters, Depth
# =============================================================================
message("🚀 Beginning full workflow!")
species_list <- c("Oak", "Beech", "Pine", "Spruce")
parameters   <- c("mean_soil_water_potential", "mean_transpiration_deficit", "mean_Quantiles")
var_names    <- c(mean_soil_water_potential = "PSI",
                  mean_transpiration_deficit  = "TDiff",
                  mean_Quantiles              = "Q")
depths <- c(150, 100, 50)

for (depth in depths) {
  message("---🌿 Depth: ", depth, "cm ---")
  load(file.path("results/Data", paste0("AllSpecies_AllMonths_depth", depth, ".RData")))
  fig_dir <- file.path("results/Figures", paste0("depth", depth))
  if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
  
  for (sp in species_list) {
    message("➡️  Species: ", sp)
    df_sp <- combined %>%
      filter(month %in% c("July","August"), species == sp)
    df_mean_sp <- df_sp %>%
      group_by(x, y, year) %>%
      summarise(
        mean_Quantiles            = mean(Quantiles, na.rm = TRUE),
        mean_soil_water_potential = mean(soil_water_potential, na.rm = TRUE),
        mean_transpiration_deficit= mean(transpiration_deficit, na.rm = TRUE),
        species = first(species)
      ) %>% ungroup()
    
    for (param in parameters) {
      var_short <- var_names[param]
      out_png   <- file.path(fig_dir, paste0(sp, "_", var_short, "_ACF_CI_Boxplot.png"))
      if (file.exists(out_png)) {
        message("⏭️  Skipping existing plot: ", out_png)
      } else {
        processDataFrame(
          df            = df_mean_sp,
          parameter_col = param,
          lag.max       = 22,
          n_boot        = 1000,
          alpha         = 0.05,
          var_short     = var_short,
          species       = sp,
          cb_palette    = cb_palette,
          fig_out_dir   = "results/Figures",
          acf_out_dir   = "results_acf",
          depth         = depth
        )
      }
    }
  }
}
message("🎉 Workflow complete! All done! 🌟")
