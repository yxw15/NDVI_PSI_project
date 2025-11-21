# =============================================================================
# Multi-Species Pixel-wise ACF Analysis
# =============================================================================

rm(list = ls())

# Set working directory (adjust if needed)
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

load("results/Data/AllSpecies_AllMonths_rootzone.RData") # provides `combined`

library(tidyverse)
library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)  # parse_number()

# Define color palette for species
cb_palette <- c(
  "Oak"    = "#E69F00",
  "Beech"  = "#0072B2",
  "Spruce" = "#009E73",
  "Pine"   = "#F0E442"
)

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

pivotTimeSeries <- function(df, parameter_col) {
  message("🚀 Pivoting dataframe to time series matrix...")
  # Ensure year is character (so column names are valid)
  ts_df <- df %>%
    mutate(year = as.character(year)) %>%
    select(x, y, year, !!sym(parameter_col)) %>%
    arrange(x, y, year)
  ts_wide <- pivot_wider(ts_df, names_from = year, values_from = !!sym(parameter_col))
  
  # If no rows, return NULL
  if (nrow(ts_wide) == 0) {
    message("⚠️ No rows in pivot_wider result.")
    return(NULL)
  }
  
  # Convert to matrix (drop x,y)
  if (ncol(ts_wide) <= 2) {
    message("⚠️ Pivot produced no year columns.")
    return(NULL)
  }
  
  ts_matrix <- as.matrix(ts_wide[, -(1:2)])
  rownames(ts_matrix) <- paste(ts_wide$x, ts_wide$y, sep = "_")
  
  # keep only rows with no missing values (complete series)
  good_idx <- apply(ts_matrix, 1, function(x) all(!is.na(x)))
  valid_ts_matrix <- ts_matrix[good_idx, , drop = FALSE]
  message("✅ Found ", nrow(valid_ts_matrix), " valid pixel series out of ", nrow(ts_matrix), " 🐾")
  return(valid_ts_matrix)
}

computeACFMatrix <- function(valid_vals, lag.max) {
  message("🔄 Computing ACF for each pixel time series...")
  if (is.null(valid_vals) || nrow(valid_vals) == 0) {
    warning("No valid time series provided to computeACFMatrix(). Returning NULL.")
    return(NULL)
  }
  # Ensure lag.max <= ncol - 1
  lag_act <- min(lag.max, ncol(valid_vals) - 1)
  n_pix <- nrow(valid_vals)
  
  acf_matrix <- matrix(NA_real_, nrow = n_pix, ncol = lag_act + 1)
  rownames(acf_matrix) <- rownames(valid_vals)
  
  for (i in seq_len(n_pix)) {
    ts <- as.numeric(valid_vals[i, ])
    # protect from constant series or errors
    acf_vals <- tryCatch({
      a <- acf(ts, plot = FALSE, lag.max = lag_act)$acf
      as.numeric(a)
    }, error = function(e) {
      # fallback: if constant series, set lag0 = 1, others 0
      rep(0, lag_act + 1)
    })
    # ensure length correct
    if (length(acf_vals) < (lag_act + 1)) {
      acf_vals <- c(acf_vals, rep(NA_real_, (lag_act + 1) - length(acf_vals)))
    }
    acf_matrix[i, ] <- acf_vals
  }
  message("✅ ACF matrix computed! (lags 0:", lag_act, ")")
  return(acf_matrix)
}

bootstrapCI <- function(valid_vals, lag.max, n_boot = 500, alpha = 0.05) {
  message("💫 Bootstrapping confidence intervals (n_boot = ", n_boot, ") ...")
  if (is.null(valid_vals) || nrow(valid_vals) == 0) {
    warning("No valid time series provided to bootstrapCI(). Returning NULL.")
    return(NULL)
  }
  set.seed(123)
  n_pix <- nrow(valid_vals)
  lag_act <- min(lag.max, ncol(valid_vals) - 1)
  
  boot_lower <- matrix(NA_real_, nrow = n_pix, ncol = lag_act + 1)
  boot_upper <- matrix(NA_real_, nrow = n_pix, ncol = lag_act + 1)
  
  pb <- txtProgressBar(min = 0, max = n_pix, style = 3)
  for (i in seq_len(n_pix)) {
    ts <- as.numeric(valid_vals[i, ])
    # Pre-allocate matrix to collect bootstrap acf results (rows = lags, cols = replicates)
    # We'll replicate with replacement (true bootstrap)
    acf_boot <- replicate(n_boot, {
      ts_boot <- sample(ts, size = length(ts), replace = TRUE)
      vals <- tryCatch({
        as.numeric(acf(ts_boot, plot = FALSE, lag.max = lag_act)$acf)
      }, error = function(e) {
        # fallback for constant series
        c(1, rep(0, lag_act))
      })
      # ensure length
      if (length(vals) < (lag_act + 1)) vals <- c(vals, rep(NA_real_, (lag_act + 1) - length(vals)))
      vals
    })
    # acf_boot: matrix (lags+1) x n_boot
    boot_lower[i, ] <- apply(acf_boot, 1, quantile, probs = alpha/2, na.rm = TRUE)
    boot_upper[i, ] <- apply(acf_boot, 1, quantile, probs = 1 - alpha/2, na.rm = TRUE)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  message("✅ Bootstrapping complete! 🎉")
  rownames(boot_lower) <- rownames(valid_vals)
  rownames(boot_upper) <- rownames(valid_vals)
  return(list(lower = boot_lower, upper = boot_upper))
}

reshapeACFData <- function(acf_matrix, boot_ci_lower, boot_ci_upper, species, var_short, out_dir_base = "results_acf") {
  message("📂 Saving CSVs to ", out_dir_base)
  if (is.null(acf_matrix) || nrow(acf_matrix) == 0) {
    warning("No ACF matrix to reshape; skipping CSV write.")
    return(NULL)
  }
  if (!dir.exists(out_dir_base)) dir.create(out_dir_base, recursive = TRUE)
  
  acf_df   <- as.data.frame(acf_matrix)
  lower_df <- as.data.frame(boot_ci_lower)
  upper_df <- as.data.frame(boot_ci_upper)
  
  acf_df$Pixel   <- rownames(acf_matrix)
  lower_df$Pixel <- rownames(acf_matrix)
  upper_df$Pixel <- rownames(acf_matrix)
  
  acf_long   <- pivot_longer(acf_df,   -Pixel, names_to = "Lag", values_to = "ACF")
  lower_long <- pivot_longer(lower_df, -Pixel, names_to = "Lag", values_to = "Lower_CI")
  upper_long <- pivot_longer(upper_df, -Pixel, names_to = "Lag", values_to = "Upper_CI")
  
  acf_combined <- left_join(acf_long, lower_long, by = c("Pixel", "Lag")) %>%
    left_join(upper_long, by = c("Pixel", "Lag"))
  
  # compute numeric lag number robustly: parse digits from column name
  acf_combined <- acf_combined %>%
    mutate(LagNum = readr::parse_number(Lag) - 1)
  
  ci_summary <- acf_combined %>%
    group_by(LagNum) %>%
    summarise(
      ci_lower = quantile(Lower_CI, 0.25, na.rm = TRUE),
      ci_upper = quantile(Upper_CI, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
  
  # write CSVs
  fn_acf <- file.path(out_dir_base, paste0(species, "_", var_short, "_acf_combined.csv"))
  fn_ci  <- file.path(out_dir_base, paste0(species, "_", var_short, "_ci_summary.csv"))
  write.csv(acf_combined, file = fn_acf, row.names = FALSE)
  write.csv(ci_summary,   file = fn_ci,  row.names = FALSE)
  message("✅ CSVs saved:", " ", fn_acf, " and ", fn_ci)
  return(list(acf_combined = acf_combined, ci_summary = ci_summary))
}

generateACFPlot <- function(acf_combined, ci_summary, species, var_short, cb_palette) {
  message("📊 Generating plot for ", species, " - ", var_short, "...")
  if (is.null(acf_combined) || nrow(acf_combined) == 0) {
    warning("No ACF combined data to plot.")
    return(NULL)
  }
  acf_combined_filtered <- acf_combined %>% filter(LagNum > 0)
  ci_summary_filtered    <- ci_summary %>% filter(LagNum > 0)
  p <- ggplot(acf_combined_filtered, aes(x = LagNum, y = ACF)) +
    geom_ribbon(data = ci_summary_filtered,
                aes(x = LagNum, ymin = ci_lower, ymax = ci_upper),
                inherit.aes = FALSE, fill = cb_palette[species], alpha = 0.25) +
    geom_boxplot(aes(group = LagNum), fill = cb_palette[species], alpha = 0.6, outlier.size = 0.5) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    scale_x_continuous(breaks = sort(unique(acf_combined_filtered$LagNum))) +
    coord_cartesian(ylim = c(-1, 1)) +
    labs(x = "Lag (years)", y = "ACF") +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.8),  # ← keeps border
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(face = "bold", size = 14),
      axis.text = element_text(color = "black", size = 12)
    )
  message("✅ Plot ready! 🖼️")
  return(p)
}

saveACFPlot <- function(plot, species, var_short, out_dir, width = 12, height = 6) {
  if (is.null(plot)) {
    warning("No plot provided to saveACFPlot(). Skipping.")
    return(NULL)
  }
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  outfile <- file.path(out_dir, paste0(species, "_", var_short, "_ACF_CI_Boxplot.png"))
  ggsave(outfile, plot = plot, width = width, height = height, dpi = 300)
  message("💾 Saved plot to ", outfile)
  return(outfile)
}

processDataFrame <- function(df, parameter_col, lag.max = 21, n_boot = 500, alpha = 0.05,
                             var_short, species, cb_palette,
                             fig_out_dir = "results_rootzone/Figures_till2022/acf",
                             acf_out_dir = "results_acf") {
  message("🌟 Starting processing for ", species, ", variable=", var_short)
  ts_mat <- pivotTimeSeries(df, parameter_col)
  if (is.null(ts_mat) || nrow(ts_mat) == 0) {
    message("⏭️ No valid time series for ", species, " / ", var_short, " — skipping.")
    return(invisible(NULL))
  }
  
  lag_act <- min(lag.max, ncol(ts_mat) - 1)
  acf_mat <- computeACFMatrix(ts_mat, lag_act)
  ci <- bootstrapCI(ts_mat, lag_act, n_boot = n_boot, alpha = alpha)
  reshaped <- reshapeACFData(acf_mat, ci$lower, ci$upper, species, var_short, out_dir_base = acf_out_dir)
  
  p <- generateACFPlot(reshaped$acf_combined, reshaped$ci_summary, species, var_short, cb_palette)
  out_png <- saveACFPlot(p, species, var_short, out_dir = fig_out_dir)
  message("🎯 Completed processing for ", species, " variable=", var_short)
  return(list(plot = p, png = out_png, acf_csv = reshaped$acf_combined))
}

# ---------------------------------------------------------------------------
# Execution
# ---------------------------------------------------------------------------

message("🚀 Beginning full workflow!")

# species_list <- c("Pine")
# parameters   <- c("mean_soil_water_potential", "mean_transpiration_deficit", "mean_Quantiles")
# var_names    <- c(mean_soil_water_potential = "PSI",
#                   mean_transpiration_deficit  = "TDiff",
#                   mean_Quantiles              = "Q")

species_list <- c("Spruce")
parameters   <- c("mean_soil_water_potential","mean_Quantiles")
var_names    <- c(mean_soil_water_potential = "PSI",
                  mean_Quantiles            = "Q")

fig_dir <- file.path("results_rootzone/Figures_till2022/acf")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
if (!dir.exists("results_acf")) dir.create("results_acf", recursive = TRUE)

# Filter original combined data
data <- combined %>%
  filter(month %in% c("July", "August")) %>%
  filter(year < 2023) %>%
  filter(Quantiles > 0) %>%
  na.omit()

for (sp in species_list) {
  message("➡️  Species: ", sp)
  df_sp <- data %>%
    filter(month %in% c("July", "August"), species == sp)
  
  df_mean_sp <- df_sp %>%
    group_by(x, y, year) %>%
    summarise(
      mean_Quantiles             = mean(Quantiles, na.rm = TRUE),
      mean_soil_water_potential  = mean(soil_water_potential, na.rm = TRUE),
      mean_transpiration_deficit = mean(transpiration_deficit, na.rm = TRUE),
      species = first(species),
      .groups = "drop"
    )
  
  for (param in parameters) {
    var_short <- var_names[param]
    out_png   <- file.path(fig_dir, paste0(sp, "_", var_short, "_ACF_CI_Boxplot.png"))
    if (file.exists(out_png)) {
      message("⏭️  Skipping existing plot: ", out_png)
    } else {
      processDataFrame(
        df            = df_mean_sp,
        parameter_col = param,
        lag.max       = 19,
        n_boot        = 1000,
        alpha         = 0.05,
        var_short     = var_short,
        species       = sp,
        cb_palette    = cb_palette,
        fig_out_dir   = fig_dir,
        acf_out_dir   = "results_acf"
      )
    }
  }
}

message("🎉 Workflow complete! All done! 🌟")
