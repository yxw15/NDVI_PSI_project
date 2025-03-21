##############################################
# NDVI/PSI Temporal Auto-correlation Analysis (FFT-based)
# Combined Script
# Author: Yixuan Wang (adapted and reorganized)
# Date: Sys.Date()
##############################################
# Description:
# This script performs temporal auto-correlation analysis on NDVI/PSI datasets 
# using both original and FFT-detrended data. The analysis supports processing 
# for both "Quantiles" and "Proportions" data types. The workflow includes:
#
#  1. FFT-based Detrending & Spectral Analysis:
#       - Detrends the time series via FFT.
#       - Computes the spectral exponent from the detrended data.
#
#  2. Composite Spectral Exponent Plotting:
#       - Loads original and detrended spectral exponent rasters.
#       - Bins the values into noise categories and generates composite maps.
#
#  3. ACF (Autocorrelation Function) Analysis:
#       - Creates ACF plots for a random pixel and boxplots for all pixels.
##############################################

#### 1. SETUP & GLOBAL CONFIGURATION ####

# Set working directory
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

# Define output folders
output_folder_raster   <- "results_temporal/FFT/Raster"
output_folder_figures  <- "results_temporal/FFT/Figures"
output_folder_data     <- "results/Data"
output_dir_acf_results <- "results_temporal/FFT/results"  # For ACF plot outputs

# Create output directories if they don't exist
dirs <- c(output_folder_raster, output_folder_figures, output_folder_data, output_dir_acf_results)
for(d in dirs){
  if (!dir.exists(d)) { dir.create(d, recursive = TRUE) }
}

# Load required packages
library(dplyr)
library(tidyr)
library(terra)
library(ggplot2)
library(rlang)
library(stats)
library(lmtest)
library(forecast)
library(viridis)
library(grid)

# Global settings and parameters
years <- as.character(2003:2024)  # Analysis years
species_list <- c("Beech", "Oak", "Spruce", "Pine")  # Species for FFT processing
species_order <- c("Oak", "Beech", "Spruce", "Pine")  # Preferred order for ACF plots

# Define color palette for species
cb_palette <- c("Oak"    = "#E69F00",  # Orange
                "Beech"  = "#0072B2",  # Deep blue
                "Spruce" = "#009E73",  # Bluish-green
                "Pine"   = "#F0E442")  # Yellowish

# Data types to process: supports both "Quantiles" and "Proportions"
data_types <- c("Quantiles", "Proportions")

# Other variables (parameters) to process (common for both data types)
other_vars <- c("soil_water_potential", "transpiration_deficit")
# The complete set for each data type will be: data_type + other_vars


#### 2. FUNCTION DEFINITIONS ####

## --- 2A. FFT-BASED DETRENDING & SPECTRAL ANALYSIS FUNCTIONS ---

fft_detrend_fun <- function(ts) {
  if (all(is.na(ts))) return(rep(NA, length(ts)))
  n <- length(ts)
  fft_vals <- fft(ts)
  freq <- (0:(n - 1)) / n
  fft_vals_filtered <- fft_vals
  cutoff <- 0.1  # Adjust cutoff frequency as needed
  fft_vals_filtered[freq < cutoff] <- 0
  detrended_ts <- Re(fft(fft_vals_filtered, inverse = TRUE)) / n
  detrended_ts
}

spectral_exponent_fun <- function(ts) {
  if (all(is.na(ts))) return(NA)
  time <- 1:length(ts)
  mod <- lm(ts ~ time)
  resid <- mod$residuals
  n <- length(resid)
  fft_vals <- fft(resid)
  power <- (Mod(fft_vals))^2
  freq <- (0:(n - 1)) / n
  pos_ind <- 2:floor(n / 2)  # Exclude zero frequency and use positive frequencies
  log_freq <- log10(freq[pos_ind])
  log_power <- log10(power[pos_ind])
  fit <- lm(log_power ~ log_freq)
  coef(fit)[2]  # Return the slope (spectral exponent)
}

create_raster_from_df <- function(df, param) {
  df_wide <- df %>%
    pivot_wider(names_from = year, values_from = !!sym(param))
  rast(df_wide, type = "xyz")
}

process_detrended_data_fft <- function(sp, param, all_results_df) {
  message("Processing FFT-detrended data for ", sp, " and parameter ", param)
  
  df_filtered <- all_results_df %>%
    filter(species == sp) %>%
    select(x, y, year, !!sym(param)) %>%
    mutate(year = as.numeric(year))
  
  r_param <- create_raster_from_df(df_filtered, param)
  
  detrended_rast <- app(r_param, fft_detrend_fun)
  out_det <- file.path(output_folder_raster, paste0(sp, "_", param, "_fft_detrended.tif"))
  writeRaster(detrended_rast, out_det, overwrite = TRUE)
  message("Saved FFT-detrended raster for ", sp, " ", param)
  
  se_det <- app(detrended_rast, spectral_exponent_fun)
  out_se <- file.path(output_folder_raster, paste0(sp, "_", param, "_fft_detrended_spectral_exponent.tif"))
  writeRaster(se_det, out_se, overwrite = TRUE)
  message("Processed FFT-detrended spectral exponent for ", sp, " ", param)
}

# Added a new argument 'data_type' to avoid overwriting results
combine_detrended_data_fft <- function(species_list, params, years, output_folder_data, data_type) {
  all_species_list <- list()
  for (sp in species_list) {
    message("Combining FFT-detrended data for species: ", sp, " (", data_type, ")")
    
    # Construct file paths dynamically for each parameter
    file_paths <- lapply(params, function(param) {
      file.path(output_folder_raster, paste0(sp, "_", param, "_fft_detrended.tif"))
    })
    names(file_paths) <- params
    
    sp_year_list <- lapply(years, function(yr) {
      df_list <- lapply(params, function(param) {
        rast_obj <- rast(file_paths[[param]])
        names(rast_obj) <- years
        layer <- rast_obj[[yr]]
        df_temp <- as.data.frame(layer, xy = TRUE)
        colnames(df_temp)[3] <- param
        df_temp
      })
      df_temp <- Reduce(function(x, y) full_join(x, y, by = c("x", "y")), df_list)
      df_temp <- df_temp %>% mutate(year = yr, species = sp)
      df_temp
    })
    
    sp_df <- do.call(rbind, sp_year_list)
    all_species_list[[sp]] <- sp_df
    
    out_csv <- file.path(output_folder_data, paste0(sp, "_fft_detrended_results_", data_type, ".csv"))
    write.csv(sp_df, out_csv, row.names = FALSE)
    message("Saved CSV for ", sp, " (FFT-detrended, ", data_type, ")")
  }
  
  final_df <- do.call(rbind, all_species_list)
  final_df_clean <- na.omit(final_df)
  final_data_filename <- file.path(output_folder_data, paste0("All_", data_type, "_PSI_TDiff_species_year_fft_detrended.RData"))
  save(final_df_clean, file = final_data_filename)
  message("Saved combined FFT-detrended data to ", final_data_filename)
  final_df_clean
}

## --- 2B. ACF PLOTTING FUNCTIONS ---

plot_acf_for_raster <- function(rast_obj, species, var, output_dir, cb_palette) {
  vals <- values(rast_obj)
  
  # (a) ACF for a random pixel
  out_rand <- file.path(output_dir, paste0(species, "_", var, "_ACF_one_pixel_fft.png"))
  png(filename = out_rand, width = 1200, height = 800)
  valid_cells <- which(apply(vals, 1, function(x) all(!is.na(x))))
  if (length(valid_cells) == 0) {
    plot.new()
    title(main = paste(species, gsub("_", " ", var), "\nNo valid pixel data"),
          col.main = cb_palette[species])
  } else {
    cell_idx <- sample(valid_cells, 1)
    ts_values <- vals[cell_idx, ]
    acf(ts_values, main = paste(species, gsub("_", " ", var),
                                "ACF for Cell", cell_idx, "(FFT-detrended)"),
        col.main = cb_palette[species])
  }
  dev.off()
  
  # (b) Boxplot of ACF for all valid pixels
  out_box <- file.path(output_dir, paste0(species, "_", var, "_ACF_all_pixels_boxplot_fft.png"))
  png(filename = out_box, width = 1200, height = 800)
  valid_pixels <- apply(vals, 1, function(x) all(!is.na(x)))
  valid_vals <- vals[valid_pixels, ]
  n_years <- ncol(valid_vals)
  lag_max <- n_years - 1
  acf_matrix <- t(apply(valid_vals, 1, function(ts) {
    acf_res <- acf(ts, plot = FALSE, lag.max = lag_max)
    as.vector(acf_res$acf)
  }))
  lag_names <- paste0("lag_", 0:lag_max)
  acf_df <- as.data.frame(acf_matrix)
  colnames(acf_df) <- lag_names
  ci_bound <- 1.96 / sqrt(n_years)
  boxplot(acf_df,
          main = paste("Distribution of", species, gsub("_", " ", var),
                       "ACF (FFT-detrended)"),
          xlab = "Lag",
          ylab = "ACF",
          col = cb_palette[species])
  abline(h = ci_bound, col = "blue", lty = 2)
  abline(h = -ci_bound, col = "blue", lty = 2)
  abline(h = 0, col = "red", lty = 2)
  dev.off()
  
  message("Saved ACF plots for ", species, " ", var, " (FFT-detrended)")
}

process_acf_for_variable_fft <- function(var, species_order, output_dir, cb_palette) {
  for (sp in species_order) {
    file_path <- file.path(output_folder_raster, paste0(sp, "_", var, "_fft_detrended.tif"))
    if (!file.exists(file_path)) {
      warning("File not found: ", file_path)
      next
    }
    rast_obj <- rast(file_path)
    names(rast_obj) <- years
    plot_acf_for_raster(rast_obj, sp, var, output_dir, cb_palette)
  }
}

generate_acf_plots <- function(data, data_type, species_order, vars, output_dir, cb_palette) {
  data$species <- factor(data$species, levels = species_order)
  for (var in vars) {
    rast_list_species <- list()
    for (sp in species_order) {
      sp_df <- data %>% filter(species == sp)
      sp_df$year <- as.character(sp_df$year)
      yrs <- sort(unique(sp_df$year))
      rast_list <- lapply(yrs, function(yr) {
        sub_df <- sp_df %>% filter(year == yr)
        rast(sub_df[, c("x", "y", var)], type = "xyz")
      })
      names(rast_list) <- yrs
      sp_rast <- rast(rast_list)
      rast_list_species[[sp]] <- sp_rast
    }
    
    # Random Pixel ACF Plot
    output_file_rand <- file.path(output_dir, paste0(var, "_ACF_one_pixel_subpanel_", data_type, ".png"))
    png(filename = output_file_rand, width = 1200, height = 800)
    par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
    for (sp in species_order) {
      sp_rast <- rast_list_species[[sp]]
      raster_values <- values(sp_rast)
      valid_cells <- which(apply(raster_values, 1, function(x) all(!is.na(x))))
      if (length(valid_cells) == 0) {
        plot.new()
        title(main = paste(sp, var, "\nNo valid pixel data"), col.main = cb_palette[sp])
        next
      }
      cell_idx <- sample(valid_cells, 1)
      ts_values <- raster_values[cell_idx, ]
      acf(ts_values, main = paste(sp, var, "ACF for Cell", cell_idx), col.main = cb_palette[sp])
    }
    dev.off()
    
    # Boxplot of ACF for All Pixels
    output_file_box <- file.path(output_dir, paste0(var, "_ACF_all_pixels_boxplot_subpanel_", data_type, ".png"))
    png(filename = output_file_box, width = 1200, height = 800)
    par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
    for (sp in species_order) {
      sp_rast <- rast_list_species[[sp]]
      vals <- values(sp_rast)
      valid_pixels <- apply(vals, 1, function(x) all(!is.na(x)))
      valid_vals <- vals[valid_pixels, ]
      n_years <- ncol(valid_vals)
      lag_max <- n_years - 1
      acf_matrix <- t(apply(valid_vals, 1, function(ts) {
        acf_res <- acf(ts, plot = FALSE, lag.max = lag_max)
        as.vector(acf_res$acf)
      }))
      lag_names <- paste0("lag_", 0:lag_max)
      acf_df <- as.data.frame(acf_matrix)
      colnames(acf_df) <- lag_names
      ci_bound <- 1.96 / sqrt(n_years)
      boxplot(acf_df,
              main = paste("Distribution of", sp, var, "ACF"),
              xlab = "Lag",
              ylab = "ACF",
              col = cb_palette[sp])
      abline(h = ci_bound, col = "blue", lty = 2)
      abline(h = -ci_bound, col = "blue", lty = 2)
      abline(h = 0, col = "red", lty = 2)
    }
    dev.off()
  }
}

## --- 2C. COMPOSITE PLOTTING FUNCTIONS FOR SPECTRAL EXPONENT ---

get_species_param_df <- function(species, param) {
  orig_file <- file.path(output_folder_raster, paste0(species, "_", param, "_original_spectral_exponent.tif"))
  det_file  <- file.path(output_folder_raster, paste0(species, "_", param, "_fft_detrended_spectral_exponent.tif"))
  
  rast_orig <- rast(orig_file)
  rast_det  <- rast(det_file)
  
  df_orig <- as.data.frame(rast_orig, xy = TRUE)
  df_det  <- as.data.frame(rast_det, xy = TRUE)
  
  colnames(df_orig)[3] <- "original"
  colnames(df_det)[3]  <- "detrended"
  
  df_merge <- inner_join(df_orig, df_det, by = c("x", "y"))
  df_merge <- df_merge %>% mutate(diff = original - detrended,
                                  param = param,
                                  species = species)
  
  df_long <- df_merge %>%
    pivot_longer(cols = c("original", "detrended", "diff"),
                 names_to = "method",
                 values_to = "value")
  df_long
}

#### 3. MAIN EXECUTION ####

# Loop over each data type ("Quantiles" and "Proportions")
for (data_type in data_types) {
  
  message("Processing data type: ", data_type)
  
  # Load the corresponding original combined data
  if (data_type == "Quantiles") {
    load("results/Data/All_Species_Quantiles_PSI_TDiff.RData")  # loads all_results_df
    current_df <- all_results_df
  } else if (data_type == "Proportions") {
    load("results/Data/All_Species_Proportions_PSI_TDiff.RData")  # loads all_results_df
    current_df <- all_results_df
  }
  
  # Define parameter list for current data type
  params_list <- c(data_type, other_vars)
  
  # (A) FFT-BASED DETRENDING & SPECTRAL EXPONENT ANALYSIS
  for (sp in species_list) {
    for (param in params_list) {
      process_detrended_data_fft(sp, param, current_df)
    }
  }
  
  # (B) Combine the FFT-detrended data across species
  combined_df_fft <- combine_detrended_data_fft(species_list, params_list, years, output_folder_data, data_type)
  
  # (C) Generate ACF plots using the saved raster files for current data type
  for (var in params_list) {
    process_acf_for_variable_fft(var, species_order, output_folder_figures, cb_palette)
  }
  
  # (D) Generate ACF plots from combined data
  generate_acf_plots(data = current_df, 
                     data_type = paste0("original_", data_type),
                     species_order = species_order,
                     vars = params_list,
                     output_dir = output_dir_acf_results,
                     cb_palette = cb_palette)
  
  generate_acf_plots(data = combined_df_fft, 
                     data_type = paste0("fft_detrended_", data_type),
                     species_order = species_order,
                     vars = params_list,
                     output_dir = output_dir_acf_results,
                     cb_palette = cb_palette)
  
  # (E) Composite Spectral Exponent Plotting for current data type
  species_list_plot <- species_order  # or adjust if needed
  params_comp <- params_list
  
  all_data <- do.call(rbind,
                      lapply(species_list_plot, function(sp) {
                        do.call(rbind, lapply(params_comp, function(p) {
                          get_species_param_df(sp, p)
                        }))
                      }))
  
  all_data$species <- factor(all_data$species, levels = species_list_plot)
  all_data$param   <- factor(all_data$param, levels = params_comp)
  all_data$method  <- factor(all_data$method, levels = c("original", "detrended", "diff"),
                             labels = c("Original", "Detrended", "Difference"))
  
  all_data <- all_data %>%
    mutate(value_group = cut(value,
                             breaks = c(-Inf, -2, -1, 0, Inf),
                             labels = c("Red noise", "Pink noise", "White noise", "High-frequency")))
  
  magma_colors <- magma(4)
  
  for(sp in species_list_plot) {
    sp_data <- filter(all_data, species == sp)
    p <- ggplot(sp_data, aes(x = x, y = y, color = value_group)) +
      geom_point(size = 0.5) +
      facet_grid(param ~ method) +
      scale_color_manual(values = setNames(magma_colors,
                                           c("Red noise", "Pink noise", "White noise", "High-frequency"))) +
      coord_quickmap() +
      labs(title = paste(sp, data_type, "Spectral Exponent"),
           x = "Longitude", y = "Latitude", color = "Spectral Exponent Category") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.key.size = unit(1.5, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        panel.border = element_rect(fill = NA, color = "black")
      ) +
      guides(color = guide_legend(override.aes = list(size = 4)))
    
    output_file <- file.path(output_dir_acf_results, paste0(sp, "_Composite_Spectral_Exponent_Binned_", data_type, ".png"))
    ggsave(output_file, p, width = 15, height = 12, dpi = 300)
    message("Saved composite spectral exponent plot for ", sp, " (", data_type, ")")
  }
}
