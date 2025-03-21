##############################################
# NDVI/PSI Temporal Auto-correlation Analysis
# Author: Yixuan Wang
# Date: Sys.Date()
##############################################

# Load required libraries
library(dplyr)
library(tidyr)
library(terra)
library(stats)

# ---- Global Settings ----
years <- as.character(2003:2024)
species_list <- c("Beech", "Oak", "Spruce", "Pine")

# Define output folders (adjust paths if needed)
output_folder_raster <- "results_temporal/Raster"
output_folder_figures <- "results_temporal/Figures"
output_folder_data    <- "results/Data"

# Create output directories if they do not exist
dir.create(output_folder_raster, recursive = TRUE, showWarnings = FALSE)
dir.create(output_folder_figures, recursive = TRUE, showWarnings = FALSE)

# Define a color palette for species
cb_palette <- c("Oak"   = "#E69F00",  # Orange
                "Beech" = "#0072B2",  # Deep blue
                "Spruce"= "#009E73",  # Bluish-green
                "Pine"  = "#F0E442")  # Yellowish

# ---- Function Definitions ----

# 1. Detrending Function (returns residuals)
detrend_fun <- function(ts) {
  if (all(is.na(ts))) return(rep(NA, length(ts)))
  time <- 1:length(ts)
  model <- lm(ts ~ time)
  model$residuals
}

# 2. Spectral Exponent Function for Detrended Data
spectral_exponent_fun <- function(ts) {
  if (all(is.na(ts))) return(NA)
  time <- 1:length(ts)
  mod <- lm(ts ~ time)
  resid <- mod$residuals
  n <- length(resid)
  fft_vals <- fft(resid)
  power <- (Mod(fft_vals))^2
  freq <- (0:(n - 1)) / n
  pos_ind <- 2:floor(n / 2)
  log_freq <- log10(freq[pos_ind])
  log_power <- log10(power[pos_ind])
  fit <- lm(log_power ~ log_freq)
  coef(fit)[2]
}

# 3. Spectral Exponent Function for Original Data
spectral_exponent_original <- function(ts) {
  if (all(is.na(ts))) return(NA)
  n <- length(ts)
  fft_vals <- fft(ts)
  power <- (Mod(fft_vals))^2
  freq <- (0:(n - 1)) / n
  pos_ind <- 2:floor(n / 2)
  log_freq <- log10(freq[pos_ind])
  log_power <- log10(power[pos_ind])
  fit <- lm(log_power ~ log_freq)
  coef(fit)[2]
}

# 4. Create a SpatRaster from a filtered data frame for a given parameter
create_raster_from_df <- function(df, param) {
  df_wide <- df %>%
    pivot_wider(names_from = year, values_from = !!sym(param))
  rast(df_wide, type = "xyz")
}

# 5. Process and Save Spectral Exponent for Original Data
process_spectral_exponent_original <- function(sp, param, all_results_df) {
  message("Processing original spectral exponent for ", sp, " and parameter ", param)
  df_filtered <- all_results_df %>%
    filter(species == sp) %>%
    select(x, y, year, !!sym(param)) %>%
    mutate(year = as.numeric(year))
  
  r_param <- create_raster_from_df(df_filtered, param)
  
  se_orig <- app(r_param, spectral_exponent_original)
  out_raster <- file.path(output_folder_raster, paste0(sp, "_", param, "_original_spectral_exponent.tif"))
  writeRaster(se_orig, out_raster, overwrite = TRUE)
  
  # Save plot
  out_plot <- file.path(output_folder_figures, paste0(sp, "_", param, "_original_spectral_exponent.png"))
  png(filename = out_plot, width = 800, height = 600)
  plot(se_orig, main = paste("Original Spectral Exponent for", sp, param))
  dev.off()
  message("Saved original spectral exponent (raster & plot) for ", sp, " ", param)
}

# 6. Process Detrended Data: Save detrended raster and compute spectral exponent (detrended)
process_detrended_data <- function(sp, param, all_results_df) {
  message("Processing detrended data for ", sp, " and parameter ", param)
  df_filtered <- all_results_df %>%
    filter(species == sp) %>%
    select(x, y, year, !!sym(param)) %>%
    mutate(year = as.numeric(year))
  
  r_param <- create_raster_from_df(df_filtered, param)
  
  # Detrend the time series per pixel
  detrended_rast <- app(r_param, detrend_fun)
  out_det <- file.path(output_folder_raster, paste0(sp, "_", param, "_detrended.tif"))
  writeRaster(detrended_rast, out_det, overwrite = TRUE)
  message("Saved detrended raster for ", sp, " ", param)
  
  # Compute spectral exponent on detrended data
  se_det <- app(r_param, spectral_exponent_fun)
  out_se <- file.path(output_folder_raster, paste0(sp, "_", param, "_detrended_spectral_exponent.tif"))
  writeRaster(se_det, out_se, overwrite = TRUE)
  
  # Save plot
  out_se_plot <- file.path(output_folder_figures, paste0(sp, "_", param, "_detrended_spectral_exponent.png"))
  png(filename = out_se_plot, width = 800, height = 600)
  plot(se_det, main = paste("Detrended Spectral Exponent for", sp, param))
  dev.off()
  message("Saved detrended spectral exponent (raster & plot) for ", sp, " ", param)
}

# 7. Combine Detrended Data for Further Analysis
combine_detrended_data <- function(species_list, params, years, output_folder_data) {
  all_species_list <- list()
  for (sp in species_list) {
    message("Combining detrended data for species: ", sp)
    # Construct file paths for each parameter
    quant_file <- file.path(output_folder_raster, paste0(sp, "_Quantiles_detrended.tif"))
    psi_file   <- file.path(output_folder_raster, paste0(sp, "_soil_water_potential_detrended.tif"))
    tdiff_file <- file.path(output_folder_raster, paste0(sp, "_transpiration_deficit_detrended.tif"))
    
    # Load the detrended rasters and assign year names
    quant_rast <- rast(quant_file); names(quant_rast) <- years
    psi_rast   <- rast(psi_file);   names(psi_rast)   <- years
    tdiff_rast <- rast(tdiff_file); names(tdiff_rast) <- years
    
    # Process each year
    sp_year_list <- lapply(years, function(yr) {
      quant_layer <- quant_rast[[yr]]
      psi_layer   <- psi_rast[[yr]]
      tdiff_layer <- tdiff_rast[[yr]]
      
      # Convert layers to data frames with coordinates
      df_quant <- as.data.frame(quant_layer, xy = TRUE)
      colnames(df_quant)[3] <- "Quantiles"
      
      df_psi <- as.data.frame(psi_layer, xy = TRUE)
      colnames(df_psi)[3] <- "soil_water_potential"
      
      df_tdiff <- as.data.frame(tdiff_layer, xy = TRUE)
      colnames(df_tdiff)[3] <- "transpiration_deficit"
      
      # Merge by coordinates and add species and year info
      df_temp <- full_join(df_quant, df_psi, by = c("x", "y")) %>%
        full_join(df_tdiff, by = c("x", "y")) %>%
        mutate(year = yr, species = sp)
      df_temp
    })
    sp_df <- do.call(rbind, sp_year_list)
    all_species_list[[sp]] <- sp_df
    
    # Save individual CSV for each species
    out_csv <- file.path(output_folder_data, paste0(sp, "_detrended_results.csv"))
    write.csv(sp_df, out_csv, row.names = FALSE)
    message("Saved CSV for ", sp)
  }
  
  final_df <- do.call(rbind, all_species_list)
  final_df_clean <- na.omit(final_df)
  final_data_filename <- file.path(output_folder_data, "All_Quantiles_PSI_TDiff_species_year_detrended.RData")
  save(final_df_clean, file = final_data_filename)
  message("Saved combined detrended data to ", final_data_filename)
  final_df_clean
}

# 8. Plot ACF (Random Pixel & Boxplot) for a given raster stack
plot_acf_for_raster <- function(rast_obj, species, var, output_dir, suffix, cb_palette) {
  # Extract raster values (each row = a pixel’s time series)
  vals <- values(rast_obj)
  
  ## (a) Random Pixel ACF
  out_rand <- file.path(output_dir, paste0(species, "_", var, "_ACF_one_pixel_", suffix, ".png"))
  png(filename = out_rand, width = 1200, height = 800)
  valid_cells <- which(apply(vals, 1, function(x) all(!is.na(x))))
  if (length(valid_cells) == 0) {
    plot.new()
    title(main = paste(species, gsub("_", " ", var), "\nNo valid pixel data"), col.main = cb_palette[species])
  } else {
    cell_idx <- sample(valid_cells, 1)
    ts_values <- vals[cell_idx, ]
    acf(ts_values, main = paste(species, gsub("_", " ", var), "ACF for Cell", cell_idx, "(", suffix, ")"),
        col.main = cb_palette[species])
  }
  dev.off()
  
  ## (b) Boxplot of ACF for All Pixels
  out_box <- file.path(output_dir, paste0(species, "_", var, "_ACF_all_pixels_boxplot_", suffix, ".png"))
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
          main = paste("Distribution of", species, gsub("_", " ", var), "ACF (", suffix, ")"),
          xlab = "Lag",
          ylab = "ACF",
          col = cb_palette[species])
  abline(h = ci_bound, col = "blue", lty = 2)
  abline(h = -ci_bound, col = "blue", lty = 2)
  abline(h = 0, col = "red", lty = 2)
  dev.off()
  
  message("Saved ACF plots for ", species, " ", var, " (", suffix, ")")
}

# 9. Process ACF for a given variable across species
process_acf_for_variable <- function(var, species_order, all_results_df, output_dir, detrended = FALSE, cb_palette) {
  rast_list_species <- list()
  
  for (sp in species_order) {
    if (detrended) {
      file_path <- file.path(output_folder_raster, paste0(sp, "_", var, "_detrended.tif"))
      if (!file.exists(file_path)) {
        warning("File not found: ", file_path)
        next
      }
      rast_obj <- rast(file_path)
      names(rast_obj) <- years
    } else {
      sp_df <- all_results_df %>% filter(species == sp)
      sp_df$year <- as.character(sp_df$year)
      unique_years <- sort(unique(sp_df$year))
      rast_list <- lapply(unique_years, function(yr) {
        sub_df <- sp_df %>% filter(year == yr)
        rast(sub_df[, c("x", "y", var)], type = "xyz")
      })
      names(rast_list) <- unique_years
      rast_obj <- rast(rast_list)
    }
    rast_list_species[[sp]] <- rast_obj
  }
  
  # Suffix for file naming based on whether data are detrended
  suffix <- ifelse(detrended, "detrended", "original")
  for (sp in species_order) {
    rast_obj <- rast_list_species[[sp]]
    if (!is.null(rast_obj)) {
      plot_acf_for_raster(rast_obj, sp, var, output_dir, suffix, cb_palette)
    }
  }
}

