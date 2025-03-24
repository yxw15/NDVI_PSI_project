library(ncdf4)
library(reshape2)
library(ggplot2)
library(terra)
library(dplyr)
library(tidyr)

transfer_tdiff_to_df <- function(nc_file, start_date) {
  # Step 1: Open the NetCDF file
  nc <- nc_open(nc_file)
  
  # Step 2: Extract variables
  time <- ncvar_get(nc, "time")
  x <- ncvar_get(nc, "x")
  y <- ncvar_get(nc, "y")
  tdiff <- ncvar_get(nc, "tdiff")
  
  # Convert time from "days since" format to actual dates
  dates <- as.Date(start_date) + time
  
  # Close the NetCDF file
  nc_close(nc)
  
  # Step 3: Assign meaningful dimension names
  dimnames(tdiff) <- list(
    x = x,
    y = y,
    time = as.character(dates)
  )
  
  # Step 4: Transform the 3D array into a data frame
  tdiff_melted <- melt(tdiff, varnames = c("x", "y", "time"), value.name = "transpiration_deficit")
  
  tdiff_melted <- na.omit(tdiff_melted)
  return(tdiff_melted)
}

save_tdiff_df <- function(tdiff_melted, file_path) {
  write.csv(tdiff_melted, file_path, row.names = FALSE)
}

filter_tdiff <- function(tdiff_melted, month_day) {
  tdiff_melted <- na.omit(tdiff_melted)
  tdiff_filter <- subset(tdiff_melted, 
                         format(as.Date(time), "%m-%d") == month_day)
  return(tdiff_filter)
}

save_tdiff_raster <- function(tdiff_melted, month_day, output_dir) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  tdiff <- filter_tdiff(tdiff_melted, month_day)
  tdiff$time <- as.Date(tdiff$time, format = "%Y-%m-%d")
  years <- 2003:2024
  for (year in years) {
    tdiff_year <- tdiff %>% filter(format(time, "%Y") == as.character(year))
    tdiff_selected <- tdiff_year[, c("x", "y", "transpiration_deficit")]
    tdiff_raster <- rast(tdiff_selected, type = "xyz")
    crs(tdiff_raster) <- "epsg:31467"
    output_path <- file.path(output_dir, paste0("tdiff_", year, ".tif"))
    writeRaster(tdiff_raster, output_path, overwrite = TRUE)
    cat("Raster saved for year:", year, "at", output_path, "\n")
  }
}

combine_NDVI_PSI_TDiff_raster_to_df <- function(species_list) {
  library(terra)
  library(dplyr)
  library(tidyr)
  
  # Extract species name
  species_name <- basename(species_list$output_dir)
  
  # Define folder paths
  NDVI_mask_folder <- file.path(species_list$output_dir, "NDVI_mask")
  PSI_mask_folder <- file.path(species_list$output_dir, "PSI_mask")
  TDiff_mask_folder <- file.path(species_list$output_dir, "TDiff_mask")
  
  # List and load raster files dynamically (supporting both Quantiles and Proportions)
  NDVI_files <- list.files(NDVI_mask_folder, pattern = "\\.tif$", full.names = TRUE)
  PSI_files <- list.files(PSI_mask_folder, pattern = "\\.tif$", full.names = TRUE)
  TDiff_files <- list.files(TDiff_mask_folder, pattern = "\\.tif$", full.names = TRUE)
  
  # Check if Quantiles or Proportions exist
  NDVI_type <- if (any(grepl("Quantiles", NDVI_files))) "Quantiles" else "Proportions"
  
  # Read rasters as stack
  NDVI_stack <- rast(NDVI_files)
  PSI_stack <- rast(PSI_files)
  TDiff_stack <- rast(TDiff_files)
  
  # Assign names
  years <- 2003:2024
  names(NDVI_stack) <- paste0(NDVI_type, "_", years)
  names(PSI_stack) <- paste0("PSI_", years)
  names(TDiff_stack) <- paste0("TDiff_", years)
  
  # Combine all layers into one stack
  Species_Quantiles_PSI_TDiff <- c(NDVI_stack, PSI_stack, TDiff_stack)
  
  # Convert to dataframe with coordinates
  Species_df <- as.data.frame(Species_Quantiles_PSI_TDiff, xy = TRUE)
  
  # Reshape to long format
  Species_df_long <- Species_df %>%
    pivot_longer(cols = -c(x, y), names_to = "variable", values_to = "value") %>%
    mutate(
      year = gsub(".*_(\\d{4})", "\\1", variable),  # Extract year
      variable = gsub("_(\\d{4})", "", variable)    # Extract variable name
    ) %>%
    pivot_wider(names_from = variable, values_from = value) %>%  # Reshape to wide format
    mutate(species = species_name)  # Add species column
  
  return(Species_df_long)
}

plot_box_TDiff_PSI <- function(data, save_path = NULL) {
  library(dplyr)
  library(ggplot2)
  
  # 1. Bin the Soil Water Potential into 100-unit ranges and compute midpoints
  bin_breaks <- seq(floor(min(data$soil_water_potential, na.rm = TRUE)), 
                    ceiling(max(data$soil_water_potential, na.rm = TRUE)), 
                    by = 100)
  
  bin_midpoints <- bin_breaks[-length(bin_breaks)] + diff(bin_breaks) / 2  # Compute midpoints
  
  data <- data %>%
    mutate(soil_bin = cut(soil_water_potential, 
                          breaks = bin_breaks, 
                          include.lowest = TRUE, 
                          right = FALSE, 
                          labels = bin_midpoints))  # Assign midpoints as labels
  
  # Convert to numeric for proper ordering
  data$soil_bin <- as.numeric(as.character(data$soil_bin))
  
  # 2. Reorder Species
  data$species <- factor(data$species, 
                         levels = c("Oak", "Beech", "Spruce", "Pine"))
  
  # 3. Plot the Box Plots with color-blind friendly palette
  cb_palette <- c("Oak"   = "#E69F00",  # Orange
                  "Beech" = "#0072B2",  # A deep blue with lower luminance than "#56B4E9"
                  "Spruce"= "#009E73",  # Bluish-green
                  "Pine"  = "#F0E442") 
  
  
  plot <- ggplot(data, aes(x = factor(soil_bin), y = transpiration_deficit, fill = species)) +
    geom_boxplot() +
    facet_wrap(~ species, scales = "free_x") +
    scale_fill_manual(values = cb_palette) +
    labs(x = "Soil Water Potential (Median of Bins)", y = "Transpiration Deficit") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white", color = "black"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.position = "top",
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    )
  
  # Save plot if path is provided
  if (!is.null(save_path)) {
    ggsave(filename = save_path, plot = plot, width = 10, height = 6, dpi = 300)
  }
  
  return(plot)
}

plot_box_merge_TDiff_PSI <- function(data, save_path = NULL) {
  library(dplyr)
  library(ggplot2)
  
  # 1. Bin the Soil Water Potential into 100-unit ranges and compute midpoints
  bin_breaks <- seq(floor(min(data$soil_water_potential, na.rm = TRUE)), 
                    ceiling(max(data$soil_water_potential, na.rm = TRUE)), 
                    by = 100)
  
  bin_midpoints <- bin_breaks[-length(bin_breaks)] + diff(bin_breaks) / 2  # Compute midpoints
  
  data <- data %>%
    mutate(soil_bin = cut(soil_water_potential, 
                          breaks = bin_breaks, 
                          include.lowest = TRUE, 
                          right = FALSE, 
                          labels = bin_midpoints))  # Assign midpoints as labels
  
  # Convert to numeric for proper ordering
  data$soil_bin <- as.numeric(as.character(data$soil_bin))
  
  # 2. Reorder Species
  data$species <- factor(data$species, 
                         levels = c("Oak", "Beech", "Spruce", "Pine"))
  
  # 3. Ensure all bins appear for each species
  complete_bins <- expand.grid(soil_bin = unique(data$soil_bin), species = levels(data$species))
  data <- full_join(complete_bins, data, by = c("soil_bin", "species"))
  
  # 4. Plot the Box Plots with color-blind friendly palette
  cb_palette <- c("Oak"   = "#E69F00",  # Orange
                  "Beech" = "#0072B2",  # A deep blue with lower luminance than "#56B4E9"
                  "Spruce"= "#009E73",  # Bluish-green
                  "Pine"  = "#F0E442") 
  
  plot <- ggplot(data, aes(x = factor(soil_bin), y = transpiration_deficit, fill = species)) +
    geom_boxplot(position = position_dodge(width = 0.75), na.rm = TRUE) +
    scale_fill_manual(values = cb_palette) +
    labs(x = "Soil Water Potential (Median of Bins)", y = "Transpiration Deficit", fill = "Species") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white", color = "black"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.position = "top",
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    )
  
  # Save plot if path is provided
  if (!is.null(save_path)) {
    ggsave(filename = save_path, plot = plot, width = 10, height = 6, dpi = 300)
  }
  
  return(plot)
}

plot_mean_box_TDiff_PSI <- function(data, save_path = NULL) {
  library(dplyr)
  library(ggplot2)
  
  # 1. Bin the Soil Water Potential into 100-unit ranges and compute midpoints
  bin_breaks <- seq(floor(min(data$soil_water_potential, na.rm = TRUE)), 
                    ceiling(max(data$soil_water_potential, na.rm = TRUE)), 
                    by = 100)
  
  bin_midpoints <- bin_breaks[-length(bin_breaks)] + diff(bin_breaks) / 2  # Compute midpoints
  
  data <- data %>%
    mutate(soil_bin = cut(soil_water_potential, 
                          breaks = bin_breaks, 
                          include.lowest = TRUE, 
                          right = FALSE, 
                          labels = bin_midpoints))  # Assign midpoints as labels
  
  # Convert to numeric for proper ordering
  data$soil_bin <- as.numeric(as.character(data$soil_bin))
  
  # 2. Reorder Species
  data$species <- factor(data$species, 
                         levels = c("Oak", "Beech", "Spruce", "Pine"))
  
  # 3. Compute mean for each bin and species
  mean_data <- data %>%
    group_by(soil_bin, species) %>%
    summarise(mean_transpiration_deficit = mean(transpiration_deficit, na.rm = TRUE), .groups = 'drop')
  
  # 4. Plot the mean points and lines with color-blind friendly palette
  cb_palette <- c("Oak"   = "#E69F00",  # Orange
                  "Beech" = "#0072B2",  # A deep blue with lower luminance than "#56B4E9"
                  "Spruce"= "#009E73",  # Bluish-green
                  "Pine"  = "#F0E442") 
  
  
  plot <- ggplot(mean_data, aes(x = factor(soil_bin), y = mean_transpiration_deficit, color = species, group = species)) +
    geom_point(size = 3) +
    geom_line() +
    scale_color_manual(values = cb_palette) +
    labs(x = "Soil Water Potential (Median of Bins)", y = "Mean Transpiration Deficit", color = "Species") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white", color = "black"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.position = "top",
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    )
  
  # Save plot if path is provided
  if (!is.null(save_path)) {
    ggsave(filename = save_path, plot = plot, width = 10, height = 6, dpi = 300)
  }
  
  return(plot)
}

plot_mean_box_TDiff_PSI_with_slope <- function(data, save_path = NULL) {
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(broom)  # For extracting regression results
  
  # Bin creation for Soil Water Potential (Panel A)
  bin_breaks <- seq(floor(min(data$soil_water_potential, na.rm = TRUE)), 
                    ceiling(max(data$soil_water_potential, na.rm = TRUE)), 
                    by = 100)
  bin_midpoints <- bin_breaks[-length(bin_breaks)] + diff(bin_breaks) / 2
  
  # Assign soil bin categories
  data <- data %>%
    mutate(soil_bin = cut(soil_water_potential, 
                          breaks = bin_breaks, 
                          include.lowest = TRUE, 
                          right = FALSE, 
                          labels = bin_midpoints))
  
  data$soil_bin <- as.numeric(as.character(data$soil_bin))
  
  data$species <- factor(data$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  
  # Compute mean transpiration deficit for each species within each soil bin (For Panel A)
  mean_data <- data %>%
    group_by(soil_bin, species) %>%
    summarise(mean_transpiration_deficit = mean(transpiration_deficit, na.rm = TRUE), .groups = 'drop')
  
  # Function to calculate slope, R², p-value, and standard error using original transpiration deficit
  calculate_slope_stats <- function(data) {
    regression_results <- data %>%
      group_by(species) %>%
      filter(n() >= 2) %>%
      summarise(
        model = list(lm(transpiration_deficit ~ soil_water_potential, data = cur_data())),
        .groups = 'drop'
      ) %>%
      rowwise() %>%
      mutate(
        slope = coef(model)[2],  # Extract slope
        p_value = summary(model)$coefficients[2, 4],  # Extract p-value
        r_squared = summary(model)$r.squared,  # Extract R²
        std_error = summary(model)$coefficients[2, 2]  # Extract standard error
      ) %>%
      ungroup() %>%
      select(species, slope, p_value, r_squared, std_error) %>%
      arrange(slope)
    
    return(regression_results)
  }
  
  slope_results <- calculate_slope_stats(data)
  
  # Color Palette
  cb_palette <- c("Oak"   = "#E69F00",  # Orange
                  "Beech" = "#0072B2",  # A deep blue with lower luminance than "#56B4E9"
                  "Spruce"= "#009E73",  # Bluish-green
                  "Pine"  = "#F0E442") 
  
  # Scatter plot (Panel A) with binned soil water potential
  plot_a <- ggplot(mean_data, aes(x = factor(soil_bin), y = mean_transpiration_deficit, color = species, group = species)) +
    geom_point(size = 3) +
    geom_line() +
    scale_color_manual(values = cb_palette) +
    labs(x = "Soil Water Potential (Median of 100 KPa Bins)", y = "Mean Transpiration Deficit", color = "Species") +
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
  
  # Bar plot (Panel B) with slopes, p-values, R², and error bars using original transpiration deficit
  plot_b <- ggplot(slope_results, aes(y = reorder(species, slope), x = slope, fill = species)) +
    geom_bar(stat = "identity", width = 0.6, color = "black") +
    geom_errorbar(aes(xmin = slope - std_error, xmax = slope + std_error), width = 0.2) +  # Error bars for slope uncertainty
    scale_fill_manual(values = cb_palette) +
    labs(x = "Slope of Linear Regression (Transpiration Deficit)") +  # Updated y-axis label
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.y = element_blank(),  # Hide y-axis label for Panel B
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
      legend.position = "none",  # Removes legend for Panel B
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    ) +
    geom_text(aes(label = ifelse(p_value < 0.001, 
                                 sprintf("p < 0.001\nR² = %.2f", r_squared), 
                                 sprintf("p = %.3f\nR² = %.2f", p_value, r_squared))), 
              position = position_stack(vjust = 0.5),  # Place text inside bars
              size = 4, color = "black")  # Adjust text size and position
  
  # Combine plots with adjusted widths
  combined_plot <- plot_a + plot_b + plot_annotation(tag_levels = 'A') + plot_layout(widths = c(2, 1))  # Adjusted layout
  
  # Save if needed
  if (!is.null(save_path)) {
    ggsave(filename = save_path, plot = combined_plot, width = 12, height = 6, dpi = 300)
  }
  
  return(combined_plot)
}

plot_mean_box_TDiff_PSI_with_slope_bin <- function(data, save_path = NULL) {
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(broom)  # For extracting regression results
  
  # Bin creation for Soil Water Potential (used for both Panel A and Panel B)
  bin_breaks <- seq(floor(min(data$soil_water_potential, na.rm = TRUE)), 
                    ceiling(max(data$soil_water_potential, na.rm = TRUE)), 
                    by = 100)
  bin_midpoints <- bin_breaks[-length(bin_breaks)] + diff(bin_breaks) / 2
  
  # Assign soil bin categories
  data <- data %>%
    mutate(soil_bin = cut(soil_water_potential, 
                          breaks = bin_breaks, 
                          include.lowest = TRUE, 
                          right = FALSE, 
                          labels = bin_midpoints))
  
  data$soil_bin <- as.numeric(as.character(data$soil_bin))
  
  data$species <- factor(data$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  
  # Compute mean transpiration deficit for each species within each soil bin (For Panel A)
  mean_data <- data %>%
    group_by(soil_bin, species) %>%
    summarise(mean_transpiration_deficit = mean(transpiration_deficit, na.rm = TRUE), .groups = 'drop')
  
  # Function to calculate slope, R², p-value, and standard error using mean transpiration deficit from binned data (For Panel B)
  calculate_slope_stats_bin <- function(mean_data) {
    regression_results <- mean_data %>%
      group_by(species) %>%
      filter(n() >= 2) %>%
      summarise(
        model = list(lm(mean_transpiration_deficit ~ soil_bin, data = cur_data())),  # Regression on binned data
        .groups = 'drop'
      ) %>%
      rowwise() %>%
      mutate(
        slope = coef(model)[2],  # Extract slope
        p_value = summary(model)$coefficients[2, 4],  # Extract p-value
        r_squared = summary(model)$r.squared,  # Extract R²
        std_error = summary(model)$coefficients[2, 2]  # Extract standard error
      ) %>%
      ungroup() %>%
      select(species, slope, p_value, r_squared, std_error) %>%
      arrange(slope)
    
    return(regression_results)
  }
  
  slope_results <- calculate_slope_stats_bin(mean_data)
  
  # Color Palette
  cb_palette <- c("Oak"   = "#E69F00",  # Orange
                  "Beech" = "#0072B2",  # A deep blue with lower luminance than "#56B4E9"
                  "Spruce"= "#009E73",  # Bluish-green
                  "Pine"  = "#F0E442") 
  
  
  # Scatter plot (Panel A) with binned soil water potential
  plot_a <- ggplot(mean_data, aes(x = factor(soil_bin), y = mean_transpiration_deficit, color = species, group = species)) +
    geom_point(size = 3) +
    geom_line() +
    scale_color_manual(values = cb_palette) +
    labs(x = "Soil Water Potential (Median of 100 KPa Bins)", y = "Mean Transpiration Deficit", color = "Species") +
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
  
  # Bar plot (Panel B) with slopes, p-values, R², and error bars using binned soil water potential
  plot_b <- ggplot(slope_results, aes(y = reorder(species, slope), x = slope, fill = species)) +
    geom_bar(stat = "identity", width = 0.6, color = "black") +
    geom_errorbar(aes(xmin = slope - std_error, xmax = slope + std_error), width = 0.2) +  # Error bars for slope uncertainty
    scale_fill_manual(values = cb_palette) +
    labs(x = "Slope of Linear Regression (Binned Transpiration Deficit)") +  # Updated y-axis label
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.y = element_blank(),  # Hide y-axis label for Panel B
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
      legend.position = "none",  # Removes legend for Panel B
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    ) +
    geom_text(aes(label = ifelse(p_value < 0.001, 
                                 sprintf("p < 0.001\nR² = %.2f", r_squared), 
                                 sprintf("p = %.3f\nR² = %.2f", p_value, r_squared))), 
              position = position_stack(vjust = 0.5),  # Place text inside bars
              size = 4, color = "black")  # Adjust text size and position
  
  # Combine plots with adjusted widths
  combined_plot <- plot_a + plot_b + plot_annotation(tag_levels = 'A') + plot_layout(widths = c(1.5, 1))  # Adjusted layout
  
  # Save if needed
  if (!is.null(save_path)) {
    ggsave(filename = save_path, plot = combined_plot, width = 12, height = 6, dpi = 300)
  }
  
  return(combined_plot)
}

plot_mean_box_TDiff_PSI_with_slope_filter <- function(data, save_path = NULL) {
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(broom)  # For extracting regression results
  
  # Bin creation for Soil Water Potential (only for Panel A)
  bin_breaks <- seq(floor(min(data$soil_water_potential, na.rm = TRUE)), 
                    ceiling(max(data$soil_water_potential, na.rm = TRUE)), 
                    by = 100)
  bin_midpoints <- bin_breaks[-length(bin_breaks)] + diff(bin_breaks) / 2
  
  # Assign soil bin categories
  data <- data %>%
    mutate(soil_bin = cut(soil_water_potential, 
                          breaks = bin_breaks, 
                          include.lowest = TRUE, 
                          right = FALSE, 
                          labels = bin_midpoints))
  
  data$soil_bin <- as.numeric(as.character(data$soil_bin))
  
  data$species <- factor(data$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  
  # Compute mean transpiration deficit for each species within each soil bin (For Panel A)
  mean_data <- data %>%
    group_by(soil_bin, species) %>%
    summarise(mean_transpiration_deficit = mean(transpiration_deficit, na.rm = TRUE), .groups = 'drop')
  
  # Function to calculate slope, R², p-value, and standard error using original soil water potential (For Panel B)
  calculate_slope_stats <- function(data) {
    data_filtered <- data %>%
      filter(transpiration_deficit >= 5 & transpiration_deficit <= 15)  # Filter for transpiration deficit 5-10
    
    regression_results <- data_filtered %>%
      group_by(species) %>%
      filter(n() >= 2) %>%
      summarise(
        model = list(lm(transpiration_deficit ~ soil_water_potential, data = cur_data())),  # Use original soil_water_potential
        .groups = 'drop'
      ) %>%
      rowwise() %>%
      mutate(
        slope = abs(coef(model)[2]),  # Extract absolute slope
        p_value = summary(model)$coefficients[2, 4],  # Extract p-value
        r_squared = summary(model)$r.squared,  # Extract R²
        std_error = summary(model)$coefficients[2, 2]  # Extract standard error
      ) %>%
      ungroup() %>%
      select(species, slope, p_value, r_squared, std_error) %>%
      arrange(slope)
    
    return(regression_results)
  }
  
  slope_results <- calculate_slope_stats(data)
  
  # Color Palette
  cb_palette <- c("Oak"   = "#E69F00",  # Orange
                  "Beech" = "#0072B2",  # A deep blue with lower luminance than "#56B4E9"
                  "Spruce"= "#009E73",  # Bluish-green
                  "Pine"  = "#F0E442") 
  
  
  # Scatter plot (Panel A) with binned soil water potential
  plot_a <- ggplot(mean_data, aes(x = factor(soil_bin), y = mean_transpiration_deficit, color = species, group = species)) +
    geom_point(size = 3) +
    geom_line() +
    scale_color_manual(values = cb_palette) +
    labs(x = "Soil Water Potential (Median of 100 KPa Bins)", y = "Mean Transpiration Deficit", color = "Species") +
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
  
  # Bar plot (Panel B) with absolute slopes, p-values, R², and error bars
  plot_b <- ggplot(slope_results, aes(y = reorder(species, slope), x = slope, fill = species)) +
    geom_bar(stat = "identity", width = 0.6, color = "black") +
    geom_errorbar(aes(xmin = pmax(0, slope - std_error),  # Ensure error bars don't go below zero
                      xmax = slope + std_error), width = 0.2) +  # Error bars
    scale_fill_manual(values = cb_palette) +
    labs(x = "Absolute Slope (Transpiration Deficit 5-15)", y = "Species") +  
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.border = element_rect(color = NA, fill = NA, linewidth = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",  # Removes legend for Panel B
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    ) +
    # Display p-values and R² inside the bars, centered
    geom_text(aes(x = slope / 2,  # Position text at the middle of the bar
                  label = ifelse(p_value < 0.001, 
                                 sprintf("p < 0.001\nR² = %.2f", r_squared), 
                                 sprintf("p = %.3f\nR² = %.2f", p_value, r_squared))),
              hjust = 0.5, vjust = 0.5, color = "black", size = 4.5)  # Center text inside bars
  
  # Combine plots with adjusted widths
  combined_plot <- plot_a + plot_b + plot_annotation(tag_levels = 'A') + plot_layout(widths = c(2, 1))  # Adjusted layout
  
  # Save if needed
  if (!is.null(save_path)) {
    ggsave(filename = save_path, plot = combined_plot, width = 12, height = 6, dpi = 300)
  }
  
  return(combined_plot)
}

plot_mean_box_TDiff_PSI_with_slope_bin_filter <- function(data, save_path = NULL) {
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(broom)  # For extracting regression results
  
  # Bin creation for Soil Water Potential (used for both Panel A and Panel B)
  bin_breaks <- seq(floor(min(data$soil_water_potential, na.rm = TRUE)), 
                    ceiling(max(data$soil_water_potential, na.rm = TRUE)), 
                    by = 100)
  bin_midpoints <- bin_breaks[-length(bin_breaks)] + diff(bin_breaks) / 2
  
  # Assign soil bin categories
  data <- data %>%
    mutate(soil_bin = cut(soil_water_potential, 
                          breaks = bin_breaks, 
                          include.lowest = TRUE, 
                          right = FALSE, 
                          labels = bin_midpoints))
  
  data$soil_bin <- as.numeric(as.character(data$soil_bin))
  
  data$species <- factor(data$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  
  # Compute mean transpiration deficit for each species within each soil bin (For Panel A)
  mean_data <- data %>%
    group_by(soil_bin, species) %>%
    summarise(mean_transpiration_deficit = mean(transpiration_deficit, na.rm = TRUE), .groups = 'drop')
  
  # Function to calculate absolute slope, R², p-value, and standard error using binned data (For Panel B)
  calculate_slope_stats_bin <- function(mean_data) {
    data_filtered <- mean_data %>%
      filter(mean_transpiration_deficit >= 0 & mean_transpiration_deficit <= 13)  # Use binned mean transpiration deficit range
    
    regression_results <- data_filtered %>%
      group_by(species) %>%
      filter(n() >= 2) %>%
      summarise(
        model = list(lm(mean_transpiration_deficit ~ soil_bin, data = cur_data())),  # Regression on binned data
        .groups = 'drop'
      ) %>%
      rowwise() %>%
      mutate(
        slope = abs(coef(model)[2]),  # Extract absolute slope
        p_value = summary(model)$coefficients[2, 4],  # Extract p-value
        r_squared = summary(model)$r.squared,  # Extract R²
        std_error = summary(model)$coefficients[2, 2]  # Extract standard error
      ) %>%
      ungroup() %>%
      select(species, slope, p_value, r_squared, std_error) %>%
      arrange(slope)
    
    return(regression_results)
  }
  
  slope_results <- calculate_slope_stats_bin(mean_data)
  
  # Color Palette
  cb_palette <- c("Oak"   = "#E69F00",  # Orange
                  "Beech" = "#0072B2",  # A deep blue with lower luminance than "#56B4E9"
                  "Spruce"= "#009E73",  # Bluish-green
                  "Pine"  = "#F0E442") 
  
  
  # Scatter plot (Panel A) with binned soil water potential
  plot_a <- ggplot(mean_data, aes(x = factor(soil_bin), y = mean_transpiration_deficit, color = species, group = species)) +
    geom_point(size = 3) +
    geom_line() +
    scale_color_manual(values = cb_palette) +
    labs(x = "Soil Water Potential (Median of 100 KPa Bins)", y = "Mean Transpiration Deficit", color = "Species") +
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
  
  # Bar plot (Panel B) with absolute slopes, p-values, R², and error bars
  plot_b <- ggplot(slope_results, aes(y = reorder(species, slope), x = slope, fill = species)) +
    geom_bar(stat = "identity", width = 0.6, color = "black") +
    geom_errorbar(aes(xmin = pmax(0, slope - std_error),  # Ensure error bars don't go below zero
                      xmax = slope + std_error), width = 0.2) +  # Error bars
    scale_fill_manual(values = cb_palette) +
    labs(x = "Absolute Slope (Transpiration Deficit 0-13)", y = "Species") +  # Updated x-axis label
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
      legend.position = "none",  # Removes legend for Panel B
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    ) +
    # Display p-values and R² inside the bars, centered
    geom_text(aes(x = slope / 2,  # Position text at the middle of the bar
                  label = ifelse(p_value < 0.001, 
                                 sprintf("p < 0.001\nR² = %.2f", r_squared), 
                                 sprintf("p = %.3f\nR² = %.2f", p_value, r_squared))),
              hjust = 0.5, vjust = 0.5, color = "black", size = 4.5)  # Center text inside bars
  
  # Combine plots with adjusted widths
  combined_plot <- plot_a + plot_b + plot_annotation(tag_levels = 'A') + plot_layout(widths = c(1.5, 1))  # Adjusted layout
  
  # Save if needed
  if (!is.null(save_path)) {
    ggsave(filename = save_path, plot = combined_plot, width = 12, height = 6, dpi = 300)
  }
  
  return(combined_plot)
}

plot_box_NDVI_TDiff <- function(data, save_path = NULL) {
  library(dplyr)
  library(ggplot2)
  
  # Identify the correct column (Quantiles or Proportions)
  value_col <- if ("Quantiles" %in% names(data)) "Quantiles" else "Proportions"
  
  # 1. Bin the Transpiration Deficit into 5-unit ranges and compute midpoints
  bin_breaks <- seq(floor(min(data$transpiration_deficit, na.rm = TRUE)), 
                    ceiling(max(data$transpiration_deficit, na.rm = TRUE)), 
                    by = 5)
  
  bin_midpoints <- bin_breaks[-length(bin_breaks)] + diff(bin_breaks) / 2  # Compute midpoints
  
  data <- data %>%
    mutate(transp_bin = cut(transpiration_deficit, 
                            breaks = bin_breaks, 
                            include.lowest = TRUE, 
                            right = FALSE, 
                            labels = bin_midpoints))  # Assign midpoints as labels
  
  # Convert to numeric for proper ordering
  data$transp_bin <- as.numeric(as.character(data$transp_bin))
  
  # Remove last bin if NA
  data <- data %>% filter(!is.na(transp_bin))
  
  # 2. Reorder Species
  data$species <- factor(data$species, 
                         levels = c("Oak", "Beech", "Spruce", "Pine"))
  
  # 3. Plot the Box Plots with color-blind friendly palette
  cb_palette <- c("Oak"   = "#E69F00",  # Orange
                  "Beech" = "#0072B2",  # A deep blue with lower luminance than "#56B4E9"
                  "Spruce"= "#009E73",  # Bluish-green
                  "Pine"  = "#F0E442") 
  
  
  plot <- ggplot(data, aes(x = factor(transp_bin), y = .data[[value_col]], fill = species)) +
    geom_boxplot() +
    facet_wrap(~ species, scales = "free_x") +
    scale_fill_manual(values = cb_palette) +
    labs(x = "Transpiration Deficit (Median of Bins)", y = value_col) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white", color = "black"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.position = "top",
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    )
  
  # Save plot if path is provided
  if (!is.null(save_path)) {
    ggsave(filename = save_path, plot = plot, width = 10, height = 6, dpi = 300)
  }
  
  return(plot)
}

plot_box_merge_NDVI_TDiff <- function(data, save_path = NULL) {
  library(dplyr)
  library(ggplot2)
  
  # Identify the correct column (Quantiles or Proportions)
  value_col <- if ("Quantiles" %in% names(data)) "Quantiles" else "Proportions"
  
  # 1. Bin the Transpiration Deficit into 5-unit ranges and compute midpoints
  bin_breaks <- seq(floor(min(data$transpiration_deficit, na.rm = TRUE)), 
                    ceiling(max(data$transpiration_deficit, na.rm = TRUE)), 
                    by = 5)
  
  bin_midpoints <- bin_breaks[-length(bin_breaks)] + diff(bin_breaks) / 2  # Compute midpoints
  
  data <- data %>%
    mutate(transp_bin = cut(transpiration_deficit, 
                            breaks = bin_breaks, 
                            include.lowest = TRUE, 
                            right = FALSE, 
                            labels = bin_midpoints))  # Assign midpoints as labels
  
  # Convert to numeric for proper ordering
  data$transp_bin <- as.numeric(as.character(data$transp_bin))
  
  # Remove last bin if NA
  data <- data %>% filter(!is.na(transp_bin))
  
  # 2. Reorder Species
  data$species <- factor(data$species, 
                         levels = c("Oak", "Beech", "Spruce", "Pine"))
  
  # 3. Ensure all bins appear for each species
  complete_bins <- expand.grid(transp_bin = unique(data$transp_bin), species = levels(data$species))
  data <- full_join(complete_bins, data, by = c("transp_bin", "species"))
  
  # 4. Plot the Box Plots with color-blind friendly palette
  cb_palette <- c("Oak"   = "#E69F00",  # Orange
                  "Beech" = "#0072B2",  # A deep blue with lower luminance than "#56B4E9"
                  "Spruce"= "#009E73",  # Bluish-green
                  "Pine"  = "#F0E442") 
  
  
  plot <- ggplot(data, aes(x = factor(transp_bin), y = .data[[value_col]], fill = species)) +
    geom_boxplot(position = position_dodge(width = 0.75), na.rm = TRUE) +
    scale_fill_manual(values = cb_palette) +
    labs(x = "Transpiration Deficit (Median of Bins)", y = value_col, fill = "Species") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white", color = "black"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.position = "top",
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    )
  
  # Save plot if path is provided
  if (!is.null(save_path)) {
    ggsave(filename = save_path, plot = plot, width = 10, height = 6, dpi = 300)
  }
  
  return(plot)
}

plot_mean_box_NDVI_TDiff <- function(data, save_path = NULL) {
  library(dplyr)
  library(ggplot2)
  
  # Identify the correct column (Quantiles or Proportions)
  value_col <- if ("Quantiles" %in% names(data)) "Quantiles" else "Proportions"
  
  # 1. Bin the Transpiration Deficit into 5-unit ranges and compute midpoints
  bin_breaks <- seq(floor(min(data$transpiration_deficit, na.rm = TRUE)), 
                    ceiling(max(data$transpiration_deficit, na.rm = TRUE)), 
                    by = 5)
  
  bin_midpoints <- bin_breaks[-length(bin_breaks)] + diff(bin_breaks) / 2  # Compute midpoints
  
  data <- data %>%
    mutate(transp_bin = cut(transpiration_deficit, 
                            breaks = bin_breaks, 
                            include.lowest = TRUE, 
                            right = FALSE, 
                            labels = bin_midpoints))  # Assign midpoints as labels
  
  # Convert to numeric for proper ordering
  data$transp_bin <- as.numeric(as.character(data$transp_bin))
  
  # Remove last bin if NA
  data <- data %>% filter(!is.na(transp_bin))
  
  # 2. Reorder Species
  data$species <- factor(data$species, 
                         levels = c("Oak", "Beech", "Spruce", "Pine"))
  
  # 3. Compute mean for each bin and species
  mean_data <- data %>%
    group_by(transp_bin, species) %>%
    summarise(mean_value = mean(.data[[value_col]], na.rm = TRUE), .groups = 'drop')
  
  # 4. Plot the mean points and lines with color-blind friendly palette
  cb_palette <- c("Oak"   = "#E69F00",  # Orange
                  "Beech" = "#0072B2",  # A deep blue with lower luminance than "#56B4E9"
                  "Spruce"= "#009E73",  # Bluish-green
                  "Pine"  = "#F0E442") 
  
  
  plot <- ggplot(mean_data, aes(x = factor(transp_bin), y = mean_value, color = species, group = species)) +
    geom_point(size = 3) +
    geom_line() +
    scale_color_manual(values = cb_palette) +
    labs(x = "Transpiration Deficit (Median of Bins)", y = paste("Mean", value_col), color = "Species") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white", color = "black"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.position = "top",
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    )
  
  # Save plot if path is provided
  if (!is.null(save_path)) {
    ggsave(filename = save_path, plot = plot, width = 10, height = 6, dpi = 300)
  }
  
  return(plot)
}

plot_mean_box_NDVI_TDiff_with_slope <- function(data, save_path = NULL) {
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  
  # Identify the correct column (Quantiles or Proportions)
  value_col <- if ("Quantiles" %in% names(data)) "Quantiles" else "Proportions"
  
  # Bin creation for Transpiration Deficit (only for Panel A)
  bin_breaks <- seq(floor(min(data$transpiration_deficit, na.rm = TRUE)), 
                    ceiling(max(data$transpiration_deficit, na.rm = TRUE)), 
                    by = 5)
  bin_midpoints <- bin_breaks[-length(bin_breaks)] + diff(bin_breaks) / 2
  
  # Assign transpiration bin categories
  data <- data %>%
    mutate(transp_bin = cut(transpiration_deficit, 
                            breaks = bin_breaks, 
                            include.lowest = TRUE, 
                            right = FALSE, 
                            labels = bin_midpoints))
  
  data$transp_bin <- as.numeric(as.character(data$transp_bin))
  
  # Remove last bin if NA
  data <- data %>% filter(!is.na(transp_bin))
  
  data$species <- factor(data$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  
  # Compute mean value (Quantiles or Proportions) for each species within each transpiration bin (For Panel A)
  mean_data <- data %>%
    group_by(transp_bin, species) %>%
    summarise(mean_value = mean(.data[[value_col]], na.rm = TRUE), .groups = 'drop')
  
  # Function to calculate slope, R², p-value, and standard error using original transpiration deficit (For Panel B)
  calculate_slope_stats <- function(data) {
    data_filtered <- data %>%
      filter(transp_bin < 20)  
    
    regression_results <- data_filtered %>%
      group_by(species) %>%
      filter(n() >= 2) %>%
      summarise(
        model = list(lm(.data[[value_col]] ~ transpiration_deficit, data = cur_data())),  # Use original transpiration_deficit
        .groups = 'drop'
      ) %>%
      rowwise() %>%
      mutate(
        slope = coef(model)[2],  # Extract slope
        p_value = summary(model)$coefficients[2, 4],  # Extract p-value
        r_squared = summary(model)$r.squared,  # Extract R²
        std_error = summary(model)$coefficients[2, 2]  # Extract standard error
      ) %>%
      ungroup() %>%
      select(species, slope, p_value, r_squared, std_error) %>%
      arrange(slope)
    
    return(regression_results)
  }
  
  slope_results <- calculate_slope_stats(data)
  
  # Color Palette
  cb_palette <- c("Oak"   = "#E69F00",  # Orange
                  "Beech" = "#0072B2",  # A deep blue with lower luminance than "#56B4E9"
                  "Spruce"= "#009E73",  # Bluish-green
                  "Pine"  = "#F0E442") 
  
  
  # Scatter plot (Panel A) with binned transpiration deficit
  plot_a <- ggplot(mean_data, aes(x = factor(transp_bin), y = mean_value, color = species, group = species)) +
    geom_point(size = 3) +
    geom_line() +
    scale_color_manual(values = cb_palette) +
    labs(x = "Transpiration Deficit (Median of 5-unit Bins)", 
         y = paste("Mean", value_col), 
         color = "Species") +
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
  
  # Bar plot (Panel B) with slopes, p-values, R², and error bars using original transpiration deficit
  plot_b <- ggplot(slope_results, aes(y = reorder(species, slope), x = slope, fill = species)) +
    geom_bar(stat = "identity", width = 0.6, color = "black") +
    geom_errorbar(aes(xmin = slope - std_error, xmax = slope + std_error), width = 0.2) +
    scale_fill_manual(values = cb_palette) +
    labs(x = "Slope (TDiff < 20)") +
    theme_minimal() +
    scale_x_reverse() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.y = element_blank(),  # Hide y-axis label for Panel B
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
      legend.position = "none",  # Removes legend for Panel B
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    ) +
    geom_text(aes(label = ifelse(p_value < 0.001, 
                                 sprintf("p < 0.001\nR² = %.2f", r_squared), 
                                 sprintf("p = %.3f\nR² = %.2f", p_value, r_squared))), 
              position = position_stack(vjust = 0.5),  # Place text inside bars
              size = 4, color = "black")  # Adjust text size and position
  
  combined_plot <- plot_a + plot_b + plot_annotation(tag_levels = 'A') + plot_layout(widths = c(2, 1))
  
  if (!is.null(save_path)) {
    ggsave(filename = save_path, plot = combined_plot, width = 12, height = 6, dpi = 300)
  }
  
  return(combined_plot)
}

plot_mean_box_NDVI_TDiff_with_slope_bin <- function(data, save_path = NULL) {
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  
  # Identify the correct column (Quantiles or Proportions)
  value_col <- if ("Quantiles" %in% names(data)) "Quantiles" else "Proportions"
  
  # Bin creation for Transpiration Deficit (used for both Panel A and Panel B)
  bin_breaks <- seq(floor(min(data$transpiration_deficit, na.rm = TRUE)), 
                    ceiling(max(data$transpiration_deficit, na.rm = TRUE)), 
                    by = 5)
  bin_midpoints <- bin_breaks[-length(bin_breaks)] + diff(bin_breaks) / 2
  
  # Assign transpiration bin categories
  data <- data %>%
    mutate(transp_bin = cut(transpiration_deficit, 
                            breaks = bin_breaks, 
                            include.lowest = TRUE, 
                            right = FALSE, 
                            labels = bin_midpoints))
  
  data$transp_bin <- as.numeric(as.character(data$transp_bin))
  
  # Remove last bin if NA
  data <- data %>% filter(!is.na(transp_bin))
  
  data$species <- factor(data$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  
  # Compute mean value (Quantiles or Proportions) for each species within each transpiration bin (For Panel A)
  mean_data <- data %>%
    group_by(transp_bin, species) %>%
    summarise(mean_value = mean(.data[[value_col]], na.rm = TRUE), .groups = 'drop')
  
  # Function to calculate slope, R², p-value, and standard error using mean values from binned data (For Panel B)
  calculate_slope_stats_bin <- function(mean_data) {
    data_filtered <- mean_data %>%
      filter(transp_bin < 20)
    
    regression_results <- data_filtered %>%
      group_by(species) %>%
      filter(n() >= 2) %>%
      summarise(
        model = list(lm(mean_value ~ transp_bin, data = cur_data())),  # Regression on binned data
        .groups = 'drop'
      ) %>%
      rowwise() %>%
      mutate(
        slope = coef(model)[2],  # Extract slope
        p_value = summary(model)$coefficients[2, 4],  # Extract p-value
        r_squared = summary(model)$r.squared,  # Extract R²
        std_error = summary(model)$coefficients[2, 2]  # Extract standard error
      ) %>%
      ungroup() %>%
      select(species, slope, p_value, r_squared, std_error) %>%
      arrange(slope)
    
    return(regression_results)
  }
  
  slope_results <- calculate_slope_stats_bin(mean_data)
  
  # Color Palette
  cb_palette <- c("Oak"   = "#E69F00",  # Orange
                  "Beech" = "#0072B2",  # A deep blue with lower luminance than "#56B4E9"
                  "Spruce"= "#009E73",  # Bluish-green
                  "Pine"  = "#F0E442") 
  
  
  # Scatter plot (Panel A) with binned transpiration deficit
  plot_a <- ggplot(mean_data, aes(x = factor(transp_bin), y = mean_value, color = species, group = species)) +
    geom_point(size = 3) +
    geom_line() +
    scale_color_manual(values = cb_palette) +
    labs(x = "Transpiration Deficit (Median of 5-unit Bins)", 
         y = paste("Mean", value_col), 
         color = "Species") +
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
  
  # Bar plot (Panel B) with slopes, p-values, R², and error bars using binned transpiration deficit
  plot_b <- ggplot(slope_results, aes(y = reorder(species, slope), x = slope, fill = species)) +
    geom_bar(stat = "identity", width = 0.6, color = "black") +
    geom_errorbar(aes(xmin = slope - std_error, xmax = slope + std_error), width = 0.2) +
    scale_fill_manual(values = cb_palette) +
    labs(x = "Slope (TDiff < 20)") +
    theme_minimal() +
    scale_x_reverse() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.y = element_blank(),  # Hide y-axis label for Panel B
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
      legend.position = "none",  # Removes legend for Panel B
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    ) +
    geom_text(aes(label = sprintf("p = %.3f\nR² = %.2f", p_value, r_squared)), 
              position = position_stack(vjust = 0.5),
              size = 4, color = "black")
  
  combined_plot <- plot_a + plot_b + plot_annotation(tag_levels = 'A') + plot_layout(widths = c(1.5, 1))
  
  if (!is.null(save_path)) {
    ggsave(filename = save_path, plot = combined_plot, width = 12, height = 6, dpi = 300)
  }
  
  return(combined_plot)
}

NDVI_TDiffbin <- function(df, bin_width = 5) {
  
  library(dplyr)
  
  # Identify the correct column (Quantiles or Proportions)
  value_col <- if ("Quantiles" %in% names(df)) "Quantiles" else "Proportions"
  
  # Compute total pixels per species
  species_totals <- df %>%
    group_by(species) %>%
    summarise(total_pixels = n(), .groups = "drop")
  
  # Define bin breaks for transpiration_deficit
  tdiff_min <- floor(min(df$transpiration_deficit, na.rm = TRUE))
  tdiff_max <- ceiling(max(df$transpiration_deficit, na.rm = TRUE))
  bin_breaks <- seq(tdiff_min, tdiff_max, by = bin_width)
  
  # Create TDiff_bin column
  df <- df %>%
    mutate(TDiff_bin = cut(transpiration_deficit,
                           breaks = bin_breaks,
                           include.lowest = TRUE,
                           right = FALSE))
  
  # Helper function to compute bin median
  get_bin_median <- function(bin_label) {
    nums <- as.numeric(strsplit(gsub("\\[|\\]|\\(|\\)", "", bin_label), ",")[[1]])
    mean(nums)
  }
  
  # Group, summarize, calculate percentage, and filter
  meanNDVI_TDiffbin_species <- df %>%
    group_by(species, TDiff_bin) %>%
    summarise(
      avg_value = mean(.data[[value_col]], na.rm = TRUE),
      count = n(),
      .groups = 'drop'
    ) %>%
    mutate(bin_median = sapply(as.character(TDiff_bin), get_bin_median)) %>%
    left_join(species_totals, by = "species") %>%
    mutate(percentage = round(count / total_pixels * 100, 2)) %>%
    filter(percentage >= 0.1, count > 2000) %>%
    select(species, TDiff_bin, bin_median, avg_value, count, total_pixels, percentage)
  
  return(meanNDVI_TDiffbin_species)
}

plot_NDVI_TDiff_poly_2_slope <- function(species_data, output_file) {
  # Load required libraries
  library(lme4)      # For mixed-effects modeling
  library(dplyr)     # For data manipulation
  library(ggplot2)   # For plotting
  library(patchwork) # For combining ggplots
  library(tidyr)     # For pivoting data
  library(tibble)    # For rownames_to_column
  
  # Identify the correct value column (Quantiles or Proportions)
  value_col <- if ("Quantiles" %in% names(species_data)) "Quantiles" else "Proportions"
  
  # Create the binned dataset and remove NA rows.
  NDVI_TDiffbin_df <- NDVI_TDiffbin(species_data)
  NDVI_TDiffbin_df <- na.omit(NDVI_TDiffbin_df)
  
  # Define custom color palette and species order.
  cb_palette <- c("Oak"   = "#E69F00",  # Orange
                  "Beech" = "#0072B2",  # Deep blue
                  "Spruce"= "#009E73",  # Bluish-green
                  "Pine"  = "#F0E442")  # Yellow
  
  species_levels <- c("Oak", "Beech", "Spruce", "Pine")
  NDVI_TDiffbin_df$species <- factor(NDVI_TDiffbin_df$species, levels = species_levels)
  
  # Set fixed threshold value to 11.5.
  threshold <- 11.5
  
  ## Panel A: Mixed-Effects Model Plot for avg_value vs bin_median
  # Fit a raw degree-2 polynomial mixed model:
  # NDVI = a + b * bin_median + c * bin_median^2 + random effects by species.
  model <- lmer(avg_value ~ poly(bin_median, 2, raw = TRUE) + 
                  (poly(bin_median, 2, raw = TRUE) | species),
                data = NDVI_TDiffbin_df)
  
  # Create prediction data for each species.
  pred_data <- NDVI_TDiffbin_df %>%
    group_by(species) %>%
    summarise(min_bin = min(bin_median),
              max_bin = max(bin_median),
              .groups = "drop") %>%
    group_by(species) %>%
    do(data.frame(bin_median = seq(.$min_bin, .$max_bin, length.out = 100))) %>%
    ungroup()
  
  # Ensure species factor order.
  pred_data$species <- factor(pred_data$species, levels = species_levels)
  
  # Generate predicted values (including random effects).
  pred_data$predicted <- predict(model, newdata = pred_data, re.form = NULL)
  
  # Create Panel A plot with fixed threshold line and annotation.
  plot_mixed <- ggplot(NDVI_TDiffbin_df, aes(x = bin_median, y = avg_value, color = species)) +
    geom_point() +
    geom_line(data = pred_data, aes(x = bin_median, y = predicted, color = species), size = 1) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "black", size = 1) +
    # Place the "median" label at the right-hand side of the hline.
    annotate("text", 
             x = max(NDVI_TDiffbin_df$bin_median, na.rm = TRUE), 
             y = threshold, 
             label = "median", 
             hjust = 1.1, vjust = -0.3, 
             fontface = "italic", size = 4) +
    scale_color_manual(values = cb_palette) +
    labs(x = "Transpiration Deficit (bin_median)",
         y = paste("Average", value_col)) +
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
  
  ## Panel B: Bar Plot of Local Slopes at threshold (11.5)
  # For each species, compute the local slope using the predicted values 
  # around the bin_median where the prediction is closest to the threshold.
  local_slope_data <- pred_data %>%
    group_by(species) %>%
    do({
      df <- .
      thresh_val <- threshold
      idx <- which.min(abs(df$predicted - thresh_val))
      win_start <- max(1, idx - 5)
      win_end <- min(nrow(df), idx + 5)
      df_window <- df[win_start:win_end, ]
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
  
  local_slope_data <- local_slope_data %>%
    mutate(slope_abs = abs(slope),
           label_text = ifelse(p_value < 0.05,
                               paste0("p<0.05\nR²: ", round(r2, 3)),
                               paste0("Slope: ", round(slope, 4), "\n",
                                      "p: ", signif(p_value, 3), "\n",
                                      "R²: ", round(r2, 3))))
  
  p_bar <- ggplot(local_slope_data, aes(x = reorder(species, slope_abs), y = slope_abs, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = slope_abs - slope_se, ymax = slope_abs + slope_se), width = 0.2, color = "black") +
    geom_text(aes(label = label_text), position = position_stack(vjust = 0.5), size = 4, color = "black") +
    scale_fill_manual(values = cb_palette) +
    labs(x = "Species",
         y = paste("Absolute Slope at", threshold)) +
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
  
  final_plot <- plot_mixed + p_bar + plot_layout(widths = c(2, 1))
  
  # Save the slope figure to file.
  ggsave(output_file, plot = final_plot, width = 10, height = 8, dpi = 300)
  
  # Return the slope plot.
  return(final_plot)
}

plot_NDVI_TDiff_poly_3_slope <- function(species_data, output_file) {
  # Load required libraries
  library(lme4)      # For mixed-effects modeling
  library(dplyr)     # For data manipulation
  library(ggplot2)   # For plotting
  library(patchwork) # For combining ggplots
  
  # Identify the correct value column (Quantiles or Proportions)
  value_col <- if ("Quantiles" %in% names(species_data)) "Quantiles" else "Proportions"
  
  # Create the binned dataset and remove NA rows.
  NDVI_TDiffbin_df <- NDVI_TDiffbin(species_data)
  NDVI_TDiffbin_df <- na.omit(NDVI_TDiffbin_df)
  
  # Define custom color palette and species order.
  cb_palette <- c("Oak"   = "#E69F00",  # Orange
                  "Beech" = "#0072B2",  # Deep blue
                  "Spruce"= "#009E73",  # Bluish-green
                  "Pine"  = "#F0E442")  # Yellow
  
  species_levels <- c("Oak", "Beech", "Spruce", "Pine")
  NDVI_TDiffbin_df$species <- factor(NDVI_TDiffbin_df$species, levels = species_levels)
  
  # Set fixed threshold value to 11.5.
  threshold <- 11.5
  
  ## Panel A: Mixed-Effects Model Plot for avg_value vs bin_median
  # Fit a cubic (degree-3) mixed-effects model:
  model <- lmer(avg_value ~ poly(bin_median, 3) + (poly(bin_median, 3) | species),
                data = NDVI_TDiffbin_df)
  
  # Create prediction data for each species.
  pred_data <- NDVI_TDiffbin_df %>%
    group_by(species) %>%
    summarise(min_bin = min(bin_median),
              max_bin = max(bin_median),
              .groups = "drop") %>%
    group_by(species) %>%
    do(data.frame(bin_median = seq(.$min_bin, .$max_bin, length.out = 100))) %>%
    ungroup()
  
  # Ensure species factor order.
  pred_data$species <- factor(pred_data$species, levels = species_levels)
  
  # Generate predicted values (including random effects).
  pred_data$predicted <- predict(model, newdata = pred_data, re.form = NULL)
  
  # Create Panel A plot with fixed threshold line and annotation.
  plot_mixed <- ggplot(NDVI_TDiffbin_df, aes(x = bin_median, y = avg_value, color = species)) +
    geom_point() +
    geom_line(data = pred_data, aes(x = bin_median, y = predicted, color = species), size = 1) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "black", size = 1) +
    # Place the "median" label at the right-hand side of the hline.
    annotate("text", 
             x = max(NDVI_TDiffbin_df$bin_median, na.rm = TRUE), 
             y = threshold, 
             label = "median", 
             hjust = 1.1, vjust = -0.3, 
             fontface = "italic", size = 4) +
    scale_color_manual(values = cb_palette) +
    labs(x = "Transpiration Deficit (bin_median)",
         y = paste("Average", value_col)) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    )
  
  ## Panel B: Bar Plot of Local Slopes at Threshold = 11.5
  # For each species, compute the local slope using the predicted values 
  # around the bin_median where the prediction is closest to the fixed threshold.
  local_slope_data <- pred_data %>%
    group_by(species) %>%
    do({
      df <- .
      thresh_val <- threshold  # fixed threshold
      idx <- which.min(abs(df$predicted - thresh_val))
      win_start <- max(1, idx - 5)
      win_end <- min(nrow(df), idx + 5)
      df_window <- df[win_start:win_end, ]
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
  
  local_slope_data <- local_slope_data %>%
    mutate(slope_abs = abs(slope),
           label_text = ifelse(p_value < 0.05,
                               paste0("p<0.05\nR²: ", round(r2, 3)),
                               paste0("Slope: ", round(slope, 4), "\n",
                                      "p: ", signif(p_value, 3), "\n",
                                      "R²: ", round(r2, 3))))
  
  p_bar <- ggplot(local_slope_data, aes(x = reorder(species, slope_abs), y = slope_abs, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = slope_abs - slope_se, ymax = slope_abs + slope_se), width = 0.2, color = "black") +
    geom_text(aes(label = label_text), position = position_stack(vjust = 0.5), size = 4, color = "black") +
    scale_fill_manual(values = cb_palette) +
    labs(x = "Species",
         y = paste("Absolute Slope at", threshold)) +
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
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
  
  final_plot <- plot_mixed + p_bar + plot_layout(widths = c(2, 1))
  
  print(final_plot)
  ggsave(output_file, plot = final_plot, width = 10, height = 8, dpi = 300)
  
  return(final_plot)
}

plot_Quantiles_TDiff_linear_slope_coeff <- function(data, coef_fig, output_figure) {
  
  # -------------------------------
  # SETUP: Load required libraries
  # -------------------------------
  require(ggplot2)
  require(nlme)
  require(dplyr)
  require(tibble)
  require(patchwork)
  require(purrr)
  require(car)      # for deltaMethod if needed
  require(broom)
  require(tidyr)
  
  # -------------------------------
  # DATA PREPARATION
  # -------------------------------
  # Use the NDVI values stored in "avg_value"
  value_col <- "avg_value"
  
  # Order species and define a color palette (ordered as Oak, Beech, Spruce, Pine)
  data$species <- factor(data$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  cb_palette <- c("Oak"   = "#E69F00",  # Orange
                  "Beech" = "#0072B2",  # Deep blue
                  "Spruce"= "#009E73",  # Bluish-green
                  "Pine"  = "#F0E442")  # Yellow
  
  # Prepare the data using your custom function (which creates bin_median, etc.)
  data <- NDVI_TDiffbin(data)
  
  # For Transpiration Deficit, use x = bin_median (assumed to be positive)
  data <- data %>% mutate(x = bin_median)
  
  # Clean data: remove rows with missing or non-finite values in avg_value and x
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  data_clean$species <- factor(data_clean$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  
  # Set the NDVI threshold at which we calculate x50 and the slope
  threshold <- 11.5
  
  # -------------------------------
  # FITTING THE LINEAR MODEL PER SPECIES (nlsList)
  # -------------------------------
  start_list_linear <- list(a = 1, b = 0.1)
  nls_models_linear <- nlsList(
    avg_value ~ a + b * x | species,
    data = data_clean,
    start = start_list_linear,
    control = nls.control()
  )
  print(summary(nls_models_linear))
  
  # -------------------------------
  # EXTRACT TIDY COEFFICIENTS FOR THE COEFFICIENT BAR PLOT (Panel: Coef)
  # -------------------------------
  tidy_linear <- map_df(nls_models_linear, tidy, .id = "species")
  tidy_linear$species <- factor(tidy_linear$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  tidy_linear <- tidy_linear %>%
    mutate(label = if_else(p.value < 0.05,
                           "*",
                           sprintf("p=%.2f", p.value)))
  
  # -------------------------------
  # COEFFICIENT BAR PLOT (Panel: Coef)
  # -------------------------------
  p_coeff <- ggplot(tidy_linear, aes(x = species, y = estimate, fill = species)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(aes(label = label, y = estimate/2),
              position = position_dodge(width = 0.9),
              vjust = 0.5, color = "black", size = 4) +
    scale_fill_manual(values = cb_palette) +
    facet_wrap(~ term, scales = "free_y") +
    labs(title = "Linear Coefficients by Species for NDVI",
         subtitle = expression(NDVI == a + b * x),
         x = "Species",
         y = "Coefficient Value") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.background = element_rect(fill = "white", color = "white"),
          panel.background = element_rect(fill = "white"),
          legend.background = element_rect(fill = "white", color = "white"),
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
          plot.subtitle = element_text(hjust = 0.5, size = 14, face = "italic", color = "black"),
          axis.title = element_text(face = "bold"),
          axis.text = element_text(color = "black"),
          panel.border = element_rect(color = NA, fill = NA, linewidth = 0),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "top",
          strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
          strip.text = element_text(face = "bold", size = 12))
  
  print(p_coeff)
  
  # Save the coefficient plot to file (separate figure)
  ggsave(filename = coef_fig,
         plot = p_coeff, device = "png", width = 8, height = 6, dpi = 300)
  
  # -------------------------------
  # PANEL A: Observed Data and Fitted Curves
  # -------------------------------
  # For the linear model, we plot the fitted line over the entire range of x.
  pred_list <- data_clean %>%
    group_by(species) %>%
    do({
      sp <- unique(.$species)
      x_seq <- seq(min(.$x, na.rm = TRUE), max(.$x, na.rm = TRUE), length.out = 100)
      sp_model <- nls_models_linear[[as.character(sp)]]
      pred <- predict(sp_model, newdata = data.frame(x = x_seq))
      data.frame(species = sp, x = x_seq, pred = pred)
    })
  pred_all <- bind_rows(pred_list) %>% mutate(pred = as.numeric(pred))
  
  p_combined <- ggplot() +
    geom_point(data = data_clean, aes(x = x, y = avg_value, color = species)) +
    geom_line(data = pred_all, aes(x = x, y = pred, color = species, group = species), size = 1) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "black", size = 1) +
    annotate("text", x = min(data_clean$x, na.rm = TRUE), y = threshold,
             label = "median", hjust = -5, vjust = -0.3, fontface = "italic", size = 4) +
    scale_color_manual(values = cb_palette) +
    labs(x = "Transpiration Deficit", y = "NDVI") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.background = element_rect(fill = "white", color = "white"),
          panel.background = element_rect(fill = "white"),
          legend.background = element_rect(fill = "white", color = "white"),
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
          axis.title = element_text(face = "bold"),
          axis.text = element_text(color = "black"),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "top") +
    labs(caption = "(a)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0),
          plot.caption.position = "plot") +
    coord_cartesian(clip = "off")
  
  # -------------------------------
  # Define calc_x50 Function for the Linear Model
  # -------------------------------
  # For a linear model: threshold = a + b*x  =>  x50 = (threshold - a) / b
  calc_x50_linear <- function(a, b, threshold, x_range) {
    x50_val <- (threshold - a) / b
    return(x50_val)
  }
  
  # -------------------------------
  # CALCULATE x50 AND SLOPE FOR EACH SPECIES
  # -------------------------------
  stats_list <- list()
  species_levels <- levels(data_clean$species)
  
  for (sp in species_levels) {
    mod <- nls_models_linear[[as.character(sp)]]
    coefs <- coef(mod)  # a and b
    sp_data <- data_clean %>% filter(species == sp)
    x_range <- c(min(sp_data$x, na.rm = TRUE), max(sp_data$x, na.rm = TRUE))
    x50_val <- calc_x50_linear(coefs["a"], coefs["b"], threshold, x_range)
    
    # In a linear model, the slope is simply b.
    slope_val <- coefs["b"]
    vcov_mat <- vcov(mod)
    slope_se <- sqrt(vcov_mat["b", "b"])
    
    # Compute p-value for the slope using a z statistic approximation
    z_val <- slope_val / slope_se
    slope_p <- 2 * (1 - pnorm(abs(z_val)))
    
    stats_list[[sp]] <- data.frame(species = sp, x50 = x50_val, slope_abs = abs(slope_val), se = slope_se, slope_p = slope_p)
  }
  
  stats_df <- do.call(rbind, stats_list)
  stats_df$species <- factor(stats_df$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  
  print(stats_df %>% select(species, x50, slope_abs, se, slope_p))
  
  # -------------------------------
  # COMPUTE R² FOR EACH SPECIES
  # -------------------------------
  rsq_df <- map_df(levels(data_clean$species), function(sp) {
    sp_data <- dplyr::filter(data_clean, species == sp)
    mod <- nls_models_linear[[sp]]
    pred <- predict(mod, newdata = sp_data)
    r2_val <- 1 - sum((sp_data[[value_col]] - pred)^2) / sum((sp_data[[value_col]] - mean(sp_data[[value_col]]))^2)
    data.frame(species = sp, r2 = r2_val)
  })
  
  # Join R² values with stats_df
  stats_df <- left_join(stats_df, rsq_df, by = "species")
  
  # Create label for Panel C based on p-value condition:
  # If slope_p < 0.05: label = "R²=...*"
  # Otherwise: label = "p=...\nR²=..."
  stats_df <- stats_df %>%
    mutate(slope_label = if_else(slope_p < 0.05,
                                 sprintf("%.2f*", r2),
                                 sprintf("p=%.2f\nR²=%.2f", slope_p, r2)))
  
  # -------------------------------
  # PANEL B: Bar Plot of x50 Values (Transpiration Deficit)
  # -------------------------------
  p_x50 <- ggplot(stats_df, aes(x = species, y = x50, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    labs(x = "Species", y = "Transpiration Deficit (x50)") +
    scale_fill_manual(values = cb_palette) +
    scale_x_discrete(limits = c("Oak", "Beech", "Spruce", "Pine")) +
    expand_limits(y = 0) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          axis.title.x = element_blank(),
          axis.title = element_text(face = "bold"),
          axis.text = element_text(color = "black"),
          plot.background = element_rect(fill = "white", color = "white"),
          panel.background = element_rect(fill = "white"),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(caption = "(b)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0),
          plot.caption.position = "plot")
  
  # -------------------------------
  # PANEL C: Bar Plot of Absolute Slope at x50 with Error Bars and Conditional Labels
  # -------------------------------
  p_slope <- ggplot(stats_df, aes(x = species, y = slope_abs, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, slope_abs - se), ymax = slope_abs + se), width = 0.2) +
    geom_text(aes(y = slope_abs/2, label = slope_label), color = "black", size = 4) +
    labs(x = "Species", y = "Absolute Slope") +
    scale_fill_manual(values = cb_palette) +
    scale_x_discrete(limits = c("Oak", "Beech", "Spruce", "Pine")) +
    theme_minimal() +
    labs(caption = "(c)") +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          axis.title = element_text(face = "bold"),
          axis.text = element_text(color = "black"),
          plot.background = element_rect(fill = "white", color = "white"),
          panel.background = element_rect(fill = "white"),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.caption = element_text(face = "bold", size = 16, hjust = 0),
          plot.caption.position = "plot")
  
  # -------------------------------
  # COMBINE PANELS AND SAVE FINAL PLOT
  # -------------------------------
  final_plot <- p_combined + (p_x50 / p_slope) + plot_layout(widths = c(2, 1))
  print(final_plot)
  
  ggsave(filename = output_figure,
         plot = final_plot, device = "png", width = 10, height = 8, dpi = 300)
}

plot_Quantiles_TDiff_exp_slope_coeff <- function(data, coef_fig, output_figure) {
  # Process data using the custom NDVI_TDiffbin function and remove missing values
  data <- NDVI_TDiffbin(data, 5)
  data <- na.omit(data)
  
  # Load required libraries (assumed to be installed)
  library(ggplot2)
  library(nlme)
  library(dplyr)
  library(tibble)
  library(patchwork)
  library(purrr)
  library(car)      # for deltaMethod
  library(broom)
  library(tidyr)
  
  # Identify the value column and set species order and color palette
  value_col <- "avg_value"
  data$species <- factor(data$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  cb_palette <- c("Oak"   = "#E69F00",  # Orange
                  "Beech" = "#0072B2",  # Deep blue
                  "Spruce"= "#009E73",  # Bluish-green
                  "Pine"  = "#F0E442")  # Yellow
  
  # For Transpiration Deficit, use x = bin_median directly (assumed positive)
  data <- data %>% mutate(x = bin_median)
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  
  # Set threshold using the median of avg_value
  # threshold <- median(data_clean$avg_value, na.rm = TRUE)
  threshold <- 11.5
  
  # NONLINEAR MODELING PER SPECIES
  start_list <- list(a = 5, b = 7, c = 0.04)
  control_params <- nls.control(maxiter = 1200, minFactor = 1e-09)
  
  nls_models <- nlsList(avg_value ~ a + b * exp(-c * x) | species,
                        data = data_clean,
                        start = start_list,
                        control = control_params)
  print(summary(nls_models))
  
  # Extract coefficient summaries for a, b, and c
  model_summary <- summary(nls_models)
  coeff_a <- model_summary$coefficients[, , "a"]
  coeff_b <- model_summary$coefficients[, , "b"]
  coeff_c <- model_summary$coefficients[, , "c"]
  
  # Create a data frame with species names, estimates, and p-values
  results_df <- data.frame(
    species    = rownames(coeff_a),
    a_estimate = coeff_a[, "Estimate"],
    a_pvalue   = coeff_a[, "Pr(>|t|)"],
    b_estimate = coeff_b[, "Estimate"],
    b_pvalue   = coeff_b[, "Pr(>|t|)"],
    c_estimate = coeff_c[, "Estimate"],
    c_pvalue   = coeff_c[, "Pr(>|t|)"]
  )
  results_df$species <- factor(results_df$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  
  # Reshape the results into long format: first for estimates, then merge p-values
  coeffs_long <- results_df %>%
    pivot_longer(
      cols = c(a_estimate, b_estimate, c_estimate),
      names_to = "Coefficient",
      values_to = "Value",
      names_pattern = "([abc])_.*"
    ) %>%
    left_join(
      results_df %>%
        pivot_longer(
          cols = c(a_pvalue, b_pvalue, c_pvalue),
          names_to = "Coefficient",
          values_to = "pvalue",
          names_pattern = "([abc])_.*"
        ),
      by = c("species", "Coefficient")
    ) %>%
    mutate(
      label = if_else(pvalue < 0.05, "*", sprintf("%.3f", pvalue))
    )
  
  # Create the faceted coefficient bar plot with p-value labels
  p_coeff <- ggplot(coeffs_long, aes(x = species, y = Value, fill = species)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(aes(label = label, y = Value/2),
              position = position_dodge(width = 0.9),
              vjust = 0.5, color = "black", size = 4) +
    scale_fill_manual(values = cb_palette) +
    facet_wrap(~ Coefficient, scales = "free_y") +
    labs(title = "Coefficients by Species for NDVI Quantiles",
         x = "Species",
         y = "Coefficient Value") +
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
  print(p_coeff)
  
  # Save the coefficient plot
  ggsave(filename = coef_fig, plot = p_coeff, device = "png", width = 8, height = 6, dpi = 300)
  
  # -------------------------------
  # Additional analysis and combined panels
  # -------------------------------
  
  # Extract coefficients for a, b, and c per species
  coef_df <- as.data.frame(coef(nls_models), optional = TRUE) %>% 
    rownames_to_column(var = "species") %>%
    filter(!is.na(a))
  coef_df$species <- factor(coef_df$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  
  # Calculate x50: x50 = -log((threshold - a)/b) / c
  coef_df <- coef_df %>%
    mutate(x50 = ifelse((threshold - a) > 0 & b > 0,
                        -log((threshold - a)/b) / c,
                        NA))
  
  # Calculate slope at x50: slope50 = -c * (threshold - a)
  coef_df <- coef_df %>% mutate(slope50 = -c * (threshold - a))
  
  # Generate predictions for each species (Panel A)
  pred_list <- data_clean %>%
    group_by(species) %>%
    do({
      sp <- unique(.$species)
      x_seq <- seq(min(.$x, na.rm = TRUE), max(.$x, na.rm = TRUE), length.out = 100)
      sp_model <- nls_models[[as.character(sp)]]
      pred <- predict(sp_model, newdata = data.frame(x = x_seq))
      data.frame(x = x_seq, pred = pred)
    })
  pred_all <- bind_rows(pred_list)
  
  p_combined <- ggplot() +
    geom_point(data = data_clean, aes(x = x, y = avg_value, color = species)) +
    geom_line(data = pred_all, aes(x = x, y = pred, color = species), size = 1) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "black", size = 1) +
    annotate("text", x = min(data_clean$x, na.rm = TRUE), y = threshold,
             label = "median", hjust = -2.1, vjust = -0.3, fontface = "italic", size = 4) +
    scale_color_manual(values = cb_palette) +
    labs(x = "Transpiration Deficit", y = "NDVI Quantiles") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    ) +
    labs(caption = "(a)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0),
          plot.caption.position = "plot") +
    coord_cartesian(clip = "off")
  
  # Panel B: Bar plot of x50 values
  p_x50 <- ggplot(coef_df, aes(x = species, y = x50, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    labs(x = "Species", y = "Transpiration Deficit") +
    scale_fill_manual(values = cb_palette) +
    expand_limits(y = 0) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(caption = "(b)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0),
          plot.caption.position = "plot")
  
  # Panel C: Bar plot of absolute slope at x50 with error bars and R² annotations
  stats_list <- lapply(levels(data_clean$species), function(sp) {
    mod <- nls_models[[sp]]
    coefs <- coef(mod)
    # Compute slope50 = -c * (threshold - a)
    slope50_est <- -coefs["c"] * (threshold - coefs["a"])
    # Compute standard error using deltaMethod for f(a, c) = c*(a - threshold)
    dm_result <- deltaMethod(mod, paste0("c*(a - ", threshold, ")"),
                             parameterNames = c("a", "b", "c"))
    se <- dm_result$SE
    # Compute degrees of freedom from the model summary
    df_mod <- summary(mod)$df[2]
    t_val <- slope50_est / se
    p_val <- 2 * (1 - pt(abs(t_val), df_mod))
    
    # Compute R² for the species model
    df_sp <- data_clean %>% filter(species == sp)
    fitted_vals <- predict(mod, newdata = df_sp)
    r_squared <- 1 - sum((df_sp[[value_col]] - fitted_vals)^2) /
      sum((df_sp[[value_col]] - mean(df_sp[[value_col]]))^2)
    
    tibble(species = sp,
           slope50 = slope50_est,
           slope_abs = abs(slope50_est),
           se = se,
           p_val = p_val,
           r_squared = r_squared)
  })
  stats_df <- bind_rows(stats_list)
  stats_df$species <- factor(stats_df$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  
  stats_df <- stats_df %>%
    mutate(label_text = sprintf("%.2f%s", r_squared, ifelse(p_val < 0.05, "*", "")))
  
  p_slope <- ggplot(stats_df, aes(x = species, y = slope_abs, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, slope_abs - se), ymax = slope_abs + se), width = 0.2) +
    geom_text(aes(y = slope_abs/2, label = label_text), color = "black", size = 4) +
    labs(x = "Species", y = "Absolute Slope") +
    scale_fill_manual(values = cb_palette) +
    expand_limits(y = 0) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(caption = "(c)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0),
          plot.caption.position = "plot")
  
  # Combine panels: Panel (a) on the left; Panel (b) (top) and Panel (c) (bottom) in right column
  final_plot <- p_combined + (p_x50 / p_slope) + plot_layout(widths = c(2, 1))
  print(final_plot)
  
  # Save the combined final plot
  ggsave(filename = output_figure, plot = final_plot, device = "png", width = 10, height = 8, dpi = 300)
}

plot_Quantiles_TDiff_poly_2_slope_coeff <- function(data, coef_fig, output_figure) {
  
  # -------------------------------
  # SETUP: Load required libraries
  # -------------------------------
  require(ggplot2)
  require(nlme)
  require(dplyr)
  require(tibble)
  require(patchwork)
  require(purrr)
  require(car)      # for deltaMethod
  require(broom)
  require(tidyr)
  
  # -------------------------------
  # DATA PREPARATION
  # -------------------------------
  # Use the NDVI values stored in "avg_value"
  value_col <- "avg_value"
  
  # Order species and define a color palette (ordered as Oak, Beech, Spruce, Pine)
  data$species <- factor(data$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  cb_palette <- c("Oak"   = "#E69F00",  # Orange
                  "Beech" = "#0072B2",  # Deep blue
                  "Spruce"= "#009E73",  # Bluish-green
                  "Pine"  = "#F0E442")  # Yellow
  
  # Prepare the data using your custom function (which creates bin_median, etc.)
  data <- NDVI_TDiffbin(data)
  
  # For Transpiration Deficit, use x = bin_median (assumed to be positive)
  data <- data %>% mutate(x = bin_median)
  
  # Clean data: remove rows with missing or non-finite values in avg_value and x
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  # Ensure species order is maintained in the cleaned data as well
  data_clean$species <- factor(data_clean$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  
  # Set the threshold (NDVI value) at which we calculate x50 and the slope
  threshold <- 11.5
  
  # -------------------------------
  # FITTING THE QUADRATIC MODEL PER SPECIES (nlsList)
  # -------------------------------
  start_list_poly <- list(a = 1, b = 0.1, c = 0.1)
  nls_models_poly <- nlsList(
    avg_value ~ a + b * x + c * I(x^2) | species,
    data = data_clean,
    start = start_list_poly,
    control = nls.control()
  )
  print(summary(nls_models_poly))
  
  # -------------------------------
  # EXTRACT TIDY COEFFICIENTS FOR THE COEFFICIENT BAR PLOT (Panel: Coef)
  # -------------------------------
  tidy_poly <- map_df(nls_models_poly, tidy, .id = "species")
  tidy_poly$species <- factor(tidy_poly$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  tidy_poly <- tidy_poly %>%
    mutate(label = if_else(p.value < 0.05,
                           "*",
                           sprintf("%.2f", p.value)))
  
  # -------------------------------
  # COEFFICIENT BAR PLOT (Panel: Coef)
  # -------------------------------
  p_coeff <- ggplot(tidy_poly, aes(x = species, y = estimate, fill = species)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(aes(label = label, y = estimate/2),
              position = position_dodge(width = 0.9),
              vjust = 0.5, color = "black", size = 4) +
    scale_fill_manual(values = cb_palette) +
    facet_wrap(~ term, scales = "free_y") +
    labs(title = "Quadratic Coefficients by Species for NDVI",
         subtitle = expression(NDVI == a + b * x + c * x^2),
         x = "Species",
         y = "Coefficient Value") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.background = element_rect(fill = "white", color = "white"),
          panel.background = element_rect(fill = "white"),
          legend.background = element_rect(fill = "white", color = "white"),
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
          plot.subtitle = element_text(hjust = 0.5, size = 14, face = "italic", color = "black"),
          axis.title = element_text(face = "bold"),
          axis.text = element_text(color = "black"),
          panel.border = element_rect(color = NA, fill = NA, linewidth = 0),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "top",
          strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
          strip.text = element_text(face = "bold", size = 12))
  
  print(p_coeff)
  
  # Save the coefficient plot to file (separate figure)
  ggsave(filename = coef_fig,
         plot = p_coeff, device = "png", width = 8, height = 6, dpi = 300)
  
  # -------------------------------
  # PANEL A: Observed Data and Fitted Curves (plot only left of the vertex)
  # -------------------------------
  pred_list <- data_clean %>%
    group_by(species) %>%
    do({
      sp <- unique(.$species)
      x_seq <- seq(min(.$x, na.rm = TRUE), max(.$x, na.rm = TRUE), length.out = 100)
      sp_model <- nls_models_poly[[as.character(sp)]]
      pred <- predict(sp_model, newdata = data.frame(x = x_seq))
      data.frame(species = sp, x = x_seq, pred = pred)
    })
  pred_all <- bind_rows(pred_list) %>% mutate(pred = as.numeric(pred))
  
  # Calculate the vertex (x-coordinate of the minimum) for each quadratic curve
  coef_df_vertex <- tidy_poly %>%
    select(species, term, estimate) %>%
    pivot_wider(names_from = term, values_from = estimate) %>%
    mutate(vertex = -b / (2 * c))
  
  pred_all <- pred_all %>%
    left_join(coef_df_vertex %>% select(species, vertex), by = "species")
  pred_all_left <- pred_all %>% filter(x <= vertex)
  
  p_combined <- ggplot() +
    geom_point(data = data_clean, aes(x = x, y = avg_value, color = species)) +
    geom_line(data = pred_all_left, aes(x = x, y = pred, color = species, group = species), size = 1) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "black", size = 1) +
    annotate("text", x = min(data_clean$x, na.rm = TRUE), y = threshold,
             label = "median", hjust = -2.1, vjust = -0.3, fontface = "italic", size = 4) +
    scale_color_manual(values = cb_palette) +
    labs(x = "Transpiration Deficit", y = "NDVI") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.background = element_rect(fill = "white", color = "white"),
          panel.background = element_rect(fill = "white"),
          legend.background = element_rect(fill = "white", color = "white"),
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
          axis.title = element_text(face = "bold"),
          axis.text = element_text(color = "black"),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "top") +
    labs(caption = "(a)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0),
          plot.caption.position = "plot") +
    coord_cartesian(clip = "off")
  
  # -------------------------------
  # Define calc_x50 Function
  # -------------------------------
  calc_x50 <- function(a, b, c, threshold, x_range) {
    discriminant <- b^2 - 4 * c * (a - threshold)
    if (discriminant < 0) return(NA)
    sol1 <- (-b + sqrt(discriminant)) / (2 * c)
    sol2 <- (-b - sqrt(discriminant)) / (2 * c)
    valid_solutions <- c(sol1, sol2)[c(sol1, sol2) >= x_range[1] & c(sol1, sol2) <= x_range[2]]
    if (length(valid_solutions) == 0) {
      all_sols <- c(sol1, sol2)
      diff <- abs(all_sols - mean(x_range))
      return(all_sols[which.min(diff)])
    } else if (length(valid_solutions) == 1) {
      return(valid_solutions)
    } else {
      return(valid_solutions[which.min(abs(valid_solutions - median(x_range)))])
    }
  }
  
  # -------------------------------
  # CALCULATE x50 AND SLOPE AT x50 FOR EACH SPECIES
  # -------------------------------
  stats_list <- list()
  species_levels <- levels(data_clean$species)
  
  for (sp in species_levels) {
    mod <- nls_models_poly[[as.character(sp)]]
    coefs <- coef(mod)  # a, b, c
    sp_data <- data_clean %>% filter(species == sp)
    x_range <- c(min(sp_data$x, na.rm = TRUE), max(sp_data$x, na.rm = TRUE))
    x50_val <- calc_x50(coefs["a"], coefs["b"], coefs["c"], threshold, x_range)
    
    disc <- coefs["b"]^2 - 4 * coefs["c"] * (coefs["a"] - threshold)
    slope_val <- if (disc < 0) NA else sqrt(disc)
    
    expr <- paste0("sqrt(b^2 - 4*c*(a - ", threshold, "))")
    vcov_mat <- vcov(mod)
    dm <- deltaMethod(coefs, expr, vcov_mat)
    slope_se <- as.numeric(dm["SE"])
    
    # Compute p-value for slope using a z statistic approximation
    z_val <- slope_val / slope_se
    slope_p <- 2 * (1 - pnorm(abs(z_val)))
    
    stats_list[[sp]] <- data.frame(species = sp, x50 = x50_val, slope_abs = slope_val, se = slope_se, slope_p = slope_p)
  }
  
  stats_df <- do.call(rbind, stats_list)
  # Ensure the species order in the stats data is consistent
  stats_df$species <- factor(stats_df$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  
  print(stats_df %>% select(species, x50, slope_abs, se, slope_p))
  
  # -------------------------------
  # COMPUTE R² FOR EACH SPECIES
  # -------------------------------
  rsq_df <- map_df(levels(data_clean$species), function(sp) {
    sp_data <- dplyr::filter(data_clean, species == sp)
    mod <- nls_models_poly[[sp]]
    pred <- predict(mod, newdata = sp_data)
    r2_val <- 1 - sum((sp_data[[value_col]] - pred)^2) / sum((sp_data[[value_col]] - mean(sp_data[[value_col]]))^2)
    data.frame(species = sp, r2 = r2_val)
  })
  
  # Join R² values with stats_df
  stats_df <- left_join(stats_df, rsq_df, by = "species")
  
  # Create label for Panel C based on p-value condition:
  # If slope_p < 0.05: label = "R²=...*"
  # Otherwise: label = "p=...\nR²=..."
  stats_df <- stats_df %>%
    mutate(slope_label = if_else(slope_p < 0.05,
                                 sprintf("%.2f*", r2),
                                 sprintf("p=%.2f\nR²=%.2f", slope_p, r2)))
  
  # -------------------------------
  # PANEL B: Bar Plot of x50 Values (Transpiration Deficit)
  # -------------------------------
  p_x50 <- ggplot(stats_df, aes(x = species, y = x50, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    labs(x = "Species", y = "Transpiration Deficit") +
    scale_fill_manual(values = cb_palette) +
    scale_x_discrete(limits = c("Oak", "Beech", "Spruce", "Pine")) +
    expand_limits(y = 0) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          axis.title.x = element_blank(),
          axis.title = element_text(face = "bold"),
          axis.text = element_text(color = "black"),
          plot.background = element_rect(fill = "white", color = "white"),
          panel.background = element_rect(fill = "white"),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(caption = "(b)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0),
          plot.caption.position = "plot")
  
  # -------------------------------
  # PANEL C: Bar Plot of Absolute Slope at x50 with Error Bars and Conditional Labels
  # -------------------------------
  p_slope <- ggplot(stats_df, aes(x = species, y = slope_abs, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, slope_abs - se), ymax = slope_abs + se), width = 0.2) +
    geom_text(aes(y = slope_abs/2, label = slope_label), color = "black", size = 4) +
    labs(x = "Species", y = "Absolute Slope") +
    scale_fill_manual(values = cb_palette) +
    scale_x_discrete(limits = c("Oak", "Beech", "Spruce", "Pine")) +
    theme_minimal() +
    labs(caption = "(c)") +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          axis.title = element_text(face = "bold"),
          axis.text = element_text(color = "black"),
          plot.background = element_rect(fill = "white", color = "white"),
          panel.background = element_rect(fill = "white"),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.caption = element_text(face = "bold", size = 16, hjust = 0),
          plot.caption.position = "plot")
  
  # -------------------------------
  # COMBINE PANELS AND SAVE FINAL PLOT
  # -------------------------------
  final_plot <- p_combined + (p_x50 / p_slope) + plot_layout(widths = c(2, 1))
  print(final_plot)
  
  ggsave(filename = output_figure,
         plot = final_plot, device = "png", width = 10, height = 8, dpi = 300)
}

plot_Quantiles_TDiff_poly_3_slope_coeff <- function(data, coef_fig, output_figure) {
  require(ggplot2)
  require(nlme)
  require(dplyr)
  require(tibble)
  require(patchwork)
  require(purrr)
  require(car)      # for deltaMethod
  require(broom)
  require(tidyr)
  
  value_col <- "avg_value"
  # Order species as "Oak", "Beech", "Spruce", "Pine"
  data$species <- factor(data$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  cb_palette <- c("Oak" = "#E69F00", "Beech" = "#0072B2", "Spruce" = "#009E73", "Pine" = "#F0E442")
  
  # Prepare data using your custom NDVI_TDiffbin function and set x = bin_median
  data <- NDVI_TDiffbin(data)
  data <- data %>% mutate(x = bin_median)
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  data_clean$species <- factor(data_clean$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  
  threshold <- 11.5
  
  # Fit cubic model using poly(x, 3, raw=TRUE) so that the model is: a + b*x + c*x^2 + d*x^3
  start_list_poly <- list(a = 1, b = 0.1, c = 0.1, d = 0.01)
  nls_models_poly <- nlsList(
    avg_value ~ a + b * poly(x, 3, raw = TRUE)[,1] + c * poly(x, 3, raw = TRUE)[,2] + d * poly(x, 3, raw = TRUE)[,3] | species,
    data = data_clean,
    start = start_list_poly,
    control = nls.control()
  )
  print(summary(nls_models_poly))
  
  tidy_poly <- map_df(nls_models_poly, tidy, .id = "species")
  tidy_poly$species <- factor(tidy_poly$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  tidy_poly <- tidy_poly %>% mutate(label = if_else(p.value < 0.05, "*", sprintf("%.2f", p.value)))
  
  p_coeff <- ggplot(tidy_poly, aes(x = species, y = estimate, fill = species)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(aes(label = label, y = estimate/2),
              position = position_dodge(width = 0.9),
              vjust = 0.5, color = "black", size = 4) +
    scale_fill_manual(values = cb_palette) +
    facet_wrap(~ term, scales = "free_y") +
    labs(title = "Cubic Coefficients by Species for NDVI",
         subtitle = expression(NDVI == a + b*x + c*x^2 + d*x^3),
         x = "Species", y = "Coefficient Value") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.background = element_rect(fill = "white", color = "white"),
          panel.background = element_rect(fill = "white"),
          legend.background = element_rect(fill = "white", color = "white"),
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 14, face = "italic"),
          axis.title = element_text(face = "bold"),
          axis.text = element_text(color = "black"),
          panel.grid = element_blank(),
          legend.position = "top",
          strip.background = element_rect(fill = "white", color = "black"),
          strip.text = element_text(face = "bold", size = 12))
  
  print(p_coeff)
  ggsave(filename = coef_fig, plot = p_coeff, device = "png", width = 8, height = 6, dpi = 300)
  
  # PANEL A: Observed Data and Fitted Curves
  pred_list <- data_clean %>%
    group_by(species) %>%
    do({
      sp <- unique(.$species)
      x_seq <- seq(min(.$x, na.rm = TRUE), max(.$x, na.rm = TRUE), length.out = 100)
      sp_model <- nls_models_poly[[as.character(sp)]]
      pred <- predict(sp_model, newdata = data.frame(x = x_seq))
      data.frame(species = sp, x = x_seq, pred = pred)
    })
  pred_all <- bind_rows(pred_list) %>% mutate(pred = as.numeric(pred))
  
  p_combined <- ggplot() +
    geom_point(data = data_clean, aes(x = x, y = avg_value, color = species)) +
    geom_line(data = pred_all, aes(x = x, y = pred, color = species), size = 1) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "black", size = 1) +
    annotate("text", x = min(data_clean$x, na.rm = TRUE), y = threshold, label = "median",
             hjust = -2.1, vjust = -0.3, fontface = "italic", size = 4) +
    scale_color_manual(values = cb_palette) +
    labs(x = "Transpiration Deficit", y = "NDVI", caption = "(a)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.background = element_rect(fill = "white", color = "white"),
          panel.background = element_rect(fill = "white"),
          legend.background = element_rect(fill = "white", color = "white"),
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
          axis.title = element_text(face = "bold"),
          axis.text = element_text(color = "black"),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "top") +
    labs(caption = "(a)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0),
          plot.caption.position = "plot") +
    coord_cartesian(clip = "off")
  
  calc_x50_exact <- function(a, b, c, d, threshold, x_range) {
    # Set up polynomial coefficients for: d*x^3 + c*x^2 + b*x + (a - threshold) = 0
    coeffs <- c(d, c, b, a - threshold)
    roots <- polyroot(coeffs)
    # Keep only roots with negligible imaginary part
    real_roots <- Re(roots[abs(Im(roots)) < 1e-6])
    
    if (length(real_roots) == 0) {
      return(NA)
    }
    
    # Find roots that lie within the observed x range
    valid_roots <- real_roots[real_roots >= x_range[1] & real_roots <= x_range[2]]
    
    if (length(valid_roots) > 0) {
      # If multiple valid roots exist, choose the one closest to the median
      x50 <- valid_roots[which.min(abs(valid_roots - median(x_range)))]
    } else {
      # Otherwise, choose the real root closest to the median of x_range
      x50 <- real_roots[which.min(abs(real_roots - median(x_range)))]
    }
    
    return(x50)
  }
  
  
  stats_df <- do.call(rbind, stats_list)
  
  print(stats_df %>% select(species, x50, slope_abs, se, slope_p))
  
  rsq_df <- map_df(levels(data_clean$species), function(sp) {
    sp_data <- filter(data_clean, species == sp)
    mod <- nls_models_poly[[sp]]
    pred <- predict(mod, newdata = sp_data)
    r2_val <- 1 - sum((sp_data[[value_col]] - pred)^2) / sum((sp_data[[value_col]] - mean(sp_data[[value_col]]))^2)
    data.frame(species = sp, r2 = r2_val)
  })
  
  stats_df <- left_join(stats_df, rsq_df, by = "species")
  stats_df <- stats_df %>%
    mutate(slope_label = if_else(slope_p < 0.05,
                                 sprintf("%.2f*", r2),
                                 sprintf("p=%.2f\nR²=%.2f", slope_p, r2)))
  
  stats_df$species <- factor(stats_df$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  
  # PANEL B: Bar Plot of x50 Values
  p_x50 <- ggplot(stats_df, aes(x = species, y = x50, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = cb_palette) +
    labs(x = "Species", y = "Transpiration Deficit", caption = "(b)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          axis.title.x = element_blank(),
          axis.title = element_text(face = "bold"),
          axis.text = element_text(color = "black"),
          plot.background = element_rect(fill = "white", color = "white"),
          panel.background = element_rect(fill = "white"),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(caption = "(b)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0),
          plot.caption.position = "plot")
  
  # PANEL C: Bar Plot of Absolute Slope at x50 with Error Bars and Labels
  p_slope <- ggplot(stats_df, aes(x = species, y = slope_abs, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, slope_abs - se), ymax = slope_abs + se), width = 0.2) +
    geom_text(aes(y = slope_abs / 2, label = slope_label), color = "black", size = 4) +
    labs(x = "Species", y = "Absolute Slope", caption = "(c)") +
    scale_fill_manual(values = cb_palette) +
    theme_minimal() +
    labs(caption = "(c)") +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          axis.title = element_text(face = "bold"),
          axis.text = element_text(color = "black"),
          plot.background = element_rect(fill = "white", color = "white"),
          panel.background = element_rect(fill = "white"),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.caption = element_text(face = "bold", size = 16, hjust = 0),
          plot.caption.position = "plot")
  
  final_plot <- p_combined + (p_x50 / p_slope) + plot_layout(widths = c(2, 1))
  print(final_plot)
  ggsave(filename = output_figure,
         plot = final_plot, device = "png", width = 10, height = 8, dpi = 300)
}

plot_TDiff_PSIbin_poly_3_slope <- function(data, coef_output, figure_output) {
  # Load required libraries
  library(lme4)      # For mixed-effects modeling
  library(dplyr)     # For data manipulation
  library(ggplot2)   # For plotting
  library(gridExtra) # For arranging multiple plots
  library(patchwork) # For combining ggplots
  library(tidyr)     # For pivoting data
  
  # Process the data and remove any rows with missing values.
  TDiff_PSIbin_df <- TDiff_PSIbin(data)
  TDiff_PSIbin_df <- na.omit(TDiff_PSIbin_df)
  
  # Define a custom color palette and set species order.
  cb_palette <- c("Oak"   = "#E69F00",   # Orange
                  "Beech" = "#0072B2",   # Deep blue
                  "Spruce"= "#009E73",   # Bluish-green
                  "Pine"  = "#F0E442")   # Yellow
  species_levels <- c("Oak", "Beech", "Spruce", "Pine")
  TDiff_PSIbin_df$species <- factor(TDiff_PSIbin_df$species, levels = species_levels)
  
  ## Panel A: Mixed-Effects Model Plot
  
  # Fit the mixed-effects model with a third-order polynomial for bin_median.
  model <- lmer(avg_transpiration_deficit ~ poly(bin_median, 3) + 
                  (poly(bin_median, 3) | species),
                data = TDiff_PSIbin_df)
  
  # Check for singularity and warn if the model is singular.
  if(isSingular(model, tol = 1e-4)) {
    warning("Model is singular (boundary fit). Consider simplifying the random effects structure.")
  }
  
  # Create prediction data for each species.
  pred_data <- TDiff_PSIbin_df %>%
    group_by(species) %>%
    summarise(min_bin = min(bin_median),
              max_bin = max(bin_median)) %>%
    group_by(species) %>%
    do(data.frame(bin_median = seq(.$min_bin, .$max_bin, length.out = 100))) %>%
    ungroup()
  pred_data$species <- factor(pred_data$species, levels = species_levels)
  
  # Predict the fitted values using the model (including random effects).
  pred_data$predicted <- predict(model, newdata = pred_data, re.form = NULL)
  
  # Calculate the overall mean and median of transpiration deficit.
  line_val <- mean(data$transpiration_deficit)
  line_median <- median(data$transpiration_deficit)
  
  # Create the mixed-effects model plot with horizontal lines for mean and median.
  plot_mixed <- ggplot(TDiff_PSIbin_df, aes(x = bin_median, y = avg_transpiration_deficit, color = species)) +
    geom_point() +
    geom_line(data = pred_data, aes(x = bin_median, y = predicted, color = species), linewidth = 1) +
    geom_hline(yintercept = line_val, linetype = "dashed", color = "black", linewidth = 1) +
    annotate("text", x = min(TDiff_PSIbin_df$bin_median, na.rm = TRUE), 
             y = line_val, 
             label = paste0("mean: ", round(line_val, 2)),
             hjust = -0.1, vjust = -0.3, fontface = "italic", size = 4) +
    geom_hline(yintercept = line_median, linetype = "dashed", color = "black", linewidth = 1) +
    annotate("text", x = min(TDiff_PSIbin_df$bin_median, na.rm = TRUE), 
             y = line_median, 
             label = paste0("median: ", round(line_median, 2)),
             hjust = -0.1, vjust = -0.3, fontface = "italic", size = 4) +
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
  
  ## Panel B: Bar Plot of x-values at Mean Transpiration Deficit
  
  x_at_mean <- pred_data %>%
    group_by(species) %>%
    summarise(x_at_mean = bin_median[which.min(abs(predicted - line_val))])
  x_at_mean$species <- factor(x_at_mean$species, levels = species_levels)
  
  p_bar_x <- ggplot(x_at_mean, aes(x = species, y = x_at_mean, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = cb_palette) +
    labs(x = "", y = "PSI at Mean TDiff") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_text(face = "bold"),
          axis.text = element_text(color = "black"),
          plot.background = element_rect(fill = "white", color = "white"),
          panel.background = element_rect(fill = "white"),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  ## Panel C: Bar Plot of Local Slopes at 30% of Maximum Transpiration Deficit
  
  local_slope_data <- pred_data %>%
    group_by(species) %>%
    do({
      df <- .
      max_pred <- max(df$predicted, na.rm = TRUE)
      thresh_val <- 0.3 * max_pred
      idx <- which.min(abs(df$predicted - thresh_val))
      win_start <- max(1, idx - 5)
      win_end <- min(nrow(df), idx + 5)
      df_window <- df[win_start:win_end, ]
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
  
  local_slope_data <- local_slope_data %>%
    mutate(slope_abs = abs(slope),
           label_text = ifelse(p_value < 0.05,
                               paste0(round(r2, 3), " *"),
                               paste0("p = ", signif(p_value, 3))))
  
  p_bar <- ggplot(local_slope_data, aes(x = species, y = slope_abs, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = slope_abs - slope_se, ymax = slope_abs + slope_se),
                  width = 0.2, color = "black") +
    geom_text(aes(y = slope_abs/2, label = label_text), size = 4, color = "black") +
    scale_fill_manual(values = cb_palette) +
    labs(x = "", y = "Absolute Slope at Mean PSI") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_text(face = "bold"),
          axis.text = element_text(color = "black"),
          plot.background = element_rect(fill = "white", color = "white"),
          panel.background = element_rect(fill = "white"),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  ## Combine Panels: Panel A on the left, Panels B (top) and C (bottom) on the right.
  right_panel <- p_bar_x / p_bar
  final_plot <- plot_mixed + right_panel + plot_layout(widths = c(2, 1))
  print(final_plot)
  
  # Ensure that the directories for the output files exist.
  if (!dir.exists(dirname(figure_output))) {
    dir.create(dirname(figure_output), recursive = TRUE)
  }
  if (!dir.exists(dirname(coef_output))) {
    dir.create(dirname(coef_output), recursive = TRUE)
  }
  
  # Save the combined plot.
  ggsave(figure_output, plot = final_plot, width = 10, height = 8, dpi = 300)
  
  ## New Figure: Grouped Bar Plot for All Model Coefficients per Species
  
  # Extract species-specific (conditional) coefficients.
  species_coef <- coef(model)$species
  
  # Convert rownames to a proper column and pivot the data to long format.
  coeff_data <- species_coef %>%
    tibble::rownames_to_column("species") %>%
    pivot_longer(cols = -species, names_to = "term", values_to = "value") %>%
    mutate(species = factor(species, levels = species_levels),
           term = case_when(
             term == "(Intercept)" ~ "Intercept",
             term == "poly(bin_median, 3)1" ~ "Linear",
             term == "poly(bin_median, 3)2" ~ "Quadratic",
             term == "poly(bin_median, 3)3" ~ "Cubic",
             TRUE ~ term
           ))
  
  # Set the order of the term factor.
  coeff_data$term <- factor(coeff_data$term, levels = c("Intercept", "Linear", "Quadratic", "Cubic"))
  
  ## Compute p-values for coefficients by fitting separate linear models for each species.
  species_list <- species_levels
  coeff_stats_list <- lapply(species_list, function(sp) {
    subdata <- subset(TDiff_PSIbin_df, species == sp)
    mod_sp <- lm(avg_transpiration_deficit ~ poly(bin_median, 3), data = subdata)
    summ <- summary(mod_sp)$coefficients
    df <- as.data.frame(summ)
    df$term <- rownames(df)
    df$species <- sp
    df
  })
  coeff_stats <- do.call(rbind, coeff_stats_list)
  coeff_stats <- coeff_stats %>%
    mutate(term = case_when(
      term == "(Intercept)" ~ "Intercept",
      term == "poly(bin_median, 3)1" ~ "Linear",
      term == "poly(bin_median, 3)2" ~ "Quadratic",
      term == "poly(bin_median, 3)3" ~ "Cubic",
      TRUE ~ term
    )) %>%
    dplyr::select(species, term, p_value = `Pr(>|t|)`)
  
  # Join the p-values with the coefficient data.
  coeff_data <- left_join(coeff_data, coeff_stats, by = c("species", "term"))
  
  # Create labels: if p < 0.05 then "*" else the p_value rounded to 2 decimals.
  coeff_data <- coeff_data %>%
    mutate(label_text = ifelse(p_value < 0.05, "*", sprintf("%.2f", p_value)))
  
  coeff_data$species <- factor(coeff_data$species, levels = species_levels)
  
  # Create the grouped bar plot with p-value labels centered in each bar.
  plot_coeff <- ggplot(coeff_data, aes(x = term, y = value, fill = species)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = label_text, y = value/2), 
              color = "black", size = 3.5, 
              position = position_dodge(width = 0.8)) +
    scale_fill_manual(values = cb_palette) +
    labs(title = "Model Coefficients", 
         subtitle = expression(hat(Y) == beta[0] + beta[1]*x + beta[2]*x^2 + beta[3]*x^3),
         x = "Coefficient Term", 
         y = "Coefficient Value") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.position = "top",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  print(plot_coeff)
  
  # Save the coefficient plot.
  ggsave(coef_output, plot = plot_coeff, width = 10, height = 8, dpi = 300)
  
  # Return both plots as a list.
  return(list(combined_plot = final_plot, coeff_plot = plot_coeff))
}

plot_Proportions_TDiff_exp_slope_coeff <- function(data, coef_fig, output_figure) {
  # Process data using the custom NDVI_TDiffbin function and remove missing values
  data <- NDVI_TDiffbin(data, 5)
  data <- na.omit(data)
  
  # Load required libraries (assumed to be installed)
  library(ggplot2)
  library(nlme)
  library(dplyr)
  library(tibble)
  library(patchwork)
  library(purrr)
  library(car)      # for deltaMethod
  library(broom)
  library(tidyr)
  
  # Identify the value column and set species order and color palette
  value_col <- "avg_value"
  data$species <- factor(data$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  cb_palette <- c("Oak"   = "#E69F00",  # Orange
                  "Beech" = "#0072B2",  # Deep blue
                  "Spruce"= "#009E73",  # Bluish-green
                  "Pine"  = "#F0E442")  # Yellow
  
  # For Transpiration Deficit, use x = bin_median directly (assumed positive)
  data <- data %>% mutate(x = bin_median)
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  
  # Set threshold using the median of avg_value
  # threshold <- median(data_clean$avg_value, na.rm = TRUE)
  threshold <- 9
  
  # NONLINEAR MODELING PER SPECIES
  start_list <- list(a = 5, b = 7, c = 0.04)
  control_params <- nls.control(maxiter = 1200, minFactor = 1e-09)
  
  nls_models <- nlsList(avg_value ~ a + b * exp(-c * x) | species,
                        data = data_clean,
                        start = start_list,
                        control = control_params)
  print(summary(nls_models))
  
  # Extract coefficient summaries for a, b, and c
  model_summary <- summary(nls_models)
  coeff_a <- model_summary$coefficients[, , "a"]
  coeff_b <- model_summary$coefficients[, , "b"]
  coeff_c <- model_summary$coefficients[, , "c"]
  
  # Create a data frame with species names, estimates, and p-values
  results_df <- data.frame(
    species    = rownames(coeff_a),
    a_estimate = coeff_a[, "Estimate"],
    a_pvalue   = coeff_a[, "Pr(>|t|)"],
    b_estimate = coeff_b[, "Estimate"],
    b_pvalue   = coeff_b[, "Pr(>|t|)"],
    c_estimate = coeff_c[, "Estimate"],
    c_pvalue   = coeff_c[, "Pr(>|t|)"]
  )
  results_df$species <- factor(results_df$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  
  # Reshape the results into long format: first for estimates, then merge p-values
  coeffs_long <- results_df %>%
    pivot_longer(
      cols = c(a_estimate, b_estimate, c_estimate),
      names_to = "Coefficient",
      values_to = "Value",
      names_pattern = "([abc])_.*"
    ) %>%
    left_join(
      results_df %>%
        pivot_longer(
          cols = c(a_pvalue, b_pvalue, c_pvalue),
          names_to = "Coefficient",
          values_to = "pvalue",
          names_pattern = "([abc])_.*"
        ),
      by = c("species", "Coefficient")
    ) %>%
    mutate(
      label = if_else(pvalue < 0.05, "*", sprintf("%.3f", pvalue))
    )
  
  # Create the faceted coefficient bar plot with p-value labels
  p_coeff <- ggplot(coeffs_long, aes(x = species, y = Value, fill = species)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(aes(label = label, y = Value/2),
              position = position_dodge(width = 0.9),
              vjust = 0.5, color = "black", size = 4) +
    scale_fill_manual(values = cb_palette) +
    facet_wrap(~ Coefficient, scales = "free_y") +
    labs(title = "Coefficients by Species for NDVI Quantiles",
         x = "Species",
         y = "Coefficient Value") +
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
  print(p_coeff)
  
  # Save the coefficient plot
  ggsave(filename = coef_fig, plot = p_coeff, device = "png", width = 8, height = 6, dpi = 300)
  
  # -------------------------------
  # Additional analysis and combined panels
  # -------------------------------
  
  # Extract coefficients for a, b, and c per species
  coef_df <- as.data.frame(coef(nls_models), optional = TRUE) %>% 
    rownames_to_column(var = "species") %>%
    filter(!is.na(a))
  coef_df$species <- factor(coef_df$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  
  # Calculate x50: x50 = -log((threshold - a)/b) / c
  coef_df <- coef_df %>%
    mutate(x50 = ifelse((threshold - a) > 0 & b > 0,
                        -log((threshold - a)/b) / c,
                        NA))
  
  # Calculate slope at x50: slope50 = -c * (threshold - a)
  coef_df <- coef_df %>% mutate(slope50 = -c * (threshold - a))
  
  # Generate predictions for each species (Panel A)
  pred_list <- data_clean %>%
    group_by(species) %>%
    do({
      sp <- unique(.$species)
      x_seq <- seq(min(.$x, na.rm = TRUE), max(.$x, na.rm = TRUE), length.out = 100)
      sp_model <- nls_models[[as.character(sp)]]
      pred <- predict(sp_model, newdata = data.frame(x = x_seq))
      data.frame(x = x_seq, pred = pred)
    })
  pred_all <- bind_rows(pred_list)
  
  p_combined <- ggplot() +
    geom_point(data = data_clean, aes(x = x, y = avg_value, color = species)) +
    geom_line(data = pred_all, aes(x = x, y = pred, color = species), size = 1) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "black", size = 1) +
    annotate("text", x = min(data_clean$x, na.rm = TRUE), y = threshold,
             # label = "median", hjust = -0.1, vjust = -0.3, fontface = "italic", size = 4) +
             label = "", hjust = -0.1, vjust = -0.3, fontface = "italic", size = 4) +
    scale_color_manual(values = cb_palette) +
    labs(x = "Transpiration Deficit", y = "NDVI Quantiles") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    ) +
    labs(caption = "(a)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0),
          plot.caption.position = "plot") +
    coord_cartesian(clip = "off")
  
  # Panel B: Bar plot of x50 values
  p_x50 <- ggplot(coef_df, aes(x = species, y = x50, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    labs(x = "Species", y = "Transpiration Deficit") +
    scale_fill_manual(values = cb_palette) +
    expand_limits(y = 0) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(caption = "(b)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0),
          plot.caption.position = "plot")
  
  # Panel C: Bar plot of absolute slope at x50 with error bars and R² annotations
  stats_list <- lapply(levels(data_clean$species), function(sp) {
    mod <- nls_models[[sp]]
    coefs <- coef(mod)
    # Compute slope50 = -c * (threshold - a)
    slope50_est <- -coefs["c"] * (threshold - coefs["a"])
    # Compute standard error using deltaMethod for f(a, c) = c*(a - threshold)
    dm_result <- deltaMethod(mod, paste0("c*(a - ", threshold, ")"),
                             parameterNames = c("a", "b", "c"))
    se <- dm_result$SE
    # Compute degrees of freedom from the model summary
    df_mod <- summary(mod)$df[2]
    t_val <- slope50_est / se
    p_val <- 2 * (1 - pt(abs(t_val), df_mod))
    
    # Compute R² for the species model
    df_sp <- data_clean %>% filter(species == sp)
    fitted_vals <- predict(mod, newdata = df_sp)
    r_squared <- 1 - sum((df_sp[[value_col]] - fitted_vals)^2) /
      sum((df_sp[[value_col]] - mean(df_sp[[value_col]]))^2)
    
    tibble(species = sp,
           slope50 = slope50_est,
           slope_abs = abs(slope50_est),
           se = se,
           p_val = p_val,
           r_squared = r_squared)
  })
  stats_df <- bind_rows(stats_list)
  stats_df$species <- factor(stats_df$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  
  stats_df <- stats_df %>%
    mutate(label_text = sprintf("%.2f%s", r_squared, ifelse(p_val < 0.05, "*", "")))
  
  p_slope <- ggplot(stats_df, aes(x = species, y = slope_abs, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, slope_abs - se), ymax = slope_abs + se), width = 0.2) +
    geom_text(aes(y = slope_abs/2, label = label_text), color = "black", size = 4) +
    labs(x = "Species", y = "Absolute Slope") +
    scale_fill_manual(values = cb_palette) +
    expand_limits(y = 0) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(caption = "(c)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0),
          plot.caption.position = "plot")
  
  # Combine panels: Panel (a) on the left; Panel (b) (top) and Panel (c) (bottom) in right column
  final_plot <- p_combined + (p_x50 / p_slope) + plot_layout(widths = c(2, 1))
  print(final_plot)
  
  # Save the combined final plot
  ggsave(filename = output_figure, plot = final_plot, device = "png", width = 10, height = 8, dpi = 300)
}

plot_Proportions_TDiff_linear_parabola_slope_coeff <- function(data, coef_fig, output_figure) {
  
  # Load required libraries (if not already loaded)
  library(ggplot2)
  library(nlme)
  library(dplyr)
  library(tibble)
  library(patchwork)
  library(purrr)
  library(car)      # for deltaMethod (if needed)
  library(broom)
  library(tidyr)
  
  # -------------------------------
  # DATA PROCESSING
  # -------------------------------
  # Process data using your custom NDVI_TDiffbin function.
  # (Assumes 'data' is the raw data; adjust if needed.)
  data <- NDVI_TDiffbin(data, 5)
  data <- na.omit(data)
  
  # Identify the correct value column
  value_col <- "avg_value"
  
  # Order species and define color palette
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  data$species <- factor(data$species, levels = species_order)
  cb_palette <- c("Oak"   = "#E69F00",  # Orange
                  "Beech" = "#0072B2",  # Deep blue
                  "Spruce"= "#009E73",  # Bluish-green
                  "Pine"  = "#F0E442")  # Yellow
  
  # For Transpiration Deficit, use x = bin_median directly (positive values)
  data <- data %>% mutate(x = bin_median)
  
  # Clean data: remove rows with missing or non-finite values
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  
  # Set threshold using the median of avg_value
  threshold <- median(data_clean$avg_value, na.rm = TRUE)
  
  # -------------------------------
  # MODEL FITTING: DIFFERENT MODELS PER SPECIES
  # -------------------------------
  # For Oak and Beech: simple linear models (avg_value ~ x)
  data_lin <- data_clean %>% filter(species %in% c("Oak", "Beech"))
  lm_models_lin <- data_lin %>% 
    group_by(species) %>% 
    group_modify(~ {
      mod <- lm(avg_value ~ x, data = .x)
      tibble(model = list(mod))
    })
  
  tidy_lin <- lm_models_lin %>% 
    mutate(tidy = map(model, broom::tidy)) %>% 
    unnest(tidy) %>% 
    mutate(model_type = "linear",
           term = case_when(term == "(Intercept)" ~ "a",
                            term == "x" ~ "b",
                            TRUE ~ term))
  
  # For Spruce and Pine: quadratic models using nlsList
  data_poly <- data_clean %>% filter(species %in% c("Spruce", "Pine"))
  start_list_poly <- list(a = 0.1, b = 0.1, c = 0.1)
  nls_models_poly <- nlsList(
    avg_value ~ a + b * x + c * I(x^2) | species,
    data = data_poly,
    start = start_list_poly,
    control = nls.control()  # default control settings
  )
  
  tidy_poly <- map_df(nls_models_poly, broom::tidy, .id = "species") %>%
    mutate(model_type = "quadratic")
  
  # Combine the tidy coefficients for coefficient plotting
  tidy_all <- bind_rows(tidy_lin, tidy_poly)
  tidy_all$species <- factor(tidy_all$species, levels = species_order)
  
  # Create a label column for coefficient plot (for reference)
  tidy_all <- tidy_all %>%
    mutate(label = if_else(p.value < 0.05,
                           sprintf("*"),
                           sprintf("p=%.3f", p.value)
    ))
  
  # -------------------------------
  # PANEL: COEFFICIENT BAR PLOT (Linear + Quadratic Coefficients)
  # -------------------------------
  p_coeff <- ggplot(tidy_all, aes(x = species, y = estimate, fill = species)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(aes(label = label, y = estimate/2),
              position = position_dodge(width = 0.9),
              vjust = 0.5, color = "black", size = 4) +
    scale_fill_manual(values = cb_palette) +
    facet_wrap(~ term, scales = "free_y") +
    labs(title = "Model Coefficients by Species for NDVI Proportions",
         x = "Species",
         y = "Coefficient Value") +
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
  print(p_coeff)
  ggsave(filename = coef_fig,
         plot = p_coeff, device = "png", width = 8, height = 6, dpi = 300)
  
  # -------------------------------
  # PANEL A: Observed Data and Fitted Curves
  # -------------------------------
  # Predictions for linear models (Oak and Beech): full x-range
  pred_lin <- data_lin %>% 
    group_by(species) %>%
    do({
      sp <- unique(.$species)
      x_seq <- seq(min(.$x, na.rm = TRUE), max(.$x, na.rm = TRUE), length.out = 100)
      model <- lm_models_lin$model[[which(lm_models_lin$species == sp)]]
      pred <- predict(model, newdata = data.frame(x = x_seq))
      data.frame(species = sp, x = x_seq, pred = pred, model_type = "linear")
    })
  
  # Predictions for quadratic models (Spruce and Pine)
  pred_poly <- data_poly %>% 
    group_by(species) %>%
    do({
      sp <- unique(.$species)
      x_seq <- seq(min(.$x, na.rm = TRUE), max(.$x, na.rm = TRUE), length.out = 100)
      model <- nls_models_poly[[as.character(sp)]]
      pred <- predict(model, newdata = data.frame(x = x_seq))
      data.frame(species = sp, x = x_seq, pred = pred, model_type = "quadratic")
    })
  pred_poly <- pred_poly %>% mutate(pred = as.numeric(pred))
  
  # For quadratic models, compute the vertex and restrict predictions.
  coef_df_vertex <- tidy_poly %>% 
    pivot_wider(id_cols = species, names_from = term, values_from = estimate) %>%
    mutate(vertex = -b / (2 * c))
  print(coef_df_vertex)  # Check vertex values
  
  x_range_df_poly <- data_poly %>%
    group_by(species) %>%
    summarize(x_max = max(x, na.rm = TRUE), .groups = "drop")
  
  pred_poly <- pred_poly %>%
    left_join(coef_df_vertex %>% select(species, vertex), by = "species") %>%
    left_join(x_range_df_poly, by = "species") %>%
    mutate(eff_vertex = if_else(vertex > x_max, x_max, vertex)) %>%
    filter(x <= eff_vertex)
  
  # Combine predictions from both model types
  pred_all <- bind_rows(pred_lin, pred_poly)
  
  # Create Panel A: Plot observed data and fitted lines
  p_combined <- ggplot() +
    geom_point(data = data_clean, aes(x = x, y = avg_value, color = species)) +
    geom_line(data = pred_all, aes(x = x, y = pred, color = species, group = species), linewidth = 1) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "black", linewidth = 1) +
    annotate("text", x = min(data_clean$x, na.rm = TRUE), y = threshold,
             label = "median", hjust = -0.1, vjust = -0.3, fontface = "italic", size = 4) +
    scale_color_manual(values = cb_palette) +
    labs(x = "Transpiration Deficit", y = "NDVI Proportions") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    ) +
    labs(caption = "(a)") +
    theme(
      plot.caption = element_text(face = "bold", size = 16, hjust = 0),
      plot.caption.position = "plot"
    ) +
    coord_cartesian(clip = "off")
  
  print(p_combined)
  
  # -------------------------------
  # PANEL B: Calculate and Plot x50 Values (x where y equals the overall median)
  # -------------------------------
  # For linear models:
  lin_coef <- tidy_lin %>%
    filter(species %in% c("Oak", "Beech")) %>%
    select(species, term, estimate) %>%
    pivot_wider(names_from = term, values_from = estimate) %>%
    mutate(x50 = (threshold - a) / b)
  
  lin_x50 <- lin_coef %>%
    select(species, x50) %>%
    left_join(lm_models_lin %>% 
                mutate(p_val = map_dbl(model, ~ summary(.x)$coefficients["x", "Pr(>|t|)"])) %>% 
                select(species, p_val), 
              by = "species") %>%
    mutate(model_type = "linear")
  
  # For quadratic models:
  calc_x50 <- function(a, b, c, threshold, x_range) {
    discriminant <- b^2 - 4 * c * (a - threshold)
    if (discriminant < 0) return(NA)
    sol1 <- (-b + sqrt(discriminant)) / (2 * c)
    sol2 <- (-b - sqrt(discriminant)) / (2 * c)
    valid_solutions <- c(sol1, sol2)[c(sol1, sol2) >= x_range[1] & c(sol1, sol2) <= x_range[2]]
    if (length(valid_solutions) == 0) {
      all_sols <- c(sol1, sol2)
      diff <- abs(all_sols - mean(x_range))
      return(all_sols[which.min(diff)])
    } else if (length(valid_solutions) == 1) {
      return(valid_solutions)
    } else {
      return(valid_solutions[which.min(abs(valid_solutions - median(x_range)))])
    }
  }
  
  coef_df_quad <- tidy_all %>% 
    filter(model_type == "quadratic") %>%
    select(species, term, estimate) %>%
    pivot_wider(names_from = term, values_from = estimate)
  
  x_range_df_all <- data_clean %>% 
    group_by(species) %>% 
    summarize(x_min = min(x, na.rm = TRUE),
              x_max = max(x, na.rm = TRUE),
              .groups = "drop")
  
  coef_df_quad <- coef_df_quad %>% 
    left_join(x_range_df_all, by = "species") %>%
    mutate(x50 = pmap_dbl(list(a, b, c, rep(threshold, n()), list(c(x_min, x_max))),
                          ~ calc_x50(..1, ..2, ..3, ..4, ..5)))
  
  quad_pvals <- tibble(species = c("Spruce", "Pine")) %>%
    rowwise() %>%
    mutate(p_val = {
      model <- nls_models_poly[[as.character(species)]]
      sum_model <- summary(model)
      sum_model$coefficients["b", "Pr(>|t|)"]
    }) %>%
    ungroup()
  
  quad_x50 <- coef_df_quad %>%
    filter(species %in% c("Spruce", "Pine")) %>%
    select(species, x50) %>%
    left_join(quad_pvals, by = "species") %>%
    mutate(model_type = "quadratic")
  
  x50_all <- bind_rows(lin_x50, quad_x50) %>%
    arrange(species)
  x50_all$species <- factor(x50_all$species, levels = species_order)
  
  x50_all <- x50_all %>%
    mutate(label = if_else(p_val < 0.05, "*", sprintf("p=%.3f", p_val)))
  
  print(x50_all)
  
  p_x50_new <- ggplot(x50_all, aes(x = species, y = x50, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = label, y = x50/2), color = "black", size = 4) +
    labs(x = "Species", y = "Transpiration Deficit") +
    scale_fill_manual(values = cb_palette) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(caption = "(b)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0),
          plot.caption.position = "plot")
  
  print(p_x50_new)
  
  # -------------------------------
  # PANEL C: Calculate and Plot Slope at Median (|deriv|) with Error Bars and R² Annotation
  # -------------------------------
  # For linear models, the slope is simply the coefficient "x"
  stats_lin_slope <- lm_models_lin %>% 
    mutate(mod_sum = map(model, summary)) %>%
    mutate(
      slope = map_dbl(mod_sum, ~ abs(.x$coefficients["x", "Estimate"])),
      se = map_dbl(mod_sum, ~ .x$coefficients["x", "Std. Error"]),
      p_val = map_dbl(mod_sum, ~ .x$coefficients["x", "Pr(>|t|)"]),
      r_square = map_dbl(mod_sum, ~ .x$r.squared)
    ) %>% 
    select(species, slope, se, p_val, r_square) %>%
    mutate(model_type = "linear")
  
  # Helper function for quadratic models
  compute_quad_slope <- function(sp) {
    x50_local <- coef_df_quad$x50[coef_df_quad$species == sp]
    model_obj <- nls_models_poly[[sp]]
    b_hat <- coef(model_obj)["b"]
    c_hat <- coef(model_obj)["c"]
    deriv <- b_hat + 2 * c_hat * x50_local
    slope <- abs(deriv)
    s_mod <- summary(model_obj)
    cov_mat <- as.matrix(s_mod$cov.unscaled * s_mod$sigma^2)
    se <- sqrt( cov_mat["b", "b"] +
                  (2 * x50_local)^2 * cov_mat["c", "c"] +
                  2 * (2 * x50_local) * cov_mat["b", "c"] )
    p_val <- 2 * pnorm(-abs(deriv / se))
    dat_sp <- data_poly %>% filter(species == sp)
    preds <- predict(model_obj, newdata = dat_sp)
    r_sq <- 1 - sum((dat_sp$avg_value - preds)^2) / sum((dat_sp$avg_value - mean(dat_sp$avg_value))^2)
    tibble(species = sp, slope = slope, se = se, p_val = p_val, r_square = r_sq, model_type = "quadratic")
  }
  
  species_char <- as.character(c("Spruce", "Pine"))
  result_list <- lapply(species_char, compute_quad_slope)
  stats_quad_slope <- bind_rows(result_list)
  
  stats_slope_all <- bind_rows(stats_lin_slope, stats_quad_slope) %>%
    arrange(species)
  stats_slope_all$species <- factor(stats_slope_all$species, levels = species_order)
  
  stats_slope_all <- stats_slope_all %>%
    mutate(label = sprintf("R²=%.2f%s", r_square, if_else(p_val < 0.05, "*", "")))
  
  print(stats_slope_all)
  
  p_slope <- ggplot(stats_slope_all, aes(x = species, y = slope, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, slope - se), ymax = slope + se), width = 0.2) +
    geom_text(aes(y = slope/2, label = label), color = "black", size = 4) +
    labs(x = "Species", y = "Absolute Slope") +
    scale_fill_manual(values = cb_palette) +
    theme_minimal() +
    labs(caption = "(c)") +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.caption = element_text(face = "bold", size = 16, hjust = 0),
      plot.caption.position = "plot"
    )
  print(p_slope)
  
  # -------------------------------
  # COMBINE PANELS (A, B, and C) AND SAVE FINAL PLOT
  # -------------------------------
  final_plot <- p_combined + (p_x50_new / p_slope) + plot_layout(widths = c(2, 1))
  print(final_plot)
  
  ggsave(filename = output_figure,
         plot = final_plot, device = "png", width = 10, height = 8, dpi = 300)
}

plot_Proportions_TDiff_poly_2_slope_coeff <- function(data, coef_fig, output_figure) {
  
  # -------------------------------
  # SETUP: Load required libraries
  # -------------------------------
  require(ggplot2)
  require(nlme)
  require(dplyr)
  require(tibble)
  require(patchwork)
  require(purrr)
  require(car)      # for deltaMethod
  require(broom)
  require(tidyr)
  
  # -------------------------------
  # DATA PREPARATION
  # -------------------------------
  # Process data using your custom binning function for proportions.
  # (Make sure that NDVI_TDiffbin() returns a column "avg_value" and "bin_median".)
  data <- NDVI_TDiffbin(data)
  data <- na.omit(data)
  
  # Identify the proportion column and set species order and color palette.
  value_col <- "avg_value"
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  data$species <- factor(data$species, levels = species_order)
  cb_palette <- c("Oak"   = "#E69F00",  # Orange
                  "Beech" = "#0072B2",  # Deep blue
                  "Spruce"= "#009E73",  # Bluish-green
                  "Pine"  = "#F0E442")  # Yellow
  
  # For Transpiration Deficit, use x = bin_median (assumed positive)
  data <- data %>% mutate(x = bin_median)
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  data_clean$species <- factor(data_clean$species, levels = species_order)
  
  # Set the threshold (NDVI value) at which we calculate x50 and the slope.
  threshold <- 1
  
  # -------------------------------
  # FITTING THE QUADRATIC MODEL PER SPECIES (nlsList)
  # -------------------------------
  start_list_poly <- list(a = 0.5, b = 0.1, c = 0.1)
  nls_models_poly <- nlsList(
    avg_value ~ a + b * x + c * I(x^2) | species,
    data = data_clean,
    start = start_list_poly,
    control = nls.control()
  )
  print(summary(nls_models_poly))
  
  # -------------------------------
  # COEFFICIENT BAR PLOT (Panel: Coef)
  # -------------------------------
  tidy_poly <- map_df(nls_models_poly, broom::tidy, .id = "species")
  tidy_poly$species <- factor(tidy_poly$species, levels = species_order)
  tidy_poly <- tidy_poly %>%
    mutate(label = if_else(p.value < 0.05, "*", sprintf("p=%.2f", p.value)))
  
  p_coeff <- ggplot(tidy_poly, aes(x = species, y = estimate, fill = species)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(aes(label = label, y = estimate/2),
              position = position_dodge(width = 0.9),
              vjust = 0.5, color = "black", size = 4) +
    scale_fill_manual(values = cb_palette) +
    facet_wrap(~ term, scales = "free_y") +
    labs(title = "Quadratic Coefficients by Species for NDVI Proportions",
         subtitle = expression(NDVI == a + b * x + c * x^2),
         x = "Species",
         y = "Coefficient Value") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.background = element_rect(fill = "white", color = "white"),
          panel.background = element_rect(fill = "white"),
          legend.background = element_rect(fill = "white", color = "white"),
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
          plot.subtitle = element_text(hjust = 0.5, size = 14, face = "italic", color = "black"),
          axis.title = element_text(face = "bold"),
          axis.text = element_text(color = "black"),
          panel.border = element_rect(color = NA, fill = NA, linewidth = 0),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "top",
          strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
          strip.text = element_text(face = "bold", size = 12))
  print(p_coeff)
  ggsave(filename = coef_fig,
         plot = p_coeff, device = "png", width = 8, height = 6, dpi = 300)
  
  # -------------------------------
  # PANEL A: Observed Data and Fitted Curves
  # -------------------------------
  pred_list <- data_clean %>%
    group_by(species) %>%
    do({
      sp <- unique(.$species)
      x_seq <- seq(min(.$x, na.rm = TRUE), max(.$x, na.rm = TRUE), length.out = 100)
      mod <- nls_models_poly[[as.character(sp)]]
      pred <- predict(mod, newdata = data.frame(x = x_seq))
      data.frame(species = sp, x = x_seq, pred = pred)
    })
  pred_all <- bind_rows(pred_list) %>% mutate(pred = as.numeric(pred))
  
  # Calculate vertex for each quadratic curve: vertex = -b / (2*c)
  coef_df_vertex <- tidy_poly %>%
    select(species, term, estimate) %>%
    pivot_wider(names_from = term, values_from = estimate) %>%
    mutate(vertex = -b / (2 * c))
  # print(coef_df_vertex)
  
  # Get observed x-range for each species
  x_range_df <- data_clean %>% group_by(species) %>%
    summarize(x_min = min(x), x_max = max(x), .groups = "drop")
  
  # Join vertex and x-range with predictions and flag if vertex is in range
  pred_all <- pred_all %>%
    left_join(coef_df_vertex %>% select(species, vertex), by = "species") %>%
    left_join(x_range_df, by = "species") %>%
    mutate(vertex_in_range = (vertex >= x_min & vertex <= x_max))
  
  # For species where vertex is in range AND vertex is greater than the minimum observed x,
  # we only plot predictions up to vertex; otherwise, we plot all predictions.
  pred_all_final <- pred_all %>%
    group_by(species) %>%
    filter(if_else(vertex_in_range & (vertex > min(x)),
                   x <= vertex,
                   TRUE)) %>%
    ungroup()
  
  p_combined <- ggplot() +
    geom_point(data = data_clean, aes(x = x, y = avg_value, color = species)) +
    geom_line(data = pred_all_final, aes(x = x, y = pred, color = species, group = species), size = 1) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "black", size = 1) +
    annotate("text", x = min(data_clean$x, na.rm = TRUE), y = threshold,
             label = "median", hjust = -2.1, vjust = -0.3, fontface = "italic", size = 4) +
    scale_color_manual(values = cb_palette) +
    labs(x = "Transpiration Deficit", y = "NDVI") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.background = element_rect(fill = "white", color = "white"),
          panel.background = element_rect(fill = "white"),
          legend.background = element_rect(fill = "white", color = "white"),
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
          axis.title = element_text(face = "bold"),
          axis.text = element_text(color = "black"),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "top") +
    labs(caption = "(a)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0)) +
    coord_cartesian(clip = "off")
  
  # -------------------------------
  # PANEL B: Bar Plot of x50 Values
  # -------------------------------
  calc_x50 <- function(a, b, c, threshold, x_range) {
    discriminant <- b^2 - 4 * c * (a - threshold)
    if (discriminant < 0) return(NA)
    sol1 <- (-b + sqrt(discriminant)) / (2 * c)
    sol2 <- (-b - sqrt(discriminant)) / (2 * c)
    valid_solutions <- c(sol1, sol2)[c(sol1, sol2) >= x_range[1] & c(sol1, sol2) <= x_range[2]]
    if (length(valid_solutions) == 0) {
      all_sols <- c(sol1, sol2)
      diff <- abs(all_sols - mean(x_range))
      return(all_sols[which.min(diff)])
    } else if (length(valid_solutions) == 1) {
      return(valid_solutions)
    } else {
      return(valid_solutions[which.min(abs(valid_solutions - median(x_range)))])
    }
  }
  
  stats_list <- list()
  for (sp in levels(data_clean$species)) {
    mod <- nls_models_poly[[sp]]
    coefs <- coef(mod)  # a, b, c
    sp_data <- data_clean %>% filter(species == sp)
    x_range <- c(min(sp_data$x, na.rm = TRUE), max(sp_data$x, na.rm = TRUE))
    x50_val <- calc_x50(coefs["a"], coefs["b"], coefs["c"], threshold, x_range)
    
    # Slope at x50: derivative = b + 2*c*x50
    slope_val <- coefs["b"] + 2 * coefs["c"] * x50_val
    expr <- paste0("b + 2*c*(", x50_val, ")")
    vcov_mat <- vcov(mod)
    dm <- deltaMethod(coefs, expr, vcov_mat)
    slope_se <- as.numeric(dm["SE"])
    
    # Compute p-value using a z statistic approximation
    z_val <- slope_val / slope_se
    slope_p <- 2 * (1 - pnorm(abs(z_val)))
    
    stats_list[[sp]] <- data.frame(species = sp, x50 = x50_val, slope_abs = abs(slope_val), se = slope_se, slope_p = slope_p)
  }
  
  stats_df <- bind_rows(stats_list)
  stats_df$species <- factor(stats_df$species, levels = species_order)
  print(stats_df %>% select(species, x50, slope_abs, se, slope_p))
  
  p_x50 <- ggplot(stats_df, aes(x = species, y = x50, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    labs(x = "Species", y = "Transpiration Deficit (x50)") +
    scale_fill_manual(values = cb_palette) +
    scale_x_discrete(limits = species_order) +
    expand_limits(y = 0) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          axis.title.x = element_blank(),
          axis.title = element_text(face = "bold"),
          axis.text = element_text(color = "black"),
          plot.background = element_rect(fill = "white", color = "white"),
          panel.background = element_rect(fill = "white"),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(caption = "(b)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0))
  # print(p_x50)
  
  # -------------------------------
  # PANEL C: Bar Plot of Absolute Slope at x50 with Error Bars and Conditional Labels
  # -------------------------------
  rsq_df <- map_df(levels(data_clean$species), function(sp) {
    sp_data <- dplyr::filter(data_clean, species == sp)
    mod <- nls_models_poly[[sp]]
    pred <- predict(mod, newdata = sp_data)
    r2_val <- 1 - sum((sp_data[[value_col]] - pred)^2) / sum((sp_data[[value_col]] - mean(sp_data[[value_col]]))^2)
    data.frame(species = sp, r2 = r2_val)
  })
  
  stats_df <- left_join(stats_df, rsq_df, by = "species")
  
  stats_df <- stats_df %>%
    mutate(slope_label = if_else(slope_p < 0.05,
                                 sprintf("%.2f*", r2),
                                 sprintf("p=%.2f\nR² = %.2f", slope_p, r2)))
  
  p_slope <- ggplot(stats_df, aes(x = species, y = slope_abs, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, slope_abs - se), ymax = slope_abs + se), width = 0.2, color = "black") +
    geom_text(aes(y = slope_abs/2, label = slope_label), color = "black", size = 4) +
    labs(x = "Species", y = "Absolute Slope at x50") +
    scale_fill_manual(values = cb_palette) +
    scale_x_discrete(limits = species_order) +
    theme_minimal() +
    labs(caption = "(c)") +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          axis.title = element_text(face = "bold"),
          axis.text = element_text(color = "black"),
          plot.background = element_rect(fill = "white", color = "white"),
          panel.background = element_rect(fill = "white"),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.caption = element_text(face = "bold", size = 16, hjust = 0),
          plot.caption.position = "plot")
  # print(p_slope)
  
  # -------------------------------
  # COMBINE PANELS AND SAVE FINAL PLOT
  # -------------------------------
  final_plot <- p_combined + (p_x50 / p_slope) + plot_layout(widths = c(2, 1))
  print(final_plot)
  
  ggsave(filename = output_figure,
         plot = final_plot, device = "png", width = 10, height = 8, dpi = 300)
}

plot_NDVI_PDM_TDiff_slope <- function(species_data, output_file) {
  # Load required libraries
  library(lme4)      # For mixed-effects modeling
  library(dplyr)     # For data manipulation
  library(ggplot2)   # For plotting
  library(patchwork) # For combining ggplots
  
  # Identify the correct value column (Quantiles or Proportions)
  value_col <- if ("Quantiles" %in% names(species_data)) "Quantiles" else "Proportions"
  
  # Create the binned dataset and remove NA rows.
  NDVI_TDiffbin_df <- NDVI_TDiffbin(species_data)
  NDVI_TDiffbin_df <- na.omit(NDVI_TDiffbin_df)
  
  # Define custom color palette and species order.
  cb_palette <- c("Oak"   = "#E69F00",  # Orange
                  "Beech" = "#0072B2",  # A deep blue with lower luminance than "#56B4E9"
                  "Spruce"= "#009E73",  # Bluish-green
                  "Pine"  = "#F0E442") 
  
  species_levels <- c("Oak", "Beech", "Spruce", "Pine")
  NDVI_TDiffbin_df$species <- factor(NDVI_TDiffbin_df$species, levels = species_levels)
  
  ## Panel A: Mixed-Effects Model Plot for avg_value (Quantiles or Proportions) vs bin_median
  model <- lmer(avg_value ~ poly(bin_median, 3) + (poly(bin_median, 3) | species),
                data = NDVI_TDiffbin_df)
  
  # Create prediction data for each species.
  pred_data <- NDVI_TDiffbin_df %>%
    group_by(species) %>%
    summarise(min_bin = min(bin_median),
              max_bin = max(bin_median),
              .groups = "drop") %>%
    group_by(species) %>%
    do(data.frame(bin_median = seq(.$min_bin, .$max_bin, length.out = 100))) %>%
    ungroup()
  
  # Ensure the species factor is in the desired order.
  pred_data$species <- factor(pred_data$species, levels = species_levels)
  
  # Generate predicted values (including random effects).
  pred_data$predicted <- predict(model, newdata = pred_data, re.form = NULL)
  
  # Calculate reference line value.
  line_val <- 0.98
  
  # Create Panel A plot.
  plot_mixed <- ggplot(NDVI_TDiffbin_df, aes(x = bin_median, y = avg_value, color = species)) +
    geom_point() +
    geom_line(data = pred_data, aes(x = bin_median, y = predicted, color = species), size = 1) +
    geom_hline(yintercept = line_val, linetype = "dashed", color = "black", size = 1) +
    scale_color_manual(values = cb_palette) +
    labs(x = "Transpiration Deficit (bin_median)",
         y = paste("Average", value_col)) +
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
  
  ## Panel B: Bar Plot of Local Slopes at 50% of Maximum Average Value
  # Compute the local slope at the point where the predicted value reaches 50% of its species-specific maximum.
  local_slope_data <- pred_data %>%
    group_by(species) %>%
    do({
      df <- .
      max_pred <- max(df$predicted, na.rm = TRUE)
      thresh_val <- 0.5 * max_pred
      idx <- which.min(abs(df$predicted - thresh_val))
      win_start <- max(1, idx - 5)
      win_end <- min(nrow(df), idx + 5)
      df_window <- df[win_start:win_end, ]
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
                               paste0("p: ", signif(p_value, 3), "\nR²: ", round(r2, 3))))
  
  # Create the bar plot with text labels for p value and R².
  p_bar <- ggplot(local_slope_data, aes(x = reorder(species, slope_abs), y = slope_abs, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = slope_abs - slope_se, ymax = slope_abs + slope_se), width = 0.2, color = "black") +
    geom_text(aes(label = label_text), position = position_stack(vjust = 0.5), size = 4, color = "black") +
    scale_fill_manual(values = cb_palette) +
    labs(x = "Species",
         y = paste("Absolute Slope at 0.98 PDM")) +
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
  
  # Combine the two panels side by side.
  final_plot <- plot_mixed + p_bar + plot_layout(widths = c(2, 1))
  
  # Print and save the final plot.
  print(final_plot)
  ggsave(output_file, plot = final_plot, width = 12, height = 6, dpi = 300)
  
  return(final_plot)
}

plot_NDVI_TDiffbin_line <- function(data, figure_output = NULL) {
  # Process the data using your NDVI_PSIbin function
  df <- NDVI_TDiffbin(data)
  
  # Ensure species are ordered as specified
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  df$species <- factor(df$species, levels = species_order)
  
  # Define the color palette for the species
  cb_palette <- c("Oak"   = "#E69F00",   # Orange
                  "Beech" = "#0072B2",   # Deep blue
                  "Spruce"= "#009E73",   # Bluish-green
                  "Pine"  = "#F0E442")   # Yellow
  
  # Create the ggplot: points and lines
  p <- ggplot(df, aes(x = bin_median, y = avg_value, color = species, group = species)) +
    geom_point() +
    geom_line() +
    scale_color_manual(values = cb_palette) +
    labs(title = "NDVI vs TDiff Bin",
         x = "Bin Median",
         y = "Avg Value") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    )
  
  # Save the figure if a file name is provided
  if (!is.null(figure_output)) {
    ggsave(filename = figure_output, plot = p)
  }
  
  # Return the plot object
  return(p)
}

plot_NDVI_TDiffbin_linear <- function(data, save_slope_fig) {
  
  # Process data with your custom NDVI_PSIbin function and remove missing values
  data <- na.omit(data)
  data <- NDVI_TDiffbin(data)
  
  # Load required libraries
  library(ggplot2)
  library(nlme)      # for lmList
  library(dplyr)
  library(tibble)
  library(patchwork)
  library(purrr)
  library(tidyr)
  
  # Identify the value column and order species; define the color palette
  value_col <- "avg_value"
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  data$species <- factor(data$species, levels = species_order)
  cb_palette <- c("Oak"   = "#E69F00",   # Orange
                  "Beech" = "#0072B2",   # Deep blue
                  "Spruce"= "#009E73",   # Bluish-green
                  "Pine"  = "#F0E442")   # Yellow
  
  # Create a positive soil water potential (x = -bin_median)
  data <- data %>% mutate(x = -bin_median)
  
  # Clean data: remove rows with missing or non-finite values
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  
  # Compute species-specific medians of avg_value
  medians_df <- data_clean %>% 
    group_by(species) %>% 
    summarize(median_y = median(avg_value))
  
  #### LINEAR MODELING PER SPECIES ####
  lm_models <- lmList(avg_value ~ x | species, data = data_clean)
  print(summary(lm_models))
  
  # Extract coefficients for each species into a data frame
  # The coefficients are stored as `(Intercept)` and `x`
  coef_df <- as.data.frame(coef(lm_models), optional = TRUE) %>% 
    rownames_to_column(var = "species") %>%
    filter(!is.na(`(Intercept)`))
  coef_df$species <- factor(coef_df$species, levels = species_order)
  
  # Merge with species-specific medians and compute x50: where fitted line equals the species median
  # For a linear model: (Intercept) + x * x50 = median_y  ->  x50 = (median_y - (Intercept)) / x
  coef_df <- coef_df %>% 
    left_join(medians_df, by = "species") %>%
    mutate(x50 = (median_y - `(Intercept)`) / x)
  
  #### PANEL A: Observed Data and Fitted Lines (Faceted by Species) ####
  # Generate predicted values for each species over a sequence of x values
  pred_list <- data_clean %>%
    group_by(species) %>%
    do({
      sp <- unique(.$species)
      x_seq <- seq(min(.$x, na.rm = TRUE), max(.$x, na.rm = TRUE), length.out = 100)
      sp_model <- lm_models[[as.character(sp)]]
      pred <- predict(sp_model, newdata = data.frame(x = x_seq))
      data.frame(x = x_seq, pred = pred, species = sp)
    })
  pred_all <- bind_rows(pred_list)
  
  # Reverse the x-axis so that the original (negative) soil water potential is shown
  x_scale <- scale_x_continuous(trans = "reverse", labels = function(x) -x)
  
  p_combined <- ggplot(data_clean, aes(x = x, y = avg_value, color = species)) +
    geom_point() +
    geom_line(data = pred_all, aes(x = x, y = pred), size = 1) +
    # Add a horizontal dashed line for each species median
    # geom_hline(data = medians_df, aes(yintercept = median_y),
    #            linetype = "dashed", color = "black", size = 1) +
    # Mark the point on the fitted line where the predicted value equals the species median
    geom_point(data = coef_df, aes(x = x50, y = median_y),
               #            shape = 17, size = 3, color = cb_palette) +
               shape = 17, size = 3, color = "black") +
    # facet_wrap(~ species) +
    scale_color_manual(values = cb_palette) +
    x_scale +
    labs(x = "Transpiration Deficit", y = "NDVI Quantiles") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    ) +
    labs(caption = "(a)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0),
          plot.caption.position = "plot") +
    coord_cartesian(clip = "off")
  
  #### PANEL B: Bar Plot of x50 Values ####
  p_x50 <- ggplot(coef_df, aes(x = species, y = -x50, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    labs(x = "Species", y = "Transpiration Deficit") +
    scale_fill_manual(values = cb_palette) +
    expand_limits(y = 0) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(caption = "(b)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0),
          plot.caption.position = "plot")
  
  #### PANEL C: Bar Plot of Absolute Slope with Error Bars and Annotations ####
  # For each species, extract the slope (coefficient for x) and its standard error and p-value from the lm model
  stats_list <- lapply(levels(data_clean$species), function(sp) {
    mod <- lm_models[[sp]]
    s <- summary(mod)$coefficients
    slope_est <- s["x", "Estimate"]
    se <- s["x", "Std. Error"]
    t_val <- slope_est / se
    df_mod <- summary(mod)$df[2]
    p_val <- 2 * (1 - pt(abs(t_val), df_mod))
    
    # Compute R² for the species model
    df_sp <- data_clean %>% filter(species == sp)
    fitted_vals <- predict(mod, newdata = df_sp)
    r_squared <- 1 - sum((df_sp[[value_col]] - fitted_vals)^2) /
      sum((df_sp[[value_col]] - mean(df_sp[[value_col]]))^2)
    
    tibble(species = sp,
           slope50 = slope_est,
           slope_abs = abs(slope_est),
           se = se,
           p_val = p_val,
           r_squared = r_squared)
  })
  stats_df <- bind_rows(stats_list)
  stats_df$species <- factor(stats_df$species, levels = species_order)
  
  # Create annotation labels for Panel C
  stats_df <- stats_df %>%
    mutate(label_text = ifelse(p_val < 0.05,
                               sprintf("%.2f*", r_squared),
                               sprintf("%.2f\np = %.3f", r_squared, p_val)))
  
  p_slope <- ggplot(stats_df, aes(x = species, y = slope_abs, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, slope_abs - se), ymax = slope_abs + se), width = 0.2) +
    geom_text(aes(y = slope_abs/2, label = label_text), color = "black", size = 4) +
    labs(x = "Species", y = "Absolute Slope") +
    scale_fill_manual(values = cb_palette) +
    expand_limits(y = 0) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(caption = "(c)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0),
          plot.caption.position = "plot")
  
  #### COMBINE PANELS A, B, and C INTO THE FINAL SLOPE FIGURE ####
  final_slope_plot <- p_combined + (p_x50 / p_slope) + plot_layout(widths = c(2, 1))
  print(final_slope_plot)
  
  # Save the slope figure to file
  dir.create(dirname(save_slope_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = save_slope_fig, plot = final_slope_plot, width = 10, height = 8, dpi = 300)
}

plot_Quantiles_TDiff_slope_each <- function(data, coef_fig, output_figure) {
  
  # Process data using the custom NDVI_TDiffbin function and remove missing values
  data <- NDVI_TDiffbin(data, 5)
  data <- na.omit(data)
  
  # Load required libraries
  library(ggplot2)
  library(nlme)
  library(dplyr)
  library(tibble)
  library(patchwork)
  library(purrr)
  library(car)      # (kept for future needs, e.g., deltaMethod)
  library(broom)
  library(tidyr)
  library(ggpattern)  # for patterned bars in Panel C
  
  # Identify the value column and set species order and color palette
  value_col <- "avg_value"
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  data$species <- factor(data$species, levels = species_order)
  cb_palette <- c("Oak"   = "#E69F00",  # Orange
                  "Beech" = "#0072B2",  # Deep blue
                  "Spruce"= "#009E73",  # Bluish-green
                  "Pine"  = "#F0E442")  # Yellow
  
  # For Transpiration Deficit, use x = bin_median directly
  data <- data %>% mutate(x = bin_median)
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  
  #### Compute species-specific medians ####
  medians_df <- data_clean %>% 
    group_by(species) %>% 
    summarise(median_value = median(.data[[value_col]]))
  
  #### Compute maximum and minimum transpiration deficit for each species ####
  max_x_df <- data_clean %>%
    group_by(species) %>%
    summarise(x_max = max(x))
  min_x_df <- data_clean %>%
    group_by(species) %>%
    summarise(x_min = min(x))
  
  #### NONLINEAR MODELING PER SPECIES ####
  start_list <- list(a = 5, b = 7, c = 0.04)
  control_params <- nls.control(maxiter = 1200, minFactor = 1e-09)
  
  nls_models <- nlsList(avg_value ~ a + b * exp(-c * x) | species,
                        data = data_clean,
                        start = start_list,
                        control = control_params)
  print(summary(nls_models))
  
  # Extract coefficients and join medians, x_max, and x_min
  coef_df <- as.data.frame(coef(nls_models), optional = TRUE) %>% 
    rownames_to_column(var = "species") %>%
    filter(!is.na(a))
  coef_df$species <- factor(coef_df$species, levels = species_order)
  coef_df <- left_join(coef_df, medians_df, by = "species")
  coef_df <- left_join(coef_df, max_x_df, by = "species")
  coef_df <- left_join(coef_df, min_x_df, by = "species")
  
  # Calculate x50: the x value corresponding to the species-specific median NDVI
  coef_df <- coef_df %>%
    mutate(x50 = ifelse((median_value - a) > 0 & b > 0,
                        -log((median_value - a)/b) / c,
                        NA))
  
  #### PANEL A: Observed Data, Fitted Nonlinear Curves, and Dashed Linear Regression Lines ####
  # Get nonlinear model predictions
  pred_list <- data_clean %>%
    group_by(species) %>%
    do({
      sp <- unique(.$species)
      x_seq <- seq(min(.$x, na.rm = TRUE), max(.$x, na.rm = TRUE), length.out = 100)
      sp_model <- nls_models[[as.character(sp)]]
      pred <- predict(sp_model, newdata = data.frame(x = x_seq))
      data.frame(x = x_seq, pred = pred)
    })
  pred_all <- bind_rows(pred_list)
  
  # Compute endpoints for median-to-maximum regression line
  linear_endpoints_max <- coef_df %>%
    mutate(x_upper = x_max,
           y_upper = a + b * exp(-c * x_max)) %>%
    select(species, x50, median_value, x_upper, y_upper)
  
  linear_lines_max <- linear_endpoints_max %>%
    pivot_longer(cols = c(x50, x_upper), names_to = "point", values_to = "x") %>%
    mutate(y = ifelse(point == "x50", median_value, y_upper))
  
  # Compute endpoints for median-to-minimum regression line
  linear_endpoints_min <- coef_df %>%
    mutate(x_lower = x_min,
           y_lower = a + b * exp(-c * x_min)) %>%
    select(species, x50, median_value, x_lower, y_lower)
  
  linear_lines_min <- linear_endpoints_min %>%
    pivot_longer(cols = c(x50, x_lower), names_to = "point", values_to = "x") %>%
    mutate(y = ifelse(point == "x50", median_value, y_lower))
  
  p_combined <- ggplot() +
    geom_point(data = data_clean, aes(x = x, y = avg_value, color = species)) +
    geom_line(data = pred_all, aes(x = x, y = pred, color = species), size = 1) +
    # Add dashed regression lines from median to maximum and median to minimum
    geom_line(data = linear_lines_max, aes(x = x, y = y, group = species), 
              color = "black", size = 0.5, linetype = "dashed") +
    geom_line(data = linear_lines_min, aes(x = x, y = y, group = species), 
              color = "black", size = 0.5, linetype = "dashed") +
    # Add a larger black symbol at the median point (x50, median NDVI)
    geom_point(data = coef_df, aes(x = x50, y = median_value), 
               shape = 10, size = 5, color = "black") +
    scale_color_manual(values = cb_palette) +
    labs(x = "Transpiration Deficit", y = "NDVI Quantiles") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    ) +
    labs(caption = "(a)") +
    theme(plot.caption = element_text(face = "bold", size = 14, hjust = 0),
          plot.caption.position = "plot") +
    coord_cartesian(clip = "off")
  
  #### PANEL B: Bar Plot of x50 Values ####
  p_x50 <- ggplot(coef_df, aes(x = species, y = x50, fill = species)) +
    geom_col(width = 0.7) +
    labs(x = "Species", y = "Transpiration Deficit") +
    scale_fill_manual(values = cb_palette) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(caption = "(b)") +
    theme(plot.caption = element_text(face = "bold", size = 14, hjust = 0),
          plot.caption.position = "plot")
  
  #### PANEL C: Combined Vertical Bar Plot for Slopes (Both Regressions) ####
  # Compute slope statistics for median-to-maximum regression
  stats_list_max <- lapply(levels(data_clean$species), function(sp) {
    sp_coef <- coef_df %>% filter(species == sp)
    x50_val <- sp_coef$x50
    x_target <- sp_coef$x_max
    if(is.na(x50_val) | is.na(x_target)) {
      return(tibble(species = sp,
                    slope_lm = NA,
                    se = NA,
                    p_val = NA))
    }
    df_sp <- data_clean %>% filter(species == sp)
    lower_bound <- min(x50_val, x_target)
    upper_bound <- max(x50_val, x_target)
    df_subset <- df_sp %>% filter(x >= lower_bound, x <= upper_bound)
    if(nrow(df_subset) < 2) {
      x_points <- c(x50_val, x_target)
      y_points <- c(sp_coef$median_value, sp_coef$a + sp_coef$b * exp(-sp_coef$c * x_target))
      lm_fit <- lm(y_points ~ x_points)
      slope_lm <- coef(lm_fit)[["x_points"]]
      se <- summary(lm_fit)$coefficients[["x_points", "Std. Error"]]
      p_val <- summary(lm_fit)$coefficients[["x_points", "Pr(>|t|)"]]
    } else {
      lm_fit <- lm(avg_value ~ x, data = df_subset)
      slope_lm <- coef(lm_fit)[["x"]]
      se <- summary(lm_fit)$coefficients[["x", "Std. Error"]]
      p_val <- summary(lm_fit)$coefficients[["x", "Pr(>|t|)"]]
    }
    tibble(species = sp,
           slope_lm = slope_lm,
           se = se,
           p_val = p_val)
  })
  stats_df_max <- bind_rows(stats_list_max)
  stats_df_max$species <- factor(stats_df_max$species, levels = species_order)
  stats_df_max <- stats_df_max %>%
    mutate(RegType = "Max",
           label_text = ifelse(p_val < 0.05, "*", sprintf("%.3f", p_val)))
  
  # Compute slope statistics for median-to-minimum regression
  stats_list_min <- lapply(levels(data_clean$species), function(sp) {
    sp_coef <- coef_df %>% filter(species == sp)
    x50_val <- sp_coef$x50
    x_target <- sp_coef$x_min
    if(is.na(x50_val) | is.na(x_target)) {
      return(tibble(species = sp,
                    slope_lm = NA,
                    se = NA,
                    p_val = NA))
    }
    df_sp <- data_clean %>% filter(species == sp)
    lower_bound <- min(x50_val, x_target)
    upper_bound <- max(x50_val, x_target)
    df_subset <- df_sp %>% filter(x >= lower_bound, x <= upper_bound)
    if(nrow(df_subset) < 2) {
      x_points <- c(x50_val, x_target)
      y_points <- c(sp_coef$median_value, sp_coef$a + sp_coef$b * exp(-sp_coef$c * x_target))
      lm_fit <- lm(y_points ~ x_points)
      slope_lm <- coef(lm_fit)[["x_points"]]
      se <- summary(lm_fit)$coefficients[["x_points", "Std. Error"]]
      p_val <- summary(lm_fit)$coefficients[["x_points", "Pr(>|t|)"]]
    } else {
      lm_fit <- lm(avg_value ~ x, data = df_subset)
      slope_lm <- coef(lm_fit)[["x"]]
      se <- summary(lm_fit)$coefficients[["x", "Std. Error"]]
      p_val <- summary(lm_fit)$coefficients[["x", "Pr(>|t|)"]]
    }
    tibble(species = sp,
           slope_lm = slope_lm,
           se = se,
           p_val = p_val)
  })
  stats_df_min <- bind_rows(stats_list_min)
  stats_df_min$species <- factor(stats_df_min$species, levels = species_order)
  stats_df_min <- stats_df_min %>%
    mutate(RegType = "Min",
           label_text = ifelse(p_val < 0.05, "*", sprintf("%.3f", p_val)))
  
  # Combine the two sets of slope statistics
  stats_df_combined <- bind_rows(stats_df_max, stats_df_min)
  
  p_slope_combined <- ggplot(stats_df_combined, aes(x = species, y = slope_lm, fill = species, pattern = RegType)) +
    geom_col_pattern(position = position_dodge(width = 0.8), width = 0.8, color = "black",
                     pattern_fill = "black", pattern_angle = 45, pattern_density = 0.05, pattern_spacing = 0.1) +
    geom_errorbar(aes(ymin = slope_lm - se, ymax = slope_lm + se),
                  position = position_dodge(width = 0.8), width = 0.2) +
    geom_label(aes(label = label_text, y = slope_lm/2), 
               position = position_dodge(width = 0.8),
               fill = "white", alpha = 0.5, color = "black", size = 3) +
    scale_fill_manual(values = cb_palette) +
    scale_pattern_manual(values = c("Max" = "none", "Min" = "stripe")) +
    labs(x = "Species", y = "Slope", pattern = "Regression Type") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(caption = "(c)") +
    theme(plot.caption = element_text(face = "bold", size = 14, hjust = 0),
          plot.caption.position = "plot")
  
  #### Combine Panels ####
  # Arrange right side: Panel B on top and Panel C below.
  right_side <- p_x50 / p_slope_combined
  final_plot <- p_combined + right_side + plot_layout(widths = c(2, 1))
  print(final_plot)
  
  # Save the combined final plot
  ggsave(filename = output_figure, plot = final_plot, device = "png", width = 10, height = 8, dpi = 300)
  
  #### Coefficient Plot (a, b, and c) AS A SEPARATE FIGURE ####
  coeff_stats_list <- lapply(levels(data_clean$species), function(sp) {
    mod <- nls_models[[sp]]
    s <- summary(mod)$coefficients
    df <- as.data.frame(s)
    df$Coefficient <- rownames(df)
    df$species <- sp
    df
  })
  coeff_stats <- bind_rows(coeff_stats_list)
  coeff_stats <- coeff_stats %>%
    mutate(label = ifelse(`Pr(>|t|)` < 0.05,
                          "*",
                          sprintf("%.3f", `Pr(>|t|)`)))
  
  coeffs_long <- coef_df %>%
    select(species, a, b, c) %>%
    pivot_longer(cols = c("a", "b", "c"),
                 names_to = "Coefficient",
                 values_to = "Value")
  coeffs_long <- left_join(coeffs_long, coeff_stats %>% select(species, Coefficient, label),
                           by = c("species", "Coefficient"))
  coeffs_long$species <- factor(coeffs_long$species, levels = species_order)
  
  p_coeffs <- ggplot(coeffs_long, aes(x = species, y = Value, fill = species)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = label, y = Value/2), 
              color = "black", size = 3.5, 
              position = position_dodge(width = 0.9)) +
    scale_fill_manual(values = cb_palette) +
    facet_wrap(~ Coefficient, scales = "free_y") +
    labs(title = "Coefficients by Species for NDVI Quantiles",
         x = "Species",
         y = "Coefficient Value") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text = element_text(color = "black", size = 10),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "black"),
      axis.title = element_text(face = "bold", size = 12),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    )
  print(p_coeffs)
  
  # Save the coefficient figure to file
  dir.create(dirname(coef_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = coef_fig, plot = p_coeffs, device = "png", width = 10, height = 8, dpi = 300)
}

plot_Proportions_TDiff_slope_linear_each <- function(data, coef_fig, output_figure) {
  
  # Process data using the custom NDVI_TDiffbin function and remove missing values
  data <- NDVI_TDiffbin(data, 5)
  data <- na.omit(data)
  
  # Load required libraries
  library(ggplot2)
  library(nlme)    # for lmList
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(patchwork)
  library(purrr)
  library(ggpattern)  # for patterned bars in Panel C
  library(broom)      # to tidy lm objects
  
  # Identify the value column and set species order and color palette
  value_col <- "avg_value"
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  data$species <- factor(data$species, levels = species_order)
  
  cb_palette <- c("Oak"   = "#E69F00",   # Orange
                  "Beech" = "#0072B2",   # Deep blue
                  "Spruce"= "#009E73",   # Bluish-green
                  "Pine"  = "#F0E442")   # Yellow
  
  # For Transpiration Deficit, use x = bin_median directly
  data <- data %>% mutate(x = bin_median)
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x)) %>% droplevels()
  
  #### Compute species-specific summary statistics ####
  medians_df <- data_clean %>% 
    group_by(species) %>% 
    summarise(median_value = median(avg_value),
              x_median = median(x))
  
  max_x_df <- data_clean %>% 
    group_by(species) %>% 
    summarise(x_max = max(x))
  
  min_x_df <- data_clean %>% 
    group_by(species) %>% 
    summarise(x_min = min(x))
  
  #### LINEAR MODELING PER SPECIES (Overall) ####
  # Fit an overall linear model for each species
  lm_models <- lmList(avg_value ~ x | species, data = data_clean)
  print(summary(lm_models))
  
  # Extract overall coefficients: intercept (a) and slope (b)
  coef_df <- as.data.frame(coef(lm_models), optional = TRUE) %>%
    rownames_to_column(var = "species") %>%
    filter(!is.na(`(Intercept)`)) %>%
    rename(a = `(Intercept)`, b = x)
  coef_df$species <- factor(coef_df$species, levels = species_order)
  
  # Join with summary statistics
  coef_df <- coef_df %>% 
    left_join(medians_df, by = "species") %>%
    left_join(max_x_df, by = "species") %>%
    left_join(min_x_df, by = "species")
  
  # Calculate x50: x value where the overall predicted value equals the species-specific median
  coef_df <- coef_df %>%
    mutate(x50 = ifelse(b != 0, (median_value - a) / b, NA))
  
  #### PANEL A: Observed Data, Fitted Overall Linear Curves, and Dashed Segment Regressions ####
  # Get predictions for the overall linear model for each species
  pred_list <- data_clean %>%
    group_by(species) %>%
    do({
      sp <- unique(.$species)
      x_seq <- seq(min(.$x, na.rm = TRUE), max(.$x, na.rm = TRUE), length.out = 100)
      sp_coef <- coef_df %>% filter(species == sp)
      pred <- sp_coef$a + sp_coef$b * x_seq
      data.frame(x = x_seq, pred = pred)
    })
  pred_all <- bind_rows(pred_list)
  
  #### Segment Regressions: Fit separate linear models on each side of the median ####
  # Max side: for x >= x50
  stats_list_max <- lapply(levels(data_clean$species), function(sp) {
    df_sp <- data_clean %>% filter(species == sp)
    x50_val <- coef_df %>% filter(species == sp) %>% pull(x50)
    df_subset <- df_sp %>% filter(x >= x50_val)
    if(nrow(df_subset) < 2) {
      row <- coef_df %>% filter(species == sp)
      x_pts <- c(row$x50, row$x_max)
      y_pts <- c(row$median_value, row$a + row$b * row$x_max)
      lm_fit <- lm(y_pts ~ x_pts)
    } else {
      lm_fit <- lm(avg_value ~ x, data = df_subset)
    }
    tibble(species = sp,
           intercept_max = coef(lm_fit)[["(Intercept)"]],
           slope_max = coef(lm_fit)[["x"]],
           se_max = summary(lm_fit)$coefficients["x", "Std. Error"],
           p_val_max = summary(lm_fit)$coefficients["x", "Pr(>|t|)"])
  })
  stats_df_max <- bind_rows(stats_list_max)
  
  # Min side: for x <= x50
  stats_list_min <- lapply(levels(data_clean$species), function(sp) {
    df_sp <- data_clean %>% filter(species == sp)
    x50_val <- coef_df %>% filter(species == sp) %>% pull(x50)
    df_subset <- df_sp %>% filter(x <= x50_val)
    if(nrow(df_subset) < 2) {
      row <- coef_df %>% filter(species == sp)
      x_pts <- c(row$x50, row$x_min)
      y_pts <- c(row$median_value, row$a + row$b * row$x_min)
      lm_fit <- lm(y_pts ~ x_pts)
    } else {
      lm_fit <- lm(avg_value ~ x, data = df_subset)
    }
    tibble(species = sp,
           intercept_min = coef(lm_fit)[["(Intercept)"]],
           slope_min = coef(lm_fit)[["x"]],
           se_min = summary(lm_fit)$coefficients["x", "Std. Error"],
           p_val_min = summary(lm_fit)$coefficients["x", "Pr(>|t|)"])
  })
  stats_df_min <- bind_rows(stats_list_min)
  
  # Create data frames for the dashed line segments using the segment regression results
  seg_max_list <- lapply(levels(data_clean$species), function(sp) {
    row_coef <- coef_df %>% filter(species == sp)
    row_seg <- stats_df_max %>% filter(species == sp)
    x_seq <- seq(row_coef$x50, row_coef$x_max, length.out = 50)
    tibble(species = sp, x = x_seq, y = row_seg$intercept_max + row_seg$slope_max * x_seq)
  })
  seg_max_df <- bind_rows(seg_max_list)
  
  seg_min_list <- lapply(levels(data_clean$species), function(sp) {
    row_coef <- coef_df %>% filter(species == sp)
    row_seg <- stats_df_min %>% filter(species == sp)
    x_seq <- seq(row_coef$x_min, row_coef$x50, length.out = 50)
    tibble(species = sp, x = x_seq, y = row_seg$intercept_min + row_seg$slope_min * x_seq)
  })
  seg_min_df <- bind_rows(seg_min_list)
  
  p_combined <- ggplot() +
    geom_point(data = data_clean, aes(x = x, y = avg_value, color = species)) +
    geom_line(data = pred_all, aes(x = x, y = pred, color = species), size = 1) +
    # Add dashed segments from the segment regressions
    geom_line(data = seg_max_df, aes(x = x, y = y, group = species),
              color = "black", size = 0.5, linetype = "dashed") +
    geom_line(data = seg_min_df, aes(x = x, y = y, group = species),
              color = "black", size = 0.5, linetype = "dashed") +
    # Emphasize the median point with a triangle marker (shape = 17)
    geom_point(data = coef_df, aes(x = x50, y = median_value), 
               shape = 10, size = 5, color = "black") +
    scale_color_manual(values = cb_palette) +
    labs(x = "Transpiration Deficit", y = "NDVI Proportions") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.background = element_rect(fill = "white", color = "white"),
          panel.background = element_rect(fill = "white"),
          legend.background = element_rect(fill = "white", color = "white"),
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
          axis.title = element_text(face = "bold"),
          axis.text = element_text(color = "black"),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "top") +
    labs(caption = "(a)") +
    theme(plot.caption = element_text(face = "bold", size = 14, hjust = 0),
          plot.caption.position = "plot") +
    coord_cartesian(clip = "off")
  
  #### PANEL B: Bar Plot of x50 Values ####
  p_x50 <- ggplot(coef_df, aes(x = species, y = x50, fill = species)) +
    geom_col(width = 0.7) +
    labs(x = "Species", y = "Transpiration Deficit (x50)") +
    scale_fill_manual(values = cb_palette) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          axis.title.x = element_blank(),
          axis.title = element_text(face = "bold"),
          axis.text = element_text(color = "black"),
          plot.background = element_rect(fill = "white", color = "white"),
          panel.background = element_rect(fill = "white"),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(caption = "(b)") +
    theme(plot.caption = element_text(face = "bold", size = 14, hjust = 0),
          plot.caption.position = "plot")
  
  #### PANEL C: Combined Vertical Bar Plot for Slopes of Segment Regressions ####
  # Combine slope stats from the two segment regressions
  stats_df_combined <- bind_rows(
    stats_df_max %>% rename(slope = slope_max, se = se_max, p_val = p_val_max) %>% mutate(RegType = "Max"),
    stats_df_min %>% rename(slope = slope_min, se = se_min, p_val = p_val_min) %>% mutate(RegType = "Min")
  ) %>%
    mutate(label_text = ifelse(p_val < 0.05, "*", sprintf("%.3f", p_val)))
  
  p_slope_combined <- ggplot(stats_df_combined, aes(x = species, y = slope, fill = species, pattern = RegType)) +
    geom_col_pattern(position = position_dodge(width = 0.8), width = 0.8, color = "black",
                     pattern_fill = "black", pattern_angle = 45, pattern_density = 0.05, pattern_spacing = 0.1) +
    geom_errorbar(aes(ymin = slope - se, ymax = slope + se),
                  position = position_dodge(width = 0.8), width = 0.2) +
    geom_label(aes(label = label_text, y = slope/2),
               position = position_dodge(width = 0.8),
               fill = "white", alpha = 0.5, color = "black", size = 3) +
    scale_fill_manual(values = cb_palette) +
    scale_pattern_manual(values = c("Max" = "none", "Min" = "stripe")) +
    labs(x = "Species", y = "Slope", pattern = "Regression Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          axis.title.x = element_blank(),
          axis.title = element_text(face = "bold"),
          axis.text = element_text(color = "black"),
          plot.background = element_rect(fill = "white", color = "white"),
          panel.background = element_rect(fill = "white"),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(caption = "(c)") +
    theme(plot.caption = element_text(face = "bold", size = 14, hjust = 0),
          plot.caption.position = "plot")
  
  #### Combine Panels ####
  right_side <- p_x50 / p_slope_combined
  final_plot <- p_combined + right_side + plot_layout(widths = c(2, 1))
  print(final_plot)
  ggsave(filename = output_figure, plot = final_plot, device = "png", width = 10, height = 8, dpi = 300)
  
  #### Coefficient Plot (Intercept and Slope) AS A SEPARATE FIGURE ####
  # Use broom to tidy the individual lm models so we can extract p-values.
  tidy_lm <- map_df(lm_models, broom::tidy, .id = "species") %>%
    mutate(species = factor(species, levels = species_order),
           Coefficient = if_else(term == "(Intercept)", "a", "b"),
           label = if_else(p.value < 0.05, "*", sprintf("%.3f", p.value)))
  
  # Create a coefficient dataframe for plotting
  coeffs_long <- coef_df %>%
    select(species, a, b) %>%
    pivot_longer(cols = c("a", "b"), names_to = "Coefficient", values_to = "Value")
  coeffs_long$species <- factor(coeffs_long$species, levels = species_order)
  
  # Join with tidy lm to get p-value labels
  coeffs_long <- left_join(coeffs_long, tidy_lm %>% select(species, Coefficient, label),
                           by = c("species", "Coefficient"))
  
  p_coeffs <- ggplot(coeffs_long, aes(x = species, y = Value, fill = species)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = label, y = Value/2),
              color = "black", size = 3.5, position = position_dodge(width = 0.9)) +
    scale_fill_manual(values = cb_palette) +
    facet_wrap(~ Coefficient, scales = "free_y") +
    labs(title = "Coefficients by Species for NDVI Proportions (Linear Model)",
         x = "Species", y = "Coefficient Value") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text = element_text(color = "black", size = 10),
          plot.background = element_rect(fill = "white", color = "white"),
          panel.background = element_rect(fill = "white"),
          legend.background = element_rect(fill = "white", color = "white"),
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "black"),
          axis.title = element_text(face = "bold", size = 12),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "top",
          strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
          strip.text = element_text(face = "bold", size = 12))
  print(p_coeffs)
  
  dir.create(dirname(coef_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = coef_fig, plot = p_coeffs, device = "png", width = 10, height = 8, dpi = 300)
}

