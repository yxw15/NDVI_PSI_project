library(ncdf4)
library(reshape2)
library(ggplot2)
library(terra)
library(tidyverse)
library(dplyr)

process_NDVI_PSI_TDiff_species_year <- function(NDVI_file, species_name, mask_file, year, output_path) {
  # Load NDVI raster
  NDVI_raster <- rast(NDVI_file)
  
  # Load species mask
  species_mask <- rast(mask_file)
  
  # Extract the NDVI layer for the specific year 
  NDVI_year <- NDVI_raster[[year - 2003 + 1]]
  
  # Load and project PSI raster for the specific year
  PSI_folder <- output_path
  PSI_year <- rast(sprintf("%s/psi_%d.tif", PSI_folder, year))
  PSI_year <- project(PSI_year, species_mask)
  
  # Load and project TDiff raster for the specific year
  TDiff_folder <- output_path
  TDiff_year <- rast(sprintf("%s/tdiff_%d.tif", TDiff_folder, year))
  TDiff_year <- project(TDiff_year, species_mask)
  
  # Mask NDVI and PSI and TDiff with the species mask
  NDVI_masked <- mask(NDVI_year, species_mask)
  PSI_masked <- mask(PSI_year, species_mask)
  TDiff_masked <- mask(TDiff_year, species_mask)
  
  # Set NA values where species mask is NA
  NDVI_masked[is.na(species_mask[])] <- NA
  PSI_masked[is.na(species_mask[])] <- NA
  TDiff_masked[is.na(species_mask[])] <- NA
  
  # Save the rasters in their respective folders
  writeRaster(NDVI_masked, sprintf("%s/NDVI_%d_mask.tif", output_path, year), overwrite = TRUE)
  writeRaster(PSI_masked, sprintf("%s/PSI_%d_mask.tif", output_path, year), overwrite = TRUE)
  writeRaster(TDiff_masked, sprintf("%s/TDiff_%d_mask.tif", output_path, year), overwrite = TRUE)
  
  # Combine NDVI and PSI and TDiff into a data frame
  NDVI_PSI_TDiff <- c(NDVI_masked, PSI_masked, TDiff_masked)
  NDVI_PSI_TDiff_df <- as.data.frame(NDVI_PSI_TDiff, xy = TRUE, na.rm = TRUE)
  NDVI_PSI_TDiff_df$year <- year
  NDVI_PSI_TDiff_df$species <- species_name
  
  # Rename 'Quantiles_1' to 'Quantiles_22' to 'Quantiles'
  quantile_columns <- grep("^Quantiles_\\d+$", names(NDVI_PSI_TDiff_df), value = TRUE)
  if (length(quantile_columns) > 0) {
    names(NDVI_PSI_TDiff_df)[names(NDVI_PSI_TDiff_df) %in% quantile_columns] <- "Quantiles"
  }
  
  return(NDVI_PSI_TDiff_df)
}

NDVI_PSI_bin_mean <- function(df_all) {
  # Determine the correct column name
  value_column <- if ("Quantiles" %in% names(df_all)) {
    "Quantiles"
  } else if ("Proportions" %in% names(df_all)) {
    "Proportions"
  } else {
    stop("Neither Quantiles nor Proportions column found.")
  }
  
  # Dynamically name the output column
  mean_column_name <- if (value_column == "Quantiles") "mean_quantiles" else "mean_proportions"
  
  # Bin the soil water potential and summarize mean quantiles/proportions by year and PSI bin
  df_summary <- df_all %>%
    mutate(PSI_bin = floor(soil_water_potential / 100) * 100) %>%
    group_by(year, PSI_bin, species) %>%
    summarise(
      !!mean_column_name := mean(.data[[value_column]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    ungroup()
  
  return(df_summary)
}

NDVI_PSI_bin_median <- function(df_all) {
  # Determine the correct column name
  value_column <- if ("Quantiles" %in% names(df_all)) {
    "Quantiles"
  } else if ("Proportions" %in% names(df_all)) {
    "Proportions"
  } else {
    stop("Neither Quantiles nor Proportions column found.")
  }
  
  # Dynamically name the output column
  median_column_name <- if (value_column == "Quantiles") "median_quantiles" else "median_proportions"
  
  # Bin the soil water potential and summarize median quantiles/proportions by year and PSI bin
  df_summary <- df_all %>%
    mutate(PSI_bin = floor(soil_water_potential / 100) * 100) %>%
    group_by(year, PSI_bin, species) %>%
    summarise(
      !!median_column_name := median(.data[[value_column]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    ungroup()
  
  return(df_summary)
}

plot_NDVI_PSI_species_mean <- function(df_summary, filename) {
  # Load required library
  library(ggplot2)
  
  # Determine the correct y-axis column name
  y_column <- if ("mean_quantiles" %in% names(df_summary)) {
    "mean_quantiles"
  } else if ("mean_proportions" %in% names(df_summary)) {
    "mean_proportions"
  } else {
    stop("Neither mean_quantiles nor mean_proportions column found in df_summary.")
  }
  
  # Set y-axis label dynamically
  y_label <- if (y_column == "mean_quantiles") "Mean Quantiles" else "Mean Proportions"
  
  # Set plot title dynamically
  plot_title <- paste(y_label, "vs PSI by Species (Binned by 100)")
  
  # Create the plot
  plot <- ggplot(df_summary, aes(x = PSI_bin, y = .data[[y_column]], color = as.factor(year), group = year)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_line(size = 1) +
    labs(
      x = "PSI (binned)",
      y = y_label,
      color = "Year",
      title = plot_title
    ) +
    ylim(1, 22) +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = "black"),
      legend.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.position = "top",
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    ) +
    scale_color_viridis_d() +
    facet_wrap(~ factor(species, levels = c("Oak", "Beech", "Spruce", "Pine")), ncol = 2)
  
  # Save the plot
  ggsave(filename, plot = plot, width = 14, height = 10, dpi = 300)
}

plot_NDVI_PSI_species_median <- function(df_summary, filename) {
  # Load required library
  library(ggplot2)
  
  # Determine the correct y-axis column name
  y_column <- if ("median_quantiles" %in% names(df_summary)) {
    "median_quantiles"
  } else if ("median_proportions" %in% names(df_summary)) {
    "median_proportions"
  } else {
    stop("Neither median_quantiles nor median_proportions column found in df_summary.")
  }
  
  # Set y-axis label dynamically
  y_label <- if (y_column == "median_quantiles") "Median Quantiles" else "Median Proportions"
  
  # Set plot title dynamically
  plot_title <- paste(y_label, "vs PSI by Species (Binned by 100)")
  
  # Create the plot
  plot <- ggplot(df_summary, aes(x = PSI_bin, y = .data[[y_column]], color = as.factor(year), group = year)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_line(size = 1) +
    labs(
      x = "PSI (binned)",
      y = y_label,
      color = "Year",
      title = plot_title
    ) +
    ylim(1, 22) +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = "black"),
      legend.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.position = "top",
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    ) +
    scale_color_viridis_d() +
    facet_wrap(~ factor(species, levels = c("Oak", "Beech", "Spruce", "Pine")), ncol = 2) # Create subplots for each species with 2 columns
  
  # Save the plot
  ggsave(filename, plot = plot, width = 14, height = 10, dpi = 300)
}

df_average_year <- function(df_all) {
  # Determine the correct column name for quantile/proportion
  value_column <- if ("Quantiles" %in% names(df_all)) {
    "Quantiles"
  } else if ("Proportions" %in% names(df_all)) {
    "Proportions"
  } else {
    stop("Neither Quantiles nor Proportions column found in df_all.")
  }
  
  # Dynamically name the output column
  avg_column_name <- if (value_column == "Quantiles") "avg_quantile" else "avg_proportion"
  
  # Calculate average quantile/proportion and soil water potential (PSI) per year
  df_avg <- df_all %>%
    group_by(year) %>%
    summarise(
      !!avg_column_name := mean(.data[[value_column]], na.rm = TRUE),
      avg_psi = mean(soil_water_potential, na.rm = TRUE),
      avg_tdiff = mean(transpiration_deficit, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(df_avg)
}

df_average_year_NDVI_TDiff <- function(df_all) {
  # Determine the correct column name for quantile/proportion
  value_column <- if ("Quantiles" %in% names(df_all)) {
    "Quantiles"
  } else if ("Proportions" %in% names(df_all)) {
    "Proportions"
  } else {
    stop("Neither Quantiles nor Proportions column found in df_all.")
  }
  
  # Dynamically name the output column
  avg_column_name <- if (value_column == "Quantiles") "avg_quantile" else "avg_proportion"
  
  # Calculate average quantile/proportion and soil water potential (PSI) per year
  df_avg <- df_all %>%
    group_by(year) %>%
    summarise(
      !!avg_column_name := mean(.data[[value_column]], na.rm = TRUE),
      avg_tdiff = mean(transpiration_deficit, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(df_avg)
}

df_average_year_species <- function(df_all) {
  # Determine the correct column name for quantile/proportion
  value_column <- if ("Quantiles" %in% names(df_all)) {
    "Quantiles"
  } else if ("Proportions" %in% names(df_all)) {
    "Proportions"
  } else {
    stop("Neither Quantiles nor Proportions column found in df_all.")
  }
  
  # Dynamically name the output column
  avg_column_name <- if (value_column == "Quantiles") "avg_quantile" else "avg_proportion"
  
  # Calculate average quantile/proportion and soil water potential (PSI) per year and species
  df_avg <- df_all %>%
    group_by(year, species) %>%
    summarise(
      !!avg_column_name := mean(.data[[value_column]], na.rm = TRUE),
      avg_psi = mean(soil_water_potential, na.rm = TRUE),
      avg_tdiff = mean(transpiration_deficit, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(df_avg)
}

df_average_year_species_NDVI_TDiff <- function(df_all) {
  # Determine the correct column name for quantile/proportion
  value_column <- if ("Quantiles" %in% names(df_all)) {
    "Quantiles"
  } else if ("Proportions" %in% names(df_all)) {
    "Proportions"
  } else {
    stop("Neither Quantiles nor Proportions column found in df_all.")
  }
  
  # Dynamically name the output column
  avg_column_name <- if (value_column == "Quantiles") "avg_quantile" else "avg_proportion"
  
  # Calculate average quantile/proportion and soil water potential (PSI) per year and species
  df_avg <- df_all %>%
    group_by(year, species) %>%
    summarise(
      !!avg_column_name := mean(.data[[value_column]], na.rm = TRUE),
      avg_tdiff = mean(transpiration_deficit, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(df_avg)
}

plot_correlation_NDVI_PSI_avg <- function(df_all, plot_corr_path) {
  # Load necessary libraries
  library(ggplot2)
  
  df_avg <- df_average_year(df_all)
  
  # Determine the correct NDVI column
  ndvi_column <- if ("avg_quantile" %in% names(df_avg)) {
    "avg_quantile"
  } else if ("avg_proportion" %in% names(df_avg)) {
    "avg_proportion"
  } else {
    stop("Neither avg_quantile nor avg_proportion column found in df_avg.")
  }
  
  # Set NDVI y-axis label dynamically
  ndvi_label <- if (ndvi_column == "avg_quantile") {
    "Average NDVI (Quantiles)"
  } else {
    "Average NDVI (Proportions)"
  }
  
  # Set dynamic plot title based on the NDVI data type
  plot_title <- if (ndvi_column == "avg_quantile") {
    "Correlation between NDVI (Quantiles) and PSI"
  } else {
    "Correlation between NDVI (Proportions) and PSI"
  }
  
  # Calculate correlation coefficient
  correlation <- cor(df_avg$avg_psi, df_avg[[ndvi_column]], use = "complete.obs")
  
  # Create correlation plot
  plot_corr <- ggplot(df_avg, aes(x = avg_psi, y = .data[[ndvi_column]])) +
    geom_point(color = "purple", size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", color = "black", se = TRUE) +
    annotate(
      "text", 
      x = min(df_avg$avg_psi, na.rm = TRUE), 
      y = max(df_avg[[ndvi_column]], na.rm = TRUE), 
      label = paste("Correlation:", round(correlation, 2)), 
      hjust = 0, vjust = 1.5, size = 5, color = "blue"
    ) +
    labs(
      title = plot_title,
      x = "Average PSI",
      y = ndvi_label
    ) +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = "black"),
      legend.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.position = "top",
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    ) +
    scale_color_viridis_d()
  
  # Save the plot
  ggsave(plot_corr_path, plot = plot_corr, width = 10, height = 8, dpi = 300)
}

plot_correlation_NDVI_TDiff_avg <- function(df_all, plot_corr_path) {
  # Load necessary libraries
  library(ggplot2)
  
  df_avg <- df_average_year_NDVI_TDiff(df_all)
  
  # Determine the correct NDVI column
  ndvi_column <- if ("avg_quantile" %in% names(df_avg)) {
    "avg_quantile"
  } else if ("avg_proportion" %in% names(df_avg)) {
    "avg_proportion"
  } else {
    stop("Neither avg_quantile nor avg_proportion column found in df_avg.")
  }
  
  # Set NDVI y-axis label dynamically
  ndvi_label <- if (ndvi_column == "avg_quantile") {
    "Average NDVI (Quantiles)"
  } else {
    "Average NDVI (Proportions)"
  }
  
  # Set dynamic plot title based on the NDVI data type
  plot_title <- if (ndvi_column == "avg_quantile") {
    "Correlation between NDVI (Quantiles) and TDiff"
  } else {
    "Correlation between NDVI (Proportions) and TDiff"
  }
  
  # Calculate correlation coefficient
  correlation <- cor(df_avg$avg_tdiff, df_avg[[ndvi_column]], use = "complete.obs")
  
  # Create correlation plot
  plot_corr <- ggplot(df_avg, aes(x = avg_tdiff, y = .data[[ndvi_column]])) +
    geom_point(color = "green4", size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", color = "black", se = TRUE) +
    annotate(
      "text", 
      x = min(df_avg$avg_tdiff, na.rm = TRUE), 
      y = max(df_avg[[ndvi_column]], na.rm = TRUE), 
      label = paste("Correlation:", round(correlation, 2)), 
      hjust = 0, vjust = 1.5, size = 5, color = "blue"
    ) +
    labs(
      title = plot_title,
      x = "Average TDiff",
      y = ndvi_label
    ) +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = "black"),
      legend.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.position = "top",
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    ) +
    scale_color_viridis_d()
  
  # Save the plot
  ggsave(plot_corr_path, plot = plot_corr, width = 10, height = 8, dpi = 300)
}

plot_series_NDVI_PSI_same_month <- function(df_all, plot_series_path) {
  # Load necessary libraries
  library(ggplot2)
  library(dplyr)
  library(patchwork) # For combining plots into subplots
  
  df_avg <- df_average_year(df_all)
  
  # Determine the correct NDVI column
  ndvi_column <- if ("avg_quantile" %in% names(df_avg)) {
    "avg_quantile"
  } else if ("avg_proportion" %in% names(df_avg)) {
    "avg_proportion"
  } else {
    stop("Neither avg_quantile nor avg_proportion column found in df_avg.")
  }
  
  # Set NDVI y-axis label dynamically
  ndvi_label <- if (ndvi_column == "avg_quantile") {
    "Average NDVI (Quantiles)"
  } else {
    "Average NDVI (Proportions)"
  }
  
  # Set dynamic NDVI plot title based on the NDVI data type
  ndvi_title <- if (ndvi_column == "avg_quantile") {
    "(a) NDVI (Quantiles) Over Time"
  } else {
    "(a) NDVI (Proportions) Over Time"
  }
  
  # Create NDVI plot with group aesthetic inside aes()
  plot_ndvi <- ggplot(df_avg, aes(x = year, y = .data[[ndvi_column]], group = 1)) +
    geom_line(color = "blue", size = 1) +
    geom_point(color = "blue", size = 2) +
    labs(
      title = ndvi_title,
      x = "Year",
      y = ndvi_label
    ) +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = "black"),
      legend.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.position = "top",
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    )
  
  # Create PSI plot with group aesthetic inside aes()
  plot_psi <- ggplot(df_avg, aes(x = year, y = avg_psi, group = 1)) +
    geom_line(color = "red", size = 1) +
    geom_point(color = "red", size = 2) +
    labs(
      title = "(b) PSI Over Time",
      x = "Year",
      y = "Average PSI"
    ) +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = "black"),
      legend.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.position = "top",
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    )
  
  # Combine plots using patchwork
  combined_plot <- plot_ndvi / plot_psi
  
  # Save the combined plot to the specified path
  ggsave(plot_series_path, plot = combined_plot, width = 12, height = 10, dpi = 300)
}

plot_series_NDVI_TDiff_same_month <- function(df_all, plot_series_path) {
  # Load necessary libraries
  library(ggplot2)
  library(dplyr)
  library(patchwork) # For combining plots into subplots
  
  df_avg <- df_average_year_NDVI_TDiff(df_all)
  df_avg$year <- as.numeric(as.character(df_avg$year))
  
  # Determine the correct NDVI column
  ndvi_column <- if ("avg_quantile" %in% names(df_avg)) {
    "avg_quantile"
  } else if ("avg_proportion" %in% names(df_avg)) {
    "avg_proportion"
  } else {
    stop("Neither avg_quantile nor avg_proportion column found in df_avg.")
  }
  
  # Set NDVI y-axis label dynamically
  ndvi_label <- if (ndvi_column == "avg_quantile") {
    "Average NDVI (Quantiles)"
  } else {
    "Average NDVI (Proportions)"
  }
  
  # Set dynamic NDVI plot title based on the NDVI data type
  ndvi_title <- if (ndvi_column == "avg_quantile") {
    "(a) NDVI (Quantiles) Over Time"
  } else {
    "(a) NDVI (Proportions) Over Time"
  }
  
  # Create NDVI plot with group aesthetic inside aes()
  plot_ndvi <- ggplot(df_avg, aes(x = year, y = .data[[ndvi_column]], group = 1)) +
    geom_line(color = "blue", linewidth = 1) +
    geom_point(color = "blue", size = 2) +
    scale_x_continuous(breaks = seq(min(df_avg$year), max(df_avg$year), by = 5)) +
    labs(
      title = ndvi_title,
      x = "Year",
      y = ndvi_label
    ) +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = "black"),
      legend.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.position = "top",
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    )
  
  # Create PSI plot with group aesthetic inside aes()
  plot_tidff <- ggplot(df_avg, aes(x = year, y = avg_tdiff, group = 1)) +
    geom_line(color = "red", size = 1) +
    geom_point(color = "red", size = 2) +
    labs(
      title = "(b) TDiff Over Time",
      x = "Year",
      y = "Average TDiff"
    ) +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = "black"),
      legend.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.position = "top",
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    )
  
  # Combine plots using patchwork
  combined_plot <- plot_ndvi / plot_tidff
  
  # Save the combined plot to the specified path
  ggsave(plot_series_path, plot = combined_plot, width = 12, height = 10, dpi = 300)
}

plot_series_NDVI_PSI_same_month_species <- function(df_all, plot_species_path) {
  # Load necessary libraries
  library(ggplot2)
  library(dplyr)
  library(patchwork) # For combining plots into subplots
  
  # Determine the correct NDVI column
  ndvi_column <- if ("Quantiles" %in% names(df_all)) {
    "Quantiles"
  } else if ("Proportions" %in% names(df_all)) {
    "Proportions"
  } else {
    stop("Neither Quantiles nor Proportions column found in df_all.")
  }
  
  # Define custom color palette for species
  cb_palette <- c("Oak"    = "#E69F00",  # Orange
                  "Beech"  = "#0072B2",  # Deep blue
                  "Spruce" = "#009E73",  # Bluish-green
                  "Pine"   = "#F0E442")  # Yellow
  
  # Dynamically set column names for output
  avg_ndvi_column <- if (ndvi_column == "Quantiles") "avg_quantile" else "avg_proportion"
  
  # Calculate average NDVI and PSI per species and year
  df_species_avg <- df_all %>%
    group_by(year, species) %>%
    summarise(
      !!avg_ndvi_column := mean(.data[[ndvi_column]], na.rm = TRUE),
      avg_psi = mean(soil_water_potential, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Convert year to numeric to ensure a continuous scale
  df_species_avg <- df_species_avg %>%
    mutate(year = as.numeric(as.character(year)))
  
  # Define species order
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  df_species_avg <- df_species_avg %>%
    mutate(species = factor(species, levels = species_order))
  
  # Set dynamic NDVI label and title
  ndvi_label <- if (ndvi_column == "Quantiles") "Average NDVI (Quantiles)" else "Average NDVI (Proportions)"
  ndvi_title <- if (ndvi_column == "Quantiles") "(a) NDVI (Quantiles) by Species" else "(a) NDVI (Proportions) by Species"
  
  # Determine x-axis breaks (every 5 years)
  x_breaks <- seq(min(df_species_avg$year), max(df_species_avg$year), by = 5)
  
  # Create NDVI plot
  plot_ndvi <- ggplot(df_species_avg, aes(x = year, y = .data[[avg_ndvi_column]], color = species, group = species)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    labs(
      title = ndvi_title,
      x = "Year",
      y = ndvi_label
    ) +
    scale_x_continuous(breaks = x_breaks) +
    facet_wrap(~ species, ncol = 2) +
    theme_minimal() +
    theme(
      plot.background    = element_rect(fill = "white", color = NA),
      panel.background   = element_rect(fill = "white", color = "black"),
      legend.background  = element_rect(fill = "white", color = NA),
      plot.title         = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title         = element_text(face = "bold"),
      axis.text          = element_text(color = "black"),
      panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.position    = "top",
      strip.background   = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text         = element_text(face = "bold", size = 12)
    ) +
    scale_color_manual(values = cb_palette)
  
  # Create PSI plot
  plot_psi <- ggplot(df_species_avg, aes(x = year, y = avg_psi, color = species, group = species)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    labs(
      title = "(b) PSI by Species",
      x = "Year",
      y = "Average PSI"
    ) +
    scale_x_continuous(breaks = x_breaks) +
    facet_wrap(~ species, ncol = 2) +
    theme_minimal() +
    theme(
      plot.background    = element_rect(fill = "white", color = NA),
      panel.background   = element_rect(fill = "white", color = "black"),
      legend.background  = element_rect(fill = "white", color = NA),
      plot.title         = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title         = element_text(face = "bold"),
      axis.text          = element_text(color = "black"),
      panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.position    = "top",
      strip.background   = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text         = element_text(face = "bold", size = 12)
    ) +
    scale_color_manual(values = cb_palette)
  
  # Combine plots using patchwork
  combined_species_plot <- plot_ndvi / plot_psi
  
  # Save the combined plot to the specified path
  ggsave(plot_species_path, plot = combined_species_plot, width = 14, height = 12, dpi = 300)
}

plot_series_NDVI_TDiff_same_month_species <- function(df_all, plot_species_path) {
  # Load necessary libraries
  library(ggplot2)
  library(dplyr)
  library(patchwork) # For combining plots into subplots
  
  # Determine the correct NDVI column
  ndvi_column <- if ("Quantiles" %in% names(df_all)) {
    "Quantiles"
  } else if ("Proportions" %in% names(df_all)) {
    "Proportions"
  } else {
    stop("Neither Quantiles nor Proportions column found in df_all.")
  }
  
  # Define custom color palette for species
  cb_palette <- c("Oak"    = "#E69F00",  # Orange
                  "Beech"  = "#0072B2",  # Deep blue
                  "Spruce" = "#009E73",  # Bluish-green
                  "Pine"   = "#F0E442")  # Yellow
  
  # Dynamically set column names for output
  avg_ndvi_column <- if (ndvi_column == "Quantiles") "avg_quantile" else "avg_proportion"
  
  # Calculate average NDVI and TDiff per species and year
  df_species_avg <- df_all %>%
    group_by(year, species) %>%
    summarise(
      !!avg_ndvi_column := mean(.data[[ndvi_column]], na.rm = TRUE),
      avg_tdiff = mean(transpiration_deficit, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Convert year to numeric to ensure a continuous scale
  df_species_avg <- df_species_avg %>%
    mutate(year = as.numeric(as.character(year)))
  
  # Define species order
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  df_species_avg <- df_species_avg %>%
    mutate(species = factor(species, levels = species_order))
  
  # Set dynamic NDVI label and title
  ndvi_label <- if (ndvi_column == "Quantiles") "Average NDVI (Quantiles)" else "Average NDVI (Proportions)"
  ndvi_title <- if (ndvi_column == "Quantiles") "(a) NDVI (Quantiles) by Species" else "(a) NDVI (Proportions) by Species"
  
  # Determine x-axis breaks (every 5 years)
  x_breaks <- seq(min(df_species_avg$year), max(df_species_avg$year), by = 5)
  
  # Create NDVI plot
  plot_ndvi <- ggplot(df_species_avg, aes(x = year, y = .data[[avg_ndvi_column]], color = species, group = species)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    labs(
      title = ndvi_title,
      x = "Year",
      y = ndvi_label
    ) +
    scale_x_continuous(breaks = x_breaks) +
    facet_wrap(~ species, ncol = 2) +
    theme_minimal() +
    theme(
      plot.background    = element_rect(fill = "white", color = NA),
      panel.background   = element_rect(fill = "white", color = "black"),
      legend.background  = element_rect(fill = "white", color = NA),
      plot.title         = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title         = element_text(face = "bold"),
      axis.text          = element_text(color = "black"),
      panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.position    = "top",
      strip.background   = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text         = element_text(face = "bold", size = 12)
    ) +
    scale_color_manual(values = cb_palette)
  
  # Create TDiff plot
  plot_tdiff <- ggplot(df_species_avg, aes(x = year, y = avg_tdiff, color = species, group = species)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    labs(
      title = "(b) TDiff by Species",
      x = "Year",
      y = "Average TDiff"
    ) +
    scale_x_continuous(breaks = x_breaks) +
    facet_wrap(~ species, ncol = 2) +
    theme_minimal() +
    theme(
      plot.background    = element_rect(fill = "white", color = NA),
      panel.background   = element_rect(fill = "white", color = "black"),
      legend.background  = element_rect(fill = "white", color = NA),
      plot.title         = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title         = element_text(face = "bold"),
      axis.text          = element_text(color = "black"),
      panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.position    = "top",
      strip.background   = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text         = element_text(face = "bold", size = 12)
    ) +
    scale_color_manual(values = cb_palette)
  
  # Combine plots using patchwork
  combined_species_plot <- plot_ndvi / plot_tdiff
  
  # Save the combined plot to the specified path
  ggsave(plot_species_path, plot = combined_species_plot, width = 14, height = 12, dpi = 300)
}

plot_correlation_NDVI_PSI_species_avg <- function(df_all, plot_corr_path) {
  # Load necessary libraries
  library(ggplot2)
  library(dplyr)
  
  df_avg_species <- df_average_year_species(df_all)
  
  # Determine the correct NDVI column
  ndvi_column <- if ("avg_quantile" %in% names(df_avg_species)) {
    "avg_quantile"
  } else if ("avg_proportion" %in% names(df_avg_species)) {
    "avg_proportion"
  } else {
    stop("Neither avg_quantile nor avg_proportion column found in df_avg_species.")
  }
  
  # Set NDVI y-axis label dynamically
  ndvi_label <- if (ndvi_column == "avg_quantile") "Average NDVI (Quantiles)" else "Average NDVI (Proportions)"
  
  # Set dynamic plot title based on the NDVI data type
  plot_title <- if (ndvi_column == "avg_quantile") {
    "Correlation between NDVI (Quantiles) and PSI by Species"
  } else {
    "Correlation between NDVI (Proportions) and PSI by Species"
  }
  
  # Calculate correlation coefficients for each species
  correlations <- df_avg_species %>%
    group_by(species) %>%
    summarize(correlation = cor(avg_psi, .data[[ndvi_column]], use = "complete.obs"))
  
  # Merge correlation data back to the original dataframe for annotation
  df_avg_species <- df_avg_species %>%
    left_join(correlations, by = "species")
  
  # Create the plot
  plot_corr <- ggplot(df_avg_species, aes(x = avg_psi, y = .data[[ndvi_column]])) +
    geom_point(color = "purple", size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", color = "black", se = TRUE) +
    facet_wrap(~ species, scales = "free") +
    geom_text(
      data = correlations,
      aes(x = -Inf, y = Inf, label = paste("Correlation:", round(correlation, 2))),
      hjust = -0.1, vjust = 1.5, size = 5, color = "blue",
      inherit.aes = FALSE
    ) +
    labs(
      title = plot_title,
      x = "Average Soil Water Potential",
      y = ndvi_label
    ) +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = "black"),
      legend.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.position = "top",
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    )
  
  # Save the plot
  ggsave(plot_corr_path, plot = plot_corr, width = 14, height = 10, dpi = 300)
}

plot_correlation_NDVI_TDiff_species_avg <- function(df_all, plot_corr_path) {
  # Load necessary libraries
  library(ggplot2)
  library(dplyr)
  
  df_avg_species <- df_average_year_species_NDVI_TDiff(df_all)
  
  # Determine the correct NDVI column
  ndvi_column <- if ("avg_quantile" %in% names(df_avg_species)) {
    "avg_quantile"
  } else if ("avg_proportion" %in% names(df_avg_species)) {
    "avg_proportion"
  } else {
    stop("Neither avg_quantile nor avg_proportion column found in df_avg_species.")
  }
  
  # Set NDVI y-axis label dynamically
  ndvi_label <- if (ndvi_column == "avg_quantile") "Average NDVI (Quantiles)" else "Average NDVI (Proportions)"
  
  # Set dynamic plot title based on the NDVI data type
  plot_title <- if (ndvi_column == "avg_quantile") {
    "Correlation between NDVI (Quantiles) and TDiff by Species"
  } else {
    "Correlation between NDVI (Proportions) and TDiff by Species"
  }
  
  # Calculate correlation coefficients for each species
  correlations <- df_avg_species %>%
    group_by(species) %>%
    summarize(correlation = cor(avg_tdiff, .data[[ndvi_column]], use = "complete.obs"))
  
  # Merge correlation data back to the original dataframe for annotation
  df_avg_species <- df_avg_species %>%
    left_join(correlations, by = "species")
  
  # Create the plot
  plot_corr <- ggplot(df_avg_species, aes(x = avg_tdiff, y = .data[[ndvi_column]])) +
    geom_point(color = "green4", size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", color = "black", se = TRUE) +
    facet_wrap(~ species, scales = "free") +
    geom_text(
      data = correlations,
      aes(x = -Inf, y = Inf, label = paste("Correlation:", round(correlation, 2))),
      hjust = -0.1, vjust = 1.5, size = 5, color = "blue",
      inherit.aes = FALSE
    ) +
    labs(
      title = plot_title,
      x = "Average Transpiration Deficit",
      y = ndvi_label
    ) +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = "black"),
      legend.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.position = "top",
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    )
  
  # Save the plot
  ggsave(plot_corr_path, plot = plot_corr, width = 14, height = 10, dpi = 300)
}

plot_correlation_NDVI_PSI_species <- function(df_all, plot_corr_path) {
  # Load necessary libraries
  library(ggplot2)
  library(dplyr)
  
  # Determine the correct NDVI column
  ndvi_column <- if ("Quantiles" %in% names(df_all)) {
    "Quantiles"
  } else if ("Proportions" %in% names(df_all)) {
    "Proportions"
  } else {
    stop("Neither Quantiles nor Proportions column found in df_all.")
  }
  
  # Set NDVI y-axis label dynamically
  ndvi_label <- if (ndvi_column == "Quantiles") "NDVI (Quantiles)" else "NDVI (Proportions)"
  
  # Calculate correlation coefficients for each species
  correlations <- df_all %>%
    group_by(species) %>%
    summarize(correlation = cor(.data[[ndvi_column]], soil_water_potential, use = "complete.obs"))
  
  # Merge correlation data back to the original dataframe for annotation
  df_all <- df_all %>%
    left_join(correlations, by = "species")
  
  # Create the plot
  plot_corr <- ggplot(df_all, aes(x = soil_water_potential, y = .data[[ndvi_column]])) +
    geom_point(color = "darkgreen", size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", color = "black", se = TRUE) +
    facet_wrap(~ species, scales = "free") +
    geom_text(
      data = correlations,
      aes(x = -Inf, y = Inf, label = paste("Correlation:", round(correlation, 2))),
      hjust = -0.1, vjust = 1.5, size = 5, color = "blue",
      inherit.aes = FALSE
    ) +
    labs(
      title = "Correlation between NDVI and Soil Water Potential by Species",
      x = "Soil Water Potential",
      y = ndvi_label
    ) +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = "black"),
      legend.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.position = "top",
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    )
  
  # Save the plot
  ggsave(plot_corr_path, plot = plot_corr, width = 14, height = 10, dpi = 300)
}

plot_density_NDVI_PSI_species <- function(data, plot_density_path) {
  library(ggplot2)
  library(viridis)
  
  # Determine the correct NDVI column
  ndvi_column <- if ("Quantiles" %in% names(data)) {
    "Quantiles"
  } else if ("Proportions" %in% names(data)) {
    "Proportions"
  } else {
    stop("Neither Quantiles nor Proportions column found in data.")
  }
  
  # Set NDVI y-axis label dynamically
  ndvi_label <- if (ndvi_column == "Quantiles") "NDVI (Quantiles)" else "NDVI (Proportions)"
  
  # Filter for specific species in the desired order
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  data <- data[data$species %in% species_order, ]
  data$species <- factor(data$species, levels = species_order)  # Ensure the correct order
  
  # Create the plot with faceting
  plot_density_species <- ggplot(data, aes(x = soil_water_potential, y = .data[[ndvi_column]])) +
    geom_density_2d_filled() +
    scale_fill_viridis_d(option = "C", direction = -1) +  # Colorblind-friendly palette
    labs(
      title = "Density Plot of NDVI vs. Soil Water Potential by Species",
      x = "Soil Water Potential",
      y = ndvi_label,
      fill = "Density"
    ) +
    theme_minimal(base_size = 14) +  # Base size for better readability
    theme(
      plot.background = element_rect(fill = "white", color = NA),    # Clean white background
      panel.background = element_rect(fill = "white", color = NA),  # Panel background
      legend.background = element_rect(fill = "white", color = NA), # Legend background
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold", color = "black"), # Centered bold title
      axis.title = element_text(face = "bold", size = 16),          # Bold axis titles
      axis.text = element_text(size = 12, color = "black"),         # Axis tick labels
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8), # Add subtle border
      legend.position = "top",                                      # Legend at the top
      legend.title = element_text(size = 14, face = "bold"),        # Bold legend title
      legend.text = element_text(size = 12),                       # Legend text size
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5), # Strip for facets
      strip.text = element_text(face = "bold", size = 14)          # Facet label styles
    ) +
    facet_wrap(~species, ncol = 2)  # Create subplots for each species in 2 columns
  
  # Save the plot to the specified path
  ggsave(filename = plot_density_path, plot = plot_density_species, width = 14, height = 10, dpi = 300)
  
  # Return the plot object for viewing
  return(plot_density_species)
}

plot_density_NDVI_PSI_species_individual <- function(data, plot_density_folder) {
  library(ggplot2)
  library(viridis)
  library(stringr)
  
  # Determine the correct NDVI column
  ndvi_column <- if ("Quantiles" %in% names(data)) {
    "Quantiles"
  } else if ("Proportions" %in% names(data)) {
    "Proportions"
  } else {
    stop("Neither Quantiles nor Proportions column found in data.")
  }
  
  # Set NDVI y-axis label dynamically
  ndvi_label <- if (ndvi_column == "Quantiles") "NDVI (Quantiles)" else "NDVI (Proportions)"
  
  # Filter for specific species in the desired order
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  data <- data[data$species %in% species_order, ]
  data$species <- factor(data$species, levels = species_order)  # Ensure the correct order
  
  # Create individual plots for each species
  for (sp in species_order) {
    species_data <- data[data$species == sp, ]  # Filter data for the current species
    
    # Create the plot using kernel density estimation
    plot_density_species <- ggplot(species_data, aes(x = soil_water_potential, y = .data[[ndvi_column]])) +
      stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +  # KDE heatmap
      scale_fill_viridis_c(option = "C", direction = -1) +  # Colorblind-friendly palette
      labs(
        title = paste("Kernel Density Plot of NDVI vs. Soil Water Potential -", sp),
        x = "Soil Water Potential",
        y = ndvi_label,
        fill = "Density"
      ) +
      theme_minimal(base_size = 14) +  # Base size for better readability
      theme(
        plot.background = element_rect(fill = "white", color = NA),    # Clean white background
        panel.background = element_rect(fill = "white", color = NA),  # Panel background
        legend.background = element_rect(fill = "white", color = NA), # Legend background
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold", color = "black"), # Centered bold title
        axis.title = element_text(face = "bold", size = 16),          # Bold axis titles
        axis.text = element_text(size = 12, color = "black"),         # Axis tick labels
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8), # Add subtle border
        legend.position = "top",                                      # Legend at the top
        legend.title = element_text(size = 14, face = "bold"),        # Bold legend title
        legend.text = element_text(size = 12)                        # Legend text size
      )
    
    # Save the plot with a unique filename
    species_file_name <- str_c(plot_density_folder, "/", sp, "_kernel_density_plot.png")
    ggsave(filename = species_file_name, plot = plot_density_species, width = 10, height = 8, dpi = 300)
    
    # Print a message indicating the plot was saved
    message("Saved kernel density plot for species: ", sp, " at ", species_file_name)
  }
  
  message("All plots saved in: ", plot_density_folder)
}

plot_points_NDVI_PSI <- function(data, save_path) {
  
  library(ggplot2)
  library(dplyr)
  
  # Determine the correct NDVI column
  ndvi_column <- if ("Quantiles" %in% names(data)) {
    "Quantiles"
  } else if ("Proportions" %in% names(data)) {
    "Proportions"
  } else {
    stop("Neither Quantiles nor Proportions column found in data.")
  }
  
  # Set NDVI y-axis label dynamically
  ndvi_label <- if (ndvi_column == "Quantiles") "NDVI (Quantiles)" else "NDVI (Proportions)"
  
  # Define the desired order of species
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  
  # Filter and arrange the data based on the defined species order
  data_filtered <- data %>% 
    filter(species %in% species_order) %>%
    mutate(species = factor(species, levels = species_order)) # Ensure order
  
  # Create the plot
  p <- ggplot(data_filtered, aes(x = soil_water_potential, y = .data[[ndvi_column]])) +
    geom_point(alpha = 0.5, color = "blue") +
    facet_wrap(~species, scales = "free") + # Ordered based on factor levels
    labs(x = "Soil Water Potential", y = ndvi_label, title = "Point Plot of NDVI vs. PSI by Species") +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),    # Clean white background
      panel.background = element_rect(fill = "white", color = NA),  # Panel background
      legend.background = element_rect(fill = "white", color = NA), # Legend background
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold", color = "black"), # Centered bold title
      axis.title = element_text(face = "bold", size = 16),          # Bold axis titles
      axis.text = element_text(size = 12, color = "black"),         # Axis tick labels
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8), # Add subtle border
      legend.position = "top",                                      # Legend at the top
      legend.title = element_text(size = 14, face = "bold"),        # Bold legend title
      legend.text = element_text(size = 12)                        # Legend text size
    )
  
  # Save the plot
  ggsave(save_path, plot = p, width = 10, height = 6, dpi = 300)
  
  # Return the plot
  return(p)
}

plot_boundaries_NDVI_PSI <- function(data, save_path) {
  
  library(ggplot2)
  library(dplyr)
  library(purrr)  # For map_dfr()
  
  # Determine the correct NDVI column
  ndvi_column <- if ("Quantiles" %in% names(data)) {
    "Quantiles"
  } else if ("Proportions" %in% names(data)) {
    "Proportions"
  } else {
    stop("Neither Quantiles nor Proportions column found in data.")
  }
  
  # Set NDVI y-axis label dynamically
  ndvi_label <- if (ndvi_column == "Quantiles") "NDVI (Quantiles)" else "NDVI (Proportions)"
  
  # Define the desired order of species
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  
  # Filter and arrange the data based on the defined species order
  data_filtered <- data %>% 
    filter(species %in% species_order) %>%
    mutate(species = factor(species, levels = species_order))  # Ensure correct order
  
  # Function to compute convex hull for each species
  compute_hull <- function(df) {
    df <- na.omit(df)  # Remove NA values
    if (nrow(df) < 3) return(NULL)  # Skip species with fewer than 3 points
    return(df[chull(df$soil_water_potential, df[[ndvi_column]]), ])
  }
  
  # Compute convex hulls, safely handling species with <3 points
  hull_data <- data_filtered %>%
    group_by(species) %>%
    group_split() %>%  # Split by species
    map_dfr(compute_hull)  # Apply hull computation and combine results
  
  # Create the plot
  p <- ggplot(data_filtered, aes(x = soil_water_potential, y = .data[[ndvi_column]])) +
    geom_point(alpha = 0.5, color = "blue") +  # Scatter plot
    geom_polygon(data = hull_data, aes(x = soil_water_potential, y = .data[[ndvi_column]], group = species), 
                 fill = NA, color = "red", linetype = "dashed") +  # Convex hull boundaries
    facet_wrap(~species, scales = "free") +  # Ordered by factor levels
    labs(x = "Soil Water Potential", y = ndvi_label, title = "Point Plot with Boundaries for NDVI PSI by Species") +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold", color = "black"),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(size = 12, color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      legend.position = "top",
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12)
    )
  
  # Save the plot
  ggsave(save_path, plot = p, width = 10, height = 6, dpi = 300)
  
  # Return the plot
  return(p)
}

calculate_no_data_area <- function(data, output_csv) {
  
  library(dplyr)
  
  # Determine the correct NDVI column
  ndvi_column <- if ("Quantiles" %in% names(data)) {
    "Quantiles"
  } else if ("Proportions" %in% names(data)) {
    "Proportions"
  } else {
    stop("Neither Quantiles nor Proportions column found in data.")
  }
  
  # Set NDVI label dynamically
  ndvi_label <- if (ndvi_column == "Quantiles") "NDVI (Quantiles)" else "NDVI (Proportions)"
  
  # Define the desired species order
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  
  # Filter for selected species
  data_filtered <- data %>% filter(species %in% species_order)
  
  # Initialize an empty dataframe to store results
  results <- data.frame(Species = character(), No_Data_Area = numeric(), stringsAsFactors = FALSE)
  
  # Loop through each species
  for (sp in species_order) {
    
    species_data <- data_filtered %>% filter(species == sp)  # Filter data for current species
    
    # Define bounding box (min/max X and Y)
    bbox <- species_data %>%
      summarise(
        min_x = min(soil_water_potential, na.rm = TRUE),
        max_x = max(soil_water_potential, na.rm = TRUE),
        min_y = min(.data[[ndvi_column]], na.rm = TRUE),
        max_y = max(.data[[ndvi_column]], na.rm = TRUE)
      )
    
    # Estimate full area of the bounding box
    full_area <- (bbox$max_x - bbox$min_x) * (bbox$max_y - bbox$min_y)
    
    # Count valid data points
    valid_points <- species_data %>% filter(!is.na(soil_water_potential) & !is.na(.data[[ndvi_column]]))
    
    # Estimate occupied area using convex hull
    if (nrow(valid_points) >= 3) {
      hull_points <- valid_points[chull(valid_points$soil_water_potential, valid_points[[ndvi_column]]), ]
      occupied_area <- abs(sum(
        hull_points$soil_water_potential[-1] * hull_points[[ndvi_column]][-length(hull_points[[ndvi_column]])] -
          hull_points$soil_water_potential[-length(hull_points$soil_water_potential)] * hull_points[[ndvi_column]][-1]
      )) / 2
    } else {
      occupied_area <- 0  # Not enough points to compute hull
    }
    
    # Calculate no-data area
    no_data_area <- full_area - occupied_area
    
    # Append results
    results <- rbind(results, data.frame(Species = sp, No_Data_Area = round(no_data_area, 2)))
  }
  
  # Save results to CSV
  write.csv(results, output_csv, row.names = FALSE)
  
  # Print a message
  message("No-Data Area calculation completed! Results saved to: ", output_csv)
}

calculate_no_data_area_reorder <- function(data, output_rdata) {
  
  library(dplyr)
  
  # Determine the correct NDVI column
  ndvi_column <- if ("Quantiles" %in% names(data)) {
    "Quantiles"
  } else if ("Proportions" %in% names(data)) {
    "Proportions"
  } else {
    stop("Neither Quantiles nor Proportions column found in data.")
  }
  
  # Set NDVI y-axis label dynamically
  ndvi_label <- if (ndvi_column == "Quantiles") "NDVI (Quantiles)" else "NDVI (Proportions)"
  
  # Define the species order
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  
  # Filter for selected species
  data_filtered <- data %>% filter(species %in% species_order)
  
  # Initialize dataframe to store 10,000 simulations
  results <- data.frame(non_data_area = numeric(), data_area = numeric())
  
  # Function to compute no-data area from shuffled data
  compute_no_data_area <- function(species_data) {
    
    # Randomly shuffle X and Y separately
    shuffled_x <- sample(species_data$soil_water_potential)
    shuffled_y <- sample(species_data[[ndvi_column]])
    
    # Create new shuffled dataset
    shuffled_data <- data.frame(soil_water_potential = shuffled_x, NDVI = shuffled_y)
    
    # Define bounding box (min/max X and Y)
    bbox <- shuffled_data %>%
      summarise(
        min_x = min(soil_water_potential, na.rm = TRUE),
        max_x = max(soil_water_potential, na.rm = TRUE),
        min_y = min(NDVI, na.rm = TRUE),
        max_y = max(NDVI, na.rm = TRUE)
      )
    
    # Compute full area of the bounding box
    full_area <- (bbox$max_x - bbox$min_x) * (bbox$max_y - bbox$min_y)
    
    # Compute occupied area using convex hull (if there are at least 3 points)
    if (nrow(shuffled_data) >= 3) {
      hull_points <- shuffled_data[chull(shuffled_data$soil_water_potential, shuffled_data$NDVI), ]
      occupied_area <- abs(sum(
        hull_points$soil_water_potential[-1] * hull_points$NDVI[-length(hull_points$NDVI)] -
          hull_points$soil_water_potential[-length(hull_points$soil_water_potential)] * hull_points$NDVI[-1]
      )) / 2
    } else {
      occupied_area <- 0  # Not enough points to compute hull
    }
    
    # Calculate no-data area
    no_data_area <- full_area - occupied_area
    
    return(c(no_data_area, occupied_area))
  }
  
  # Run the process 10,000 times
  for (i in 1:10000) {
    # Compute no-data area for shuffled dataset
    area_values <- compute_no_data_area(data_filtered)
    
    # Store the results in the dataframe
    results <- rbind(results, data.frame(non_data_area = round(area_values[1], 2),
                                         data_area = round(area_values[2], 2)))
  }
  
  # Save results as an RData file
  save(results, file = output_rdata)
  
  # Print completion message
  message("No-Data Area simulation completed! Results saved to: ", output_rdata)
}

plot_box_NDVI_PSI <- function(data, save_path = NULL) {
  library(dplyr)
  library(ggplot2)
  
  # Determine the correct NDVI column
  ndvi_column <- if ("Quantiles" %in% names(data)) {
    "Quantiles"
  } else if ("Proportions" %in% names(data)) {
    "Proportions"
  } else {
    stop("Neither Quantiles nor Proportions column found in data.")
  }
  
  # Set NDVI y-axis label dynamically
  ndvi_label <- if (ndvi_column == "Quantiles") "NDVI (Quantiles)" else "NDVI (Proportions)"
  
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
  
  
  plot <- ggplot(data, aes(x = factor(soil_bin), y = .data[[ndvi_column]], fill = species)) +
    geom_boxplot() +
    facet_wrap(~ species, scales = "free_x") +
    scale_fill_manual(values = cb_palette) +
    labs(x = "Soil Water Potential (Median of Bins)", y = ndvi_label) +
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

plot_box_merge_NDVI_PSI <- function(data, save_path = NULL) {
  library(dplyr)
  library(ggplot2)
  
  # Determine the correct NDVI column
  ndvi_column <- if ("Quantiles" %in% names(data)) {
    "Quantiles"
  } else if ("Proportions" %in% names(data)) {
    "Proportions"
  } else {
    stop("Neither Quantiles nor Proportions column found in data.")
  }
  
  # Set NDVI y-axis label dynamically
  ndvi_label <- if (ndvi_column == "Quantiles") "NDVI (Quantiles)" else "NDVI (Proportions)"
  
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
  
  plot <- ggplot(data, aes(x = factor(soil_bin), y = .data[[ndvi_column]], fill = species)) +
    geom_boxplot(position = position_dodge(width = 0.75), na.rm = TRUE) +
    scale_fill_manual(values = cb_palette) +
    labs(x = "Soil Water Potential (Median of Bins)", y = ndvi_label, fill = "Species") +
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

plot_mean_box_NDVI_PSI <- function(data, save_path = NULL) {
  library(dplyr)
  library(ggplot2)
  
  # Determine the correct NDVI column
  ndvi_column <- if ("Quantiles" %in% names(data)) {
    "Quantiles"
  } else if ("Proportions" %in% names(data)) {
    "Proportions"
  } else {
    stop("Neither Quantiles nor Proportions column found in data.")
  }
  
  # Set NDVI y-axis label dynamically
  ndvi_label <- if (ndvi_column == "Quantiles") "Mean NDVI (Quantiles)" else "Mean NDVI (Proportions)"
  
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
    summarise(mean_value = mean(.data[[ndvi_column]], na.rm = TRUE), .groups = 'drop')
  
  # 4. Plot the mean points and lines with color-blind friendly palette
  cb_palette <- c("Oak"   = "#E69F00",  # Orange
                  "Beech" = "#0072B2",  # A deep blue with lower luminance than "#56B4E9"
                  "Spruce"= "#009E73",  # Bluish-green
                  "Pine"  = "#F0E442") 
  
  
  plot <- ggplot(mean_data, aes(x = factor(soil_bin), y = mean_value, color = species, group = species)) +
    geom_point(size = 3) +
    geom_line() +
    scale_color_manual(values = cb_palette) +
    labs(x = "Soil Water Potential (Median of Bins)", y = ndvi_label, color = "Species") +
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

plot_mean_box_NDVI_PSI_with_slope <- function(data, save_path = NULL) {
  
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  
  data <- na.omit(data)
  
  # Determine the correct NDVI column dynamically
  ndvi_column <- if ("Quantiles" %in% names(data)) {
    "Quantiles"
  } else if ("Proportions" %in% names(data)) {
    "Proportions"
  } else {
    stop("Neither Quantiles nor Proportions column found in data.")
  }
  
  # Set NDVI y-axis label dynamically
  ndvi_label <- if (ndvi_column == "Quantiles") "Mean NDVI (Quantiles)" else "Mean NDVI (Proportions)"
  
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
  
  # Compute mean NDVI for each species within each soil bin (For Panel A)
  mean_data <- data %>%
    group_by(soil_bin, species) %>%
    summarise(mean_value = mean(.data[[ndvi_column]], na.rm = TRUE), .groups = 'drop')
  
  # Function to calculate slope, R², p-value, and standard error using original soil water potential (For Panel B)
  calculate_slope_stats <- function(data, ndvi_column) {
    # Set filtering range dynamically
    if (ndvi_column == "Quantiles") {
      data_filtered <- data %>% filter(.data[[ndvi_column]] >= 7 & .data[[ndvi_column]] <= 14)
    } else { # Proportions case
      data_filtered <- data %>% filter(.data[[ndvi_column]] >= 0.9 & .data[[ndvi_column]] <= 1.3)
    }
    
    regression_results <- data_filtered %>%
      group_by(species) %>%
      filter(n() >= 2) %>%
      summarise(
        model = list(lm(.data[[ndvi_column]] ~ soil_water_potential, data = cur_data())),  # Use original soil_water_potential
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
  
  slope_results <- calculate_slope_stats(data, ndvi_column)
  
  # Color Palette
  cb_palette <- c("Oak"   = "#E69F00",  # Orange
                  "Beech" = "#0072B2",  # A deep blue with lower luminance than "#56B4E9"
                  "Spruce"= "#009E73",  # Bluish-green
                  "Pine"  = "#F0E442") 
  
  
  # Scatter plot (Panel A) with binned soil water potential
  plot_a <- ggplot(mean_data, aes(x = factor(soil_bin), y = mean_value, color = species, group = species)) +
    geom_point(size = 3) +
    geom_line() +
    scale_color_manual(values = cb_palette) +
    labs(x = "Soil Water Potential (Median of 100 KPa Bins)", y = ndvi_label, color = "Species") +
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
  
  # Bar plot (Panel B) with slopes, p-values, R², and error bars using original soil water potential
  plot_b <- ggplot(slope_results, aes(y = reorder(species, slope), x = slope, fill = species)) +
    geom_bar(stat = "identity", width = 0.6, color = "black") +
    geom_errorbar(aes(xmin = slope - std_error, xmax = slope + std_error), width = 0.2) +  # Error bars for slope uncertainty
    scale_fill_manual(values = cb_palette) +
    labs(x = paste("Slope (", ndvi_column, " range applied)")) +  # Removed y-axis label
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
  
  # Combine plots with adjusted widths
  combined_plot <- plot_a + plot_b + plot_annotation(tag_levels = 'A') + plot_layout(widths = c(2, 1))  # Adjusted layout
  
  # Save if needed
  if (!is.null(save_path)) {
    ggsave(filename = save_path, plot = combined_plot, width = 12, height = 6, dpi = 300)
  }
  
  return(combined_plot)
}

plot_mean_box_NDVI_PSI_with_slope_bin <- function(data, save_path = NULL) {
  
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  
  data <- na.omit(data)
  
  # Determine the correct NDVI column dynamically
  ndvi_column <- if ("Quantiles" %in% names(data)) {
    "Quantiles"
  } else if ("Proportions" %in% names(data)) {
    "Proportions"
  } else {
    stop("Neither Quantiles nor Proportions column found in data.")
  }
  
  # Set NDVI y-axis label dynamically
  ndvi_label <- if (ndvi_column == "Quantiles") "Mean NDVI (Quantiles)" else "Mean NDVI (Proportions)"
  
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
  
  # Compute mean NDVI for each species within each soil bin (For Panel A)
  mean_data <- data %>%
    group_by(soil_bin, species) %>%
    summarise(mean_value = mean(.data[[ndvi_column]], na.rm = TRUE), .groups = 'drop')
  
  # Function to calculate slope, R², p-value, and standard error using mean NDVI from binned data (For Panel B)
  calculate_slope_stats_bin <- function(mean_data, ndvi_column) {
    # Set filtering range dynamically
    if (ndvi_column == "Quantiles") {
      data_filtered <- mean_data %>% filter(mean_value >= 7 & mean_value <= 14)
    } else { # Proportions case
      data_filtered <- mean_data %>% filter(mean_value >= 0.9 & mean_value <= 1.3)
    }
    
    regression_results <- data_filtered %>%
      group_by(species) %>%
      filter(n() >= 2) %>%
      summarise(
        model = list(lm(mean_value ~ soil_bin, data = cur_data())),  # Regression on binned data
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
  
  slope_results <- calculate_slope_stats_bin(mean_data, ndvi_column)
  
  # Color Palette
  cb_palette <- c("Oak"   = "#E69F00",  # Orange
                  "Beech" = "#0072B2",  # A deep blue with lower luminance than "#56B4E9"
                  "Spruce"= "#009E73",  # Bluish-green
                  "Pine"  = "#F0E442") 
  
  
  # Scatter plot (Panel A) with binned soil water potential
  plot_a <- ggplot(mean_data, aes(x = factor(soil_bin), y = mean_value, color = species, group = species)) +
    geom_point(size = 3) +
    geom_line() +
    scale_color_manual(values = cb_palette) +
    labs(x = "Soil Water Potential (Median of 100 KPa Bins)", y = ndvi_label, color = "Species") +
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
    labs(x = paste("Slope (Binned", ndvi_column, "Range Applied)")) +  # Dynamic label
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
  
  # Combine plots with adjusted widths
  combined_plot <- plot_a + plot_b + plot_annotation(tag_levels = 'A') + plot_layout(widths = c(1.5, 1))  # Adjusted layout
  
  # Save if needed
  if (!is.null(save_path)) {
    ggsave(filename = save_path, plot = combined_plot, width = 12, height = 6, dpi = 300)
  }
  
  return(combined_plot)
}

#--------------------------------------------------------------------------
# Function: Time Series Plot (Vertical Panels)
# Plots time series for NDVI (quantiles/proportions), soil water potential, and transpiration deficit.
#--------------------------------------------------------------------------
plot_time_series_NDVI_PSI_TDiff_avg <- function(df_all, plot_ts_path) {
  
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  
  # Define the custom theme
  custom_theme <- theme_minimal() +
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
      legend.position = "top",
      plot.caption = element_text(face = "bold", size = 16, hjust = 0),
      plot.caption.position = "plot"
    )
  
  # Get averaged data (assumed to be used for both NDVI and PSI/TDiff)
  df_avg <- df_average_year(df_all)
  
  # Determine the correct NDVI column from df_avg
  ndvi_column <- if ("avg_quantile" %in% names(df_avg)) {
    "avg_quantile"
  } else if ("avg_proportion" %in% names(df_avg)) {
    "avg_proportion"
  } else {
    stop("Neither avg_quantile nor avg_proportion column found in df_avg.")
  }
  
  # Set NDVI label based on the data type
  ndvi_label <- if (ndvi_column == "avg_quantile") {
    "Average NDVI (Quantiles)"
  } else {
    "Average NDVI (Proportions)"
  }
  
  # Panel (a): NDVI Time Series over Years
  p1 <- ggplot(df_avg, aes(x = as.numeric(year), y = .data[[ndvi_column]])) +
    geom_point(color = "orange", size = 3, alpha = 0.7) +
    geom_line(color = "orange") +
    labs(title = "NDVI over Years",
         x = "Year", 
         y = ndvi_label) +
    custom_theme
  
  # Panel (b): Soil Water Potential (PSI) Time Series over Years
  p2 <- ggplot(df_avg, aes(x = as.numeric(year), y = avg_psi)) +
    geom_point(color = "purple", size = 3, alpha = 0.7) +
    geom_line(color = "purple") +
    labs(title = "Soil Water Potential over Years",
         x = "Year", 
         y = "Average PSI") +
    custom_theme
  
  # Panel (c): Transpiration Deficit (TDiff) Time Series over Years
  p3 <- ggplot(df_avg, aes(x = as.numeric(year), y = avg_tdiff)) +
    geom_point(color = "green4", size = 3, alpha = 0.7) +
    geom_line(color = "green4") +
    labs(title = "Transpiration Deficit over Years",
         x = "Year", 
         y = "Average TDiff") +
    custom_theme
  
  # Combine the three plots vertically (one row per panel)
  combined_ts <- p1 / p2 / p3
  
  # Ensure the output directory exists
  out_dir <- dirname(plot_ts_path)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  # Save the time series plot
  ggsave(filename = plot_ts_path, plot = combined_ts, width = 10, height = 15, dpi = 300, device = "png")
}

#--------------------------------------------------------------------------
# Functions: Correlation Plot (Scatter Plots with Regression Lines)
# Plots correlations between NDVI and PSI, and NDVI and TDiff.
#--------------------------------------------------------------------------
plot_correlation_NDVI_PSI_TDiff_corr <- function(df_all, plot_corr_path) {
  
  library(ggplot2)
  library(patchwork)
  
  # For correlation with PSI, use averaged data from df_average_year()
  df_avg <- df_average_year(df_all)
  ndvi_column <- if ("avg_quantile" %in% names(df_avg)) {
    "avg_quantile"
  } else if ("avg_proportion" %in% names(df_avg)) {
    "avg_proportion"
  } else {
    stop("Neither avg_quantile nor avg_proportion column found in df_avg.")
  }
  
  ndvi_label <- if (ndvi_column == "avg_quantile") {
    "Average NDVI (Quantiles)"
  } else {
    "Average NDVI (Proportions)"
  }
  
  # Calculate correlation coefficient between NDVI and Soil Water Potential (PSI)
  correlation_psi <- cor(df_avg$avg_psi, df_avg[[ndvi_column]], use = "complete.obs")
  
  plot_psi <- ggplot(df_avg, aes(x = avg_psi, y = .data[[ndvi_column]])) +
    geom_point(color = "purple", size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", color = "black", se = TRUE) +
    annotate("text", 
             x = min(df_avg$avg_psi, na.rm = TRUE), 
             y = max(df_avg[[ndvi_column]], na.rm = TRUE), 
             label = paste("Correlation:", round(correlation_psi, 2)),
             hjust = 0, vjust = 1.5, size = 5, color = "blue") +
    labs(title = "NDVI vs Soil Water Potential",
         x = "Average PSI", 
         y = ndvi_label) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  
  # For correlation with TDiff, use averaged data from df_average_year_NDVI_TDiff()
  df_avg_tdiff <- df_average_year_NDVI_TDiff(df_all)
  
  correlation_tdiff <- cor(df_avg_tdiff$avg_tdiff, df_avg_tdiff[[ndvi_column]], use = "complete.obs")
  
  plot_tdiff <- ggplot(df_avg_tdiff, aes(x = avg_tdiff, y = .data[[ndvi_column]])) +
    geom_point(color = "green4", size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", color = "black", se = TRUE) +
    annotate("text", 
             x = min(df_avg_tdiff$avg_tdiff, na.rm = TRUE), 
             y = max(df_avg_tdiff[[ndvi_column]], na.rm = TRUE), 
             label = paste("Correlation:", round(correlation_tdiff, 2)),
             hjust = 0, vjust = 1.5, size = 5, color = "blue") +
    labs(title = "NDVI vs Transpiration Deficit",
         x = "Average TDiff", 
         y = ndvi_label) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  
  # Combine the two correlation plots side-by-side
  combined_corr <- plot_psi + plot_tdiff + plot_layout(ncol = 2, guides = "collect")
  
  # Ensure the output directory exists
  out_dir <- dirname(plot_corr_path)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  # Save the correlation plot figure
  ggsave(filename = plot_corr_path, plot = combined_corr, width = 15, height = 7, dpi = 300)
}

plot_time_series_NDVI_PSI_TDiff_species_avg <- function(df_all, output_path) {
  library(ggplot2)
  library(patchwork)
  
  # Get species-level averaged data
  df_species <- df_average_year_species(df_all)
  
  # Set species order and assign factor levels
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  df_species$species <- factor(df_species$species, levels = species_order)
  
  # Define the color palette for the species
  cb_palette <- c("Oak"   = "#E69F00", 
                  "Beech" = "#0072B2",
                  "Spruce"= "#009E73", 
                  "Pine"  = "#F0E442")
  
  # Determine the NDVI column and label
  ndvi_column <- if ("avg_quantile" %in% names(df_species)) {
    "avg_quantile"
  } else if ("avg_proportion" %in% names(df_species)) {
    "avg_proportion"
  } else {
    stop("Neither avg_quantile nor avg_proportion column found in species data.")
  }
  
  ndvi_label <- if (ndvi_column == "avg_quantile") {
    "Average NDVI (Quantiles)"
  } else {
    "Average NDVI (Proportions)"
  }
  
  # Define the custom theme
  custom_theme <- theme_minimal() +
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
      legend.position = "top",
      plot.caption = element_text(face = "bold", size = 16, hjust = 0),
      plot.caption.position = "plot"
    )
  
  # Panel (a): NDVI Time Series over Years (per species)
  p1 <- ggplot(df_species, aes(x = as.numeric(year), y = .data[[ndvi_column]], color = species)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_line() +
    labs(title = "NDVI over Years",
         x = "Year", 
         y = ndvi_label,
         color = "Species") +
    scale_color_manual(values = cb_palette) +
    custom_theme
  
  # Panel (b): Soil Water Potential (PSI) Time Series over Years (per species)
  p2 <- ggplot(df_species, aes(x = as.numeric(year), y = avg_psi, color = species)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_line() +
    labs(title = "Soil Water Potential over Years",
         x = "Year", 
         y = "Average PSI",
         color = "Species") +
    scale_color_manual(values = cb_palette) +
    custom_theme
  
  # Panel (c): Transpiration Deficit (TDiff) Time Series over Years (per species)
  p3 <- ggplot(df_species, aes(x = as.numeric(year), y = avg_tdiff, color = species)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_line() +
    labs(title = "Transpiration Deficit over Years",
         x = "Year", 
         y = "Average TDiff",
         color = "Species") +
    scale_color_manual(values = cb_palette) +
    custom_theme
  
  # Combine the three panels vertically into one figure
  combined_ts <- p1 / p2 / p3
  
  # Save the combined figure to the specified output path
  ggsave(filename = output_path, plot = combined_ts, width = 10, height = 15, dpi = 300)
}

plot_correlation_NDVI_PSI_TDiff_species_corr <- function(df_all, out_dir) {
  library(ggplot2)
  library(patchwork)
  
  # Use species-level averaged data (assumes this data frame includes NDVI, PSI, TDiff, and species columns)
  df_all_avg <- df_average_year_species(df_all)
  
  unique_species <- unique(df_all_avg$species)
  
  # Ensure the output directory exists
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  for (sp in unique_species) {
    # Subset data for the current species
    df_sp <- subset(df_all_avg, species == sp)
    
    # Determine NDVI column and label
    ndvi_column <- if ("avg_quantile" %in% names(df_sp)) {
      "avg_quantile"
    } else if ("avg_proportion" %in% names(df_sp)) {
      "avg_proportion"
    } else {
      stop("Neither avg_quantile nor avg_proportion column found in species data.")
    }
    
    ndvi_label <- if (ndvi_column == "avg_quantile") {
      "Average NDVI (Quantiles)"
    } else {
      "Average NDVI (Proportions)"
    }
    
    # Correlation plot for NDVI vs Soil Water Potential (PSI)
    correlation_psi <- cor(df_sp$avg_psi, df_sp[[ndvi_column]], use = "complete.obs")
    plot_psi <- ggplot(df_sp, aes(x = avg_psi, y = .data[[ndvi_column]])) +
      geom_point(color = "purple", size = 3, alpha = 0.7) +
      geom_smooth(method = "lm", color = "black", se = TRUE) +
      annotate("text", 
               x = min(df_sp$avg_psi, na.rm = TRUE), 
               y = max(df_sp[[ndvi_column]], na.rm = TRUE), 
               label = paste("Correlation:", round(correlation_psi, 2)),
               hjust = 0, vjust = 1.5, size = 5, color = "blue") +
      labs(title = paste("NDVI vs Soil Water Potential for", sp),
           x = "Average PSI", y = ndvi_label) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    
    # Correlation plot for NDVI vs Transpiration Deficit (TDiff)
    correlation_tdiff <- cor(df_sp$avg_tdiff, df_sp[[ndvi_column]], use = "complete.obs")
    plot_tdiff <- ggplot(df_sp, aes(x = avg_tdiff, y = .data[[ndvi_column]])) +
      geom_point(color = "green4", size = 3, alpha = 0.7) +
      geom_smooth(method = "lm", color = "black", se = TRUE) +
      annotate("text", 
               x = min(df_sp$avg_tdiff, na.rm = TRUE), 
               y = max(df_sp[[ndvi_column]], na.rm = TRUE), 
               label = paste("Correlation:", round(correlation_tdiff, 2)),
               hjust = 0, vjust = 1.5, size = 5, color = "blue") +
      labs(title = paste("NDVI vs Transpiration Deficit for", sp),
           x = "Average TDiff", y = ndvi_label) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    
    # Combine the two correlation plots side-by-side
    combined_corr <- plot_psi + plot_tdiff + plot_layout(ncol = 2, guides = "collect")
    
    # Save the correlation plot for the current species
    file_name <- file.path(out_dir, paste0("correlation_species_", sp, ".png"))
    ggsave(filename = file_name, plot = combined_corr, width = 15, height = 7, dpi = 300)
  }
}

plot_time_series_and_correlation_combined <- function(df_all, output_path) {
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(cowplot)  # For get_legend and plot_grid
  
  # Get species-level averaged data
  df_species <- df_average_year_species(df_all)
  
  # Set species order and assign factor levels
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  df_species$species <- factor(df_species$species, levels = species_order)
  
  # Define the color palette for the species
  cb_palette <- c("Oak"   = "#E69F00", 
                  "Beech" = "#0072B2",
                  "Spruce"= "#009E73", 
                  "Pine"  = "#F0E442")
  
  # Determine NDVI column and label
  ndvi_column <- if ("avg_quantile" %in% names(df_species)) {
    "avg_quantile"
  } else if ("avg_proportion" %in% names(df_species)) {
    "avg_proportion"
  } else {
    stop("Neither avg_quantile nor avg_proportion column found in species data.")
  }
  
  ndvi_label <- if (ndvi_column == "avg_quantile") {
    "NDVI quantiles"
  } else {
    "NDVI proportions"
  }
  
  # Define a custom theme with consistent text sizes
  custom_theme <- theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      legend.text = element_text(color = "black", size = 14),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      plot.caption = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.caption.position = "plot",
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(color = "black", size = 14),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  # Create an x-axis scale that shows breaks every 2 years from 2003 to 2024
  x_scale <- scale_x_continuous(breaks = seq(2003, 2024, by = 2), limits = c(2003, 2024))
  
  ### TIME SERIES PLOTS ###
  # Panel (a): NDVI Time Series – keep the legend (with horizontal layout) for extraction later
  p1 <- ggplot(df_species, aes(x = as.numeric(year), y = .data[[ndvi_column]], color = species)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_line() +
    labs(title = "", x = "", y = expression(atop(bold("NDVI quantiles"), bold("(rank)"))), color = "Species") +
    scale_color_manual(values = cb_palette, name = "", guide = guide_legend(nrow = 1)) +
    custom_theme +
    labs(caption = "(a)") +
    theme(legend.position = "bottom", legend.direction = "horizontal") +
    x_scale
  
  # Panel (b): soil water potential Time Series – remove legend
  p2 <- ggplot(df_species, aes(x = as.numeric(year), y = avg_psi, color = species)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_line() +
    labs(title = "", x = "", y = expression(atop(bold("soil water potential"), bold("(kPa)"))), color = "Species") +
    scale_color_manual(values = cb_palette, name = "") +
    custom_theme +
    guides(color = "none") +
    labs(caption = "(b)") +
    x_scale
  
  # Panel (c): transpiration deficit Time Series – remove legend
  p3 <- ggplot(df_species, aes(x = as.numeric(year), y = avg_tdiff, color = species)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_line() +
    labs(title = "", x = "", y = expression(atop(bold("transpiration deficit"), bold("(mm)"))), color = "Species") +
    scale_color_manual(values = cb_palette, name = "") +
    custom_theme +
    guides(color = "none") +
    labs(caption = "(c)") +
    x_scale
  
  # Combine time series panels vertically
  ts_plot <- p1 / p2 / p3
  
  ### CORRELATION PLOTS ###
  # Compute correlations for panel (d): NDVI vs soil water potential
  annotations_d <- df_species %>%
    group_by(species) %>%
    summarize(corr = round(cor(avg_psi, .data[[ndvi_column]], use = "complete.obs"), 2))
  
  # Compute correlations for panel (e): NDVI vs transpiration deficit
  annotations_e <- df_species %>%
    group_by(species) %>%
    summarize(corr = round(cor(avg_tdiff, .data[[ndvi_column]], use = "complete.obs"), 2))
  
  # Create annotation text (ordered by species_order)
  text_d <- paste0(species_order, ": r = ", 
                   annotations_d$corr[match(species_order, annotations_d$species)],
                   collapse = "\n")
  
  text_e <- paste0(species_order, ": r = ", 
                   annotations_e$corr[match(species_order, annotations_e$species)],
                   collapse = "\n")
  
  # Panel (d): NDVI vs soil water potential – remove legend
  p_d <- ggplot(df_species, aes(x = avg_psi, y = .data[[ndvi_column]], color = species)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(x = "soil water potential (kPa)", y = expression(atop(bold("NDVI quantiles"), bold("(rank)"))), caption = "(d)", color = "Species") +
    scale_color_manual(values = cb_palette, name = "") +
    custom_theme +
    guides(color = "none") +
    annotate("text", x = -Inf, y = Inf, label = text_d, hjust = 0, vjust = 1, size = 5, color = "black")
  
  # Panel (e): NDVI vs transpiration deficit – remove legend
  p_e <- ggplot(df_species, aes(x = avg_tdiff, y = .data[[ndvi_column]], color = species)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(x = "transpiration deficit (mm)", y = "", caption = "(e)", color = "Species") +
    scale_color_manual(values = cb_palette, name = "") +
    custom_theme +
    guides(color = "none") +
    annotate("text", x = Inf, y = Inf, label = text_e, hjust = 1, vjust = 1, size = 5, color = "black")
  
  # Combine correlation panels side-by-side
  corr_plot <- p_d | p_e
  
  # Combine time series on top and correlations on bottom, and remove legends from this combination
  combined_main <- ts_plot / corr_plot
  combined_main <- combined_main & theme(legend.position = "top")
  
  # Extract the legend from p1 (now a single horizontal guide)
  legend <- suppressWarnings(cowplot::get_legend(p1))
  
  # Combine the main plot and the legend with the legend underneath
  final_plot <- cowplot::plot_grid(combined_main, legend, ncol = 1, rel_heights = c(1, 0.1))
  
  print(final_plot)
  ggsave(filename = output_path, plot = final_plot, width = 12, height = 14, dpi = 300)
}
