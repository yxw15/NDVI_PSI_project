library(ncdf4)
library(reshape2)
library(ggplot2)
library(terra)

### Transfer Netcdf from DWD to Dataframe ### 
transfer_psi_to_df <- function(nc_file, start_date) {
  
  # Step 1: Open the NetCDF file
  nc <- nc_open(nc_file)
  
  # Step 2: Extract variables
  time <- ncvar_get(nc, "time")
  depth <- ncvar_get(nc, "depth")
  x <- ncvar_get(nc, "x")
  y <- ncvar_get(nc, "y")
  psi <- ncvar_get(nc, "psi")
  
  # Convert time from "days since" format to actual dates
  dates <- as.Date(start_date) + time
  
  # Close the NetCDF file
  nc_close(nc)
  
  # Step 3: Assign meaningful dimension names
  dimnames(psi) <- list(
    x = x,
    y = y,
    depth = depth,
    time = as.character(dates)
  )
  
  # Step 4: Transform the 4D array into a data frame
  psi_melted <- melt(psi, varnames = c("x", "y", "depth", "time"), value.name = "soil_water_potential")
  
  psi_melted <- na.omit(psi_melted)
  return(psi_melted)
}

save_psi_melted <- function(psi_melted, file_path) {
  write.csv(psi_melted, file_path, row.names = FALSE)
}

filter_psi <- function(psi_melted, month_day, depth) {
  psi_melted <- na.omit(psi_melted)
  psi_filter <- subset(psi_melted, 
                       format(as.Date(time), "%m-%d") == month_day & depth == depth)
  return(psi_filter)
}

save_psi_raster <- function(psi_melted, month_day, depth, output_dir) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  psi <- filter_psi(psi_melted, month_day, depth)
  psi$time <- as.Date(psi$time, format = "%Y-%m-%d")
  years <- 2003:2024
  for (year in years) {
    psi_year <- psi %>% filter(format(time, "%Y") == as.character(year))
    psi_selected <- psi_year[, c("x", "y", "soil_water_potential")]
    psi_raster <- rast(psi_selected, type = "xyz")
    crs(psi_raster) <- "epsg:31467"
    output_path <- file.path(output_dir, paste0("psi_", year, ".tif"))
    writeRaster(psi_raster, output_path, overwrite = TRUE)
    cat("Raster saved for year:", year, "at", output_path, "\n")
  }
}

create_plot <- function(filtered_data, output_path) {
  filtered_data$time <- as.Date(filtered_data$time)
  filtered_data$year <- format(filtered_data$time, "%Y")
  
  plot <- ggplot(filtered_data, aes(x = x, y = y, fill = soil_water_potential)) +
    geom_tile() +
    scale_fill_viridis_c() +
    labs(
      title = "Soil Water Potential on July 28th for Each Year (2003-2024)",
      x = "X Coordinate",
      y = "Y Coordinate",
      fill = "Soil Water Potential (kPa)"
    ) +
    theme_minimal() +
    facet_wrap(~ year, ncol = 4)
  
  ggsave(output_path, plot, width = 20, height = 15, dpi = 300)
  plot
}

plot_time_series <- function(psi_df, pixels) {
  psi_df <- na.omit(psi_df)
  
  # Step 1: Select 5000 unique pixels
  set.seed(123) # For reproducibility
  unique_pixels <- unique(psi_df[, c("x", "y")]) # Extract unique (x, y) combinations
  selected_pixels <- unique_pixels[sample(nrow(unique_pixels), pixels), ]
  
  # Step 2: Filter data for selected pixels
  selected_data <- merge(psi_df, selected_pixels, by = c("x", "y"))
  
  # Step 3: Aggregate mean soil water potential for each time
  mean_time_series <- aggregate(
    soil_water_potential ~ time, 
    data = selected_data, 
    FUN = mean
  )
  
  # Step 4: Plot time-series of mean soil water potential
  ggplot(mean_time_series, aes(x = as.Date(time), y = soil_water_potential)) +
    geom_line() +
    labs(
      title = "Time-Series of Mean Soil Water Potential",
      x = "Date",
      y = "Mean Soil Water Potential (kPa)"
    ) +
    theme_minimal()
}

NDVI_PSIbin <- function(df, bin_width = 50) {
  
  library(dplyr)
  
  # Identify correct value column
  value_column <- if ("Quantiles" %in% names(df)) "Quantiles" else "Proportions"
  
  # Total pixel count per species
  species_totals <- df %>%
    group_by(species) %>%
    summarise(total_pixels = n(), .groups = "drop")
  
  # Define bin breaks
  psi_min <- floor(min(df$soil_water_potential, na.rm = TRUE))
  psi_max <- ceiling(max(df$soil_water_potential, na.rm = TRUE))
  bin_breaks <- seq(psi_min, psi_max, by = bin_width)
  
  # Bin the soil water potential values
  df <- df %>%
    mutate(PSI_bin = cut(soil_water_potential, breaks = bin_breaks, include.lowest = TRUE, right = FALSE))
  
  # Compute statistics for each species and bin, including filtering out bins with count < 2000
  meanNDVI_PSIbin_species <- df %>%
    group_by(species, PSI_bin) %>%
    summarise(
      avg_value = mean(.data[[value_column]], na.rm = TRUE),
      count = n(),
      .groups = 'drop'
    ) %>%
    mutate(
      bin_median = sapply(as.character(PSI_bin), function(bin_label) {
        nums <- as.numeric(strsplit(gsub("\\[|\\]|\\(|\\)", "", bin_label), ",")[[1]])
        mean(nums)
      })
    ) %>%
    left_join(species_totals, by = "species") %>%
    mutate(percentage = count / total_pixels) %>%
    filter(percentage >= 0.001) %>%
    select(species, PSI_bin, bin_median, avg_value, count, total_pixels, percentage)
  
  return(meanNDVI_PSIbin_species)
}

NDVI_PSIbin_log <- function(df, bin_width_log = 0.5) {
  library(dplyr)
  
  # Identify correct value column
  value_column <- if ("Quantiles" %in% names(df)) "Quantiles" else "Proportions"
  
  # Filter to valid, non-zero potentials
  df <- df %>% filter(!is.na(soil_water_potential), soil_water_potential != 0)
  
  # Log-transform magnitude of soil water potential
  df <- df %>% mutate(log_PSI = log(abs(soil_water_potential)))
  
  # Total pixel count per species
  species_totals <- df %>%
    group_by(species) %>%
    summarise(total_pixels = n(), .groups = "drop")
  
  # Define log-scale bin breaks
  psi_min <- floor(min(df$log_PSI, na.rm = TRUE))
  psi_max <- ceiling(max(df$log_PSI, na.rm = TRUE))
  bin_breaks <- seq(psi_min, psi_max, by = bin_width_log)
  
  # Bin log(soil_water_potential magnitude)
  df <- df %>%
    mutate(PSI_bin = cut(log_PSI, breaks = bin_breaks, include.lowest = TRUE, right = FALSE))
  
  # Compute stats per species and bin
  result <- df %>%
    group_by(species, PSI_bin) %>%
    summarise(
      avg_value = mean(.data[[value_column]], na.rm = TRUE),
      count = n(),
      .groups = 'drop'
    ) %>%
    mutate(
      bin_median = sapply(as.character(PSI_bin), function(bin_label) {
        edges <- as.numeric(strsplit(gsub("\\[|\\]|\\(|\\)", "", bin_label), ",")[[1]])
        mean(edges)
      })
    ) %>%
    left_join(species_totals, by = "species") %>%
    mutate(percentage = count / total_pixels) %>%
    filter(percentage >= 0.01) %>%
    select(species, PSI_bin, bin_median, avg_value, count, total_pixels, percentage)
  
  return(result)
}

TDiff_PSIbin <- function(df, bin_width = 50) {
  
  library(dplyr)
  
  # Identify the correct column
  value_column <- if ("transpiration_deficit" %in% names(df)) "transpiration_deficit"
  
  # Total pixel count per species
  species_totals <- df %>%
    group_by(species) %>%
    summarise(total_pixels = n(), .groups = "drop")
  
  # Define bin breaks dynamically
  psi_min <- floor(min(df$soil_water_potential, na.rm = TRUE))
  psi_max <- ceiling(max(df$soil_water_potential, na.rm = TRUE))
  bin_breaks <- seq(psi_min, psi_max, by = bin_width)
  
  df <- df %>%
    mutate(PSI_bin = cut(soil_water_potential, breaks = bin_breaks, include.lowest = TRUE, right = FALSE))
  
  # Compute mean transpiration_deficit and count per bin per species
  meanTDiff_PSIbin_species <- df %>%
    group_by(species, PSI_bin) %>%
    summarise(
      avg_transpiration_deficit = mean(.data[[value_column]], na.rm = TRUE),
      count = n(),
      .groups = 'drop'
    ) %>%
    mutate(
      bin_median = sapply(as.character(PSI_bin), function(bin_label) {
        nums <- as.numeric(strsplit(gsub("\\[|\\]|\\(|\\)", "", bin_label), ",")[[1]])
        mean(nums)
      })
    ) %>%
    left_join(species_totals, by = "species") %>%
    mutate(percentage = count / total_pixels) %>%
    filter(percentage > 0.001) %>%
    select(species, PSI_bin, bin_median, avg_transpiration_deficit, count, total_pixels, percentage)
  
  return(meanTDiff_PSIbin_species)
}

plot_NDVI_Q_PSIbin_exp_slope <- function(data, save_coeff_fig, save_slope_fig) {
  
  # Process data with your custom NDVI_PSIbin function and remove missing values
  data <- NDVI_PSIbin(data)
  data <- na.omit(data)
  
  # Load required libraries
  library(ggplot2)
  library(nlme)
  library(dplyr)
  library(tibble)
  library(patchwork)
  library(purrr)
  library(car)      # for deltaMethod
  library(tidyr)    # for pivot_longer
  
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
  
  # Set threshold manually to 11.5 instead of computing the median from the data
  threshold <- 11.5
  
  #### NONLINEAR MODELING PER SPECIES ####
  # Starting parameters and control parameters for the nls models
  start_list <- list(a = 5, b = 3, c = 0.001)
  control_params <- nls.control(maxiter = 200, minFactor = 1e-4)
  
  nls_models <- nlsList(avg_value ~ a + b * exp(-c * x) | species,
                        data = data_clean,
                        start = start_list,
                        control = control_params)
  print(summary(nls_models))
  
  # Extract coefficients for each species into a data frame and ensure species order
  coef_df <- as.data.frame(coef(nls_models), optional = TRUE) %>% 
    rownames_to_column(var = "species") %>%
    filter(!is.na(a))
  coef_df$species <- factor(coef_df$species, levels = species_order)
  
  # Calculate x50 and the slope at x50 (for Panel B)
  coef_df <- coef_df %>%
    mutate(x50 = ifelse((threshold - a) > 0 & b > 0,
                        -log((threshold - a)/b) / c,
                        NA),
           slope50 = -c * (threshold - a))
  
  #### PANEL A: Observed Data and Fitted Curves ####
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
  
  # Reverse the x-axis so soil water potential is positive
  x_scale <- scale_x_continuous(trans = "reverse", labels = function(x) -x)
  
  p_combined <- ggplot() +
    geom_point(data = data_clean, aes(x = x, y = avg_value, color = species)) +
    geom_line(data = pred_all, aes(x = x, y = pred, color = species), linewidth = 1) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "black", linewidth = 1) +
    annotate("text", x = 2000, y = threshold, label = "median", 
             hjust = -0.1, vjust = -0.3, fontface = "italic", size = 4) +
    scale_color_manual(values = cb_palette) +
    x_scale +
    labs(x = "Soil Water Potential", y = "NDVI Quantiles") +
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
    labs(x = "", y = "Soil Water Potential") +
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
  
  #### PANEL C: Bar Plot of Absolute Slope at x50 with Error Bars and Annotations ####
  stats_list <- lapply(levels(data_clean$species), function(sp) {
    mod <- nls_models[[sp]]
    coefs <- coef(mod)
    slope50_est <- -coefs["c"] * (threshold - coefs["a"])
    
    # Compute standard error using deltaMethod for c*(a - threshold)
    dm_result <- deltaMethod(mod, paste0("c*(a - ", threshold, ")"),
                             parameterNames = c("a", "b", "c"))
    se <- dm_result$SE
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
  stats_df$species <- factor(stats_df$species, levels = species_order)
  
  # Create annotation labels for Panel C:
  # If p < 0.05, display R² with an appended star; otherwise, display R² and p-value on a new line.
  stats_df <- stats_df %>%
    mutate(label_text = ifelse(p_val < 0.05,
                               sprintf("%.2f*", r_squared),
                               sprintf("%.2f\np = %.2f", r_squared, p_val)))
  
  p_slope <- ggplot(stats_df, aes(x = species, y = slope_abs, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, slope_abs - se), ymax = slope_abs + se), width = 0.2) +
    geom_text(aes(y = slope_abs/2, label = label_text), color = "black", size = 4) +
    labs(x = "", y = "Absolute Slope") +
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
  
  #### COMBINE PANELS A, B, and C INTO THE SLOPE FIGURE ####
  final_slope_plot <- p_combined + (p_x50 / p_slope) + plot_layout(widths = c(2, 1))
  print(final_slope_plot)
  
  # Save the slope figure to file
  dir.create(dirname(save_slope_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = save_slope_fig, plot = final_slope_plot, width = 10, height = 8, dpi = 300)
  
  #### Coefficient Plot (a, b, and c) AS A SEPARATE FIGURE ####
  # Get coefficient summary statistics (including p-values) from each nls model
  coeff_stats_list <- lapply(levels(data_clean$species), function(sp) {
    mod <- nls_models[[sp]]
    s <- summary(mod)$coefficients  # matrix with Estimate, Std. Error, t value, Pr(>|t|)
    df <- as.data.frame(s)
    df$Coefficient <- rownames(df)
    df$species <- sp
    df
  })
  coeff_stats <- bind_rows(coeff_stats_list)
  
  # Create a label: if p < 0.05, label is "*"; else "p = <pvalue>"
  coeff_stats <- coeff_stats %>%
    mutate(label = ifelse(`Pr(>|t|)` < 0.05,
                          "*",
                          sprintf("%.2f", `Pr(>|t|)`)))
  
  # Pivot the coefficient data (from coef_df) into long format and merge the labels
  coeffs_long <- coef_df %>%
    select(species, a, b, c) %>%
    pivot_longer(cols = c("a", "b", "c"),
                 names_to = "Coefficient",
                 values_to = "Value")
  
  coeffs_long <- left_join(coeffs_long, coeff_stats %>% select(species, Coefficient, label),
                           by = c("species", "Coefficient"))
  
  coeffs_long$species <- factor(coeffs_long$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  
  p_coeffs <- ggplot(coeffs_long, aes(x = species, y = Value, fill = species)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = label, y = Value/2),
              color = "black", size = 4,
              position = position_dodge(width = 0.9)) +
    scale_fill_manual(values = cb_palette) +
    facet_wrap(~ Coefficient, scales = "free_y") +
    labs(title = "Coefficients by Species for NDVI Quantiles with PSI",
         subtitle = expression(italic(NDVI) == a + b * e^(-c * italic(x))),
         x = "",
         y = "Coefficient Value") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      plot.subtitle = element_text(hjust = 0.5),
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
  
  print(p_coeffs)
  
  # Save the coefficient figure to file
  dir.create(dirname(save_coeff_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = save_coeff_fig, plot = p_coeffs, device = "png", width = 10, height = 8, dpi = 300)
}

plot_NDVI_Q_PSIbin_exp_slope_negPSI <- function(data, save_coeff_fig, save_slope_fig) {
  
  # Process data with your custom NDVI_PSIbin function and remove missing values
  data <- NDVI_PSIbin(data)
  data <- na.omit(data)
  
  # Load required libraries
  library(ggplot2)
  library(nlme)
  library(dplyr)
  library(tibble)
  library(patchwork)
  library(purrr)
  library(car)      # for deltaMethod
  library(tidyr)    # for pivot_longer
  
  # Identify the value column and order species; define the color palette
  value_col <- "avg_value"
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  data$species <- factor(data$species, levels = species_order)
  cb_palette <- c("Oak"   = "#E69F00",   # Orange
                  "Beech" = "#0072B2",   # Deep blue
                  "Spruce"= "#009E73",   # Bluish-green
                  "Pine"  = "#F0E442")   # Yellow
  
  # Keep original (negative) soil water potential
  data <- data %>% mutate(x = bin_median)
  
  # Clean data: remove rows with missing or non-finite values
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  
  # Set threshold manually to 11.5 instead of computing the median from the data
  threshold <- 11.5
  
  #### NONLINEAR MODELING PER SPECIES ####
  start_list <- list(a = 5, b = 3, c = 0.001)
  control_params <- nls.control(maxiter = 200, minFactor = 1e-4)
  
  nls_models <- nlsList(avg_value ~ a + b * exp(c * x) | species,
                        data = data_clean,
                        start = start_list,
                        control = control_params)
  print(summary(nls_models))
  
  # Extract coefficients
  coef_df <- as.data.frame(coef(nls_models), optional = TRUE) %>% 
    rownames_to_column(var = "species") %>%
    filter(!is.na(a))
  coef_df$species <- factor(coef_df$species, levels = species_order)
  
  # Compute x50 and slope at x50
  coef_df <- coef_df %>%
    mutate(x50 = ifelse((threshold - a) > 0 & b > 0,
                        log((threshold - a)/b) / c,
                        NA),
           slope50 = c * (threshold - a))
  
  #### PANEL A: Observed Data and Fitted Curves ####
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
  
  # No transformation — use natural negative axis
  x_scale <- scale_x_continuous()
  
  p_combined <- ggplot() +
    geom_point(data = data_clean, aes(x = x, y = avg_value, color = species)) +
    geom_line(data = pred_all, aes(x = x, y = pred, color = species), linewidth = 1) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "black", linewidth = 1) +
    annotate("text", x = min(data_clean$x), y = threshold, label = "median", 
             hjust = -0.1, vjust = -0.3, fontface = "italic", size = 6) +
    scale_color_manual(values = cb_palette) +
    x_scale +
    labs(x = "Soil Water Potential", y = "NDVI Quantiles") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(color = "black", size = 14),
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
  p_x50 <- ggplot(coef_df, aes(x = species, y = x50, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    labs(x = "", y = "Soil Water Potential") +
    scale_fill_manual(values = cb_palette) +
    expand_limits(y = 0) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank(),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(color = "black", size = 14),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(caption = "(b)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0),
          plot.caption.position = "plot")
  
  #### PANEL C: Bar Plot of Absolute Slope at x50 with Error Bars and Annotations ####
  stats_list <- lapply(levels(data_clean$species), function(sp) {
    mod <- nls_models[[sp]]
    coefs <- coef(mod)
    slope50_est <- -coefs["c"] * (threshold - coefs["a"])
    
    dm_result <- deltaMethod(mod, paste0("c*(a - ", threshold, ")"),
                             parameterNames = c("a", "b", "c"))
    se <- dm_result$SE
    df_mod <- summary(mod)$df[2]
    t_val <- slope50_est / se
    p_val <- 2 * (1 - pt(abs(t_val), df_mod))
    
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
  stats_df$species <- factor(stats_df$species, levels = species_order)
  
  stats_df <- stats_df %>%
    mutate(label_text = ifelse(p_val < 0.05,
                               sprintf("%.2f*", r_squared),
                               sprintf("%.2f\np = %.2f", r_squared, p_val)))
  
  p_slope <- ggplot(stats_df, aes(x = species, y = slope_abs, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, slope_abs - se), ymax = slope_abs + se), width = 0.2) +
    geom_text(aes(y = slope_abs/2, label = label_text), color = "black", size = 6) +
    labs(x = "", y = "Absolute Slope") +
    scale_fill_manual(values = cb_palette) +
    expand_limits(y = 0) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(color = "black", size = 14),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(caption = "(c)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0),
          plot.caption.position = "plot")
  
  #### COMBINE PANELS ####
  final_slope_plot <- p_combined + (p_x50 / p_slope) + plot_layout(widths = c(2, 1))
  print(final_slope_plot)
  dir.create(dirname(save_slope_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = save_slope_fig, plot = final_slope_plot, width = 10, height = 8, dpi = 300)
  
  #### Coefficient Plot ####
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
                          sprintf("%.2f", `Pr(>|t|)`)))
  
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
              color = "black", size = 4,
              position = position_dodge(width = 0.9)) +
    scale_fill_manual(values = cb_palette) +
    facet_wrap(~ Coefficient, scales = "free_y") +
    labs(title = "",
         # subtitle = expression(italic(NDVI) == a + b * e^{c * italic(x)}),
         x = "",
         y = "Coefficient Value") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      plot.subtitle = element_text(hjust = 0.5),
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
  
  print(p_coeffs)
  dir.create(dirname(save_coeff_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = save_coeff_fig, plot = p_coeffs, device = "png", width = 10, height = 8, dpi = 300)
}

plot_NDVI_Q_PSIbin_linear_exp_slope <- function(data, save_coeff_fig, save_slope_fig) {
  data <- NDVI_PSIbin(data)
  data <- na.omit(data)
  
  library(ggplot2)
  library(nlme)
  library(dplyr)
  library(tibble)
  library(patchwork)
  library(purrr)
  library(car)
  library(tidyr)
  
  value_col <- "avg_value"
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  data$species <- factor(data$species, levels = species_order)
  cb_palette <- c("Oak"   = "#E69F00", "Beech" = "#0072B2",
                  "Spruce"= "#009E73", "Pine"  = "#F0E442")
  
  data <- data %>% mutate(x = -bin_median)
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  threshold <- 11.5
  
  #### MODEL FITTING ####
  model_list <- list()
  for (sp in levels(data_clean$species)) {
    df_sp <- data_clean %>% filter(species == sp)
    if (sp == "Oak") {
      model_list[[sp]] <- lm(avg_value ~ x, data = df_sp)
    } else {
      model_list[[sp]] <- nls(avg_value ~ a + b * exp(-c * x),
                              data = df_sp,
                              start = list(a = 5, b = 3, c = 0.001),
                              control = nls.control(maxiter = 200, minFactor = 1e-4))
    }
  }
  
  #### COEFFICIENT EXTRACTION ####
  coef_df <- bind_rows(lapply(names(model_list), function(sp) {
    mod <- model_list[[sp]]
    if (inherits(mod, "lm")) {
      tibble(species = sp, a = coef(mod)[1], b = coef(mod)[2], c = NA)
    } else {
      coefs <- coef(mod)
      tibble(species = sp, a = coefs["a"], b = coefs["b"], c = coefs["c"])
    }
  }))
  coef_df$species <- factor(coef_df$species, levels = species_order)
  
  coef_df <- coef_df %>%
    mutate(x50 = ifelse(species == "Oak", (threshold - a) / b,
                        ifelse((threshold - a) > 0 & b > 0,
                               -log((threshold - a)/b) / c,
                               NA)),
           slope50 = ifelse(species == "Oak", b, -c * (threshold - a)))
  
  #### PANEL A: Observed Data and Fitted Curves ####
  pred_list <- lapply(levels(data_clean$species), function(sp) {
    df_sp <- data_clean %>% filter(species == sp)
    x_seq <- seq(min(df_sp$x, na.rm = TRUE), max(df_sp$x, na.rm = TRUE), length.out = 100)
    mod <- model_list[[sp]]
    pred <- predict(mod, newdata = data.frame(x = x_seq))
    data.frame(species = sp, x = x_seq, pred = pred)
  })
  pred_all <- bind_rows(pred_list)
  pred_all$species <- factor(pred_all$species, levels = species_order)
  
  x_scale <- scale_x_continuous(trans = "reverse", labels = function(x) -x)
  
  p_combined <- ggplot() +
    geom_point(data = data_clean, aes(x = x, y = avg_value, color = species)) +
    geom_line(data = pred_all, aes(x = x, y = pred, color = species), linewidth = 1) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "black", linewidth = 1) +
    annotate("text", x = 2000, y = threshold, label = "median", 
             hjust = -0.1, vjust = -0.3, fontface = "italic", size = 4) +
    scale_color_manual(values = cb_palette) +
    x_scale +
    labs(x = "Soil Water Potential", y = "NDVI Quantiles") +
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
      legend.position = "top",
      plot.caption = element_text(face = "bold", size = 16, hjust = 0),
      plot.caption.position = "plot"
    ) +
    coord_cartesian(clip = "off") +
    labs(caption = "(a)")
  
  #### PANEL B: Bar Plot of x50 ####
  p_x50 <- ggplot(coef_df, aes(x = species, y = -x50, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    labs(x = "", y = "Soil Water Potential") +
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
      panel.grid.minor = element_blank(),
      plot.caption = element_text(face = "bold", size = 16, hjust = 0),
      plot.caption.position = "plot"
    ) +
    labs(caption = "(b)")
  
  #### PANEL C: Slope Bar Plot ####
  stats_list <- lapply(levels(data_clean$species), function(sp) {
    mod <- model_list[[sp]]
    df_sp <- data_clean %>% filter(species == sp)
    
    if (inherits(mod, "lm")) {
      slope <- coef(mod)["x"]
      se <- summary(mod)$coefficients["x", "Std. Error"]
      p_val <- summary(mod)$coefficients["x", "Pr(>|t|)"]
      r_squared <- summary(mod)$r.squared
    } else {
      coefs <- coef(mod)
      slope <- -coefs["c"] * (threshold - coefs["a"])
      dm_result <- deltaMethod(mod, paste0("c*(a - ", threshold, ")"),
                               parameterNames = c("a", "b", "c"))
      se <- dm_result$SE
      df_mod <- summary(mod)$df[2]
      t_val <- slope / se
      p_val <- 2 * (1 - pt(abs(t_val), df_mod))
      fitted_vals <- predict(mod, newdata = df_sp)
      r_squared <- 1 - sum((df_sp[[value_col]] - fitted_vals)^2) /
        sum((df_sp[[value_col]] - mean(df_sp[[value_col]]))^2)
    }
    
    tibble(species = sp,
           slope50 = slope,
           slope_abs = abs(slope),
           se = se,
           p_val = p_val,
           r_squared = r_squared)
  })
  stats_df <- bind_rows(stats_list)
  stats_df$species <- factor(stats_df$species, levels = species_order)
  stats_df <- stats_df %>%
    mutate(label_text = ifelse(p_val < 0.05,
                               sprintf("%.2f*", r_squared),
                               sprintf("%.2f\np = %.2f", r_squared, p_val)))
  
  p_slope <- ggplot(stats_df, aes(x = species, y = slope_abs, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, slope_abs - se), ymax = slope_abs + se), width = 0.2) +
    geom_text(aes(y = slope_abs/2, label = label_text), size = 4, color = "black") +
    labs(x = "", y = "Absolute Slope") +
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
      panel.grid.minor = element_blank(),
      plot.caption = element_text(face = "bold", size = 16, hjust = 0),
      plot.caption.position = "plot"
    ) +
    labs(caption = "(c)")
  
  final_slope_plot <- p_combined + (p_x50 / p_slope) + plot_layout(widths = c(2, 1))
  print(final_slope_plot)
  dir.create(dirname(save_slope_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = save_slope_fig, plot = final_slope_plot, width = 10, height = 8, dpi = 300)
}

plot_NDVI_Q_PSIbin_log <- function(data, save_coeff_fig, save_slope_fig) {
  
  data <- NDVI_PSIbin_log(data)
  data <- na.omit(data)
  
  library(ggplot2)
  library(nlme)
  library(dplyr)
  library(tibble)
  library(patchwork)
  library(purrr)
  library(car)
  library(tidyr)
  
  value_col <- "avg_value"
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  data$species <- factor(data$species, levels = species_order)
  cb_palette <- c("Oak" = "#E69F00", "Beech" = "#0072B2", "Spruce" = "#009E73", "Pine" = "#F0E442")
  
  data <- data %>% mutate(x = -bin_median + 8) %>% filter(x > 0)
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x)) %>% mutate(percent_pixels = percentage * 100)
  threshold <- 11.5
  
  nls_models <- nlsList(avg_value ~ a * log(x) + b | species,
                        data = data_clean,
                        start = list(a = -1, b = 12),
                        control = nls.control(maxiter = 1000))
  
  coef_df <- coef(nls_models) %>% as.data.frame() %>% 
    rownames_to_column(var = "species") %>% 
    filter(!is.na(a))
  coef_df$species <- factor(coef_df$species, levels = species_order)
  
  coef_df <- coef_df %>% mutate(x50 = exp((threshold - b) / a), slope50 = a / x50)
  
  pred_list <- data_clean %>%
    group_by(species) %>%
    do({
      sp <- unique(.$species)
      x_seq <- seq(min(.$x), max(.$x), length.out = 100)
      pred <- predict(nls_models[[as.character(sp)]], newdata = data.frame(x = x_seq))
      data.frame(x = x_seq, pred = pred)
    })
  pred_all <- bind_rows(pred_list)
  
  stats_list <- lapply(levels(data_clean$species), function(sp) {
    mod <- nls_models[[sp]]
    coefs <- coef(mod)
    slope50_est <- coefs["a"] / exp((threshold - coefs["b"]) / coefs["a"])
    dm_result <- deltaMethod(mod, paste0("a / exp((", threshold, " - b)/a)"), parameterNames = c("a", "b"))
    se <- dm_result$SE
    df_mod <- summary(mod)$df[2]
    p_val <- 2 * (1 - pt(abs(slope50_est / se), df_mod))
    df_sp <- data_clean %>% filter(species == sp)
    fitted_vals <- predict(mod, newdata = df_sp)
    r_squared <- 1 - sum((df_sp[[value_col]] - fitted_vals)^2) / sum((df_sp[[value_col]] - mean(df_sp[[value_col]]))^2)
    tibble(species = sp, slope50 = slope50_est, slope_abs = abs(slope50_est), se = se, p_val = p_val, r_squared = r_squared)
  })
  stats_df <- bind_rows(stats_list)
  stats_df$species <- factor(stats_df$species, levels = species_order)
  
  # Panel A
  p_combined <- ggplot() +
    geom_point(data = data_clean, aes(x = x, y = avg_value, color = species, size = percent_pixels)) +
    geom_line(data = pred_all, aes(x = x, y = pred, color = species), linewidth = 1) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "black", linewidth = 1) +
    scale_color_manual(values = cb_palette, name = "") +
    scale_size_continuous(name = "Pixel %", range = c(1, 6)) +
    labs(x = "transformed soil water potential", y = "NDVI quantiles (rank)") +
    ggtitle("(a)") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0, vjust = 1, size = 16, face = "bold"),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(color = "black", size = 14),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.text = element_text(size = 14)
    )
  
  # Panel B
  p_x50 <- ggplot(coef_df, aes(x = species, y = x50, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = cb_palette, name = "") +
    guides(fill = "none") +
    labs(x = "", y = "transformed soil water potential") +
    ggtitle("(b)") +
    expand_limits(y = 0) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0, vjust = 1, size = 16, face = "bold"),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(color = "black", size = 14),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  # Panel C
  p_slope <- ggplot(stats_df, aes(x = species, y = slope_abs, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, slope_abs - se), ymax = slope_abs + se), width = 0.2) +
    geom_text(aes(label = if_else(p_val < 0.05, sprintf("%.2f*", r_squared), sprintf("%.2f", r_squared)), y = slope_abs/2), color = "black", size = 5) +
    scale_fill_manual(values = cb_palette, name = "") +
    guides(fill = "none") +
    labs(x = "", y = "absolute slope") +
    ggtitle("(c)") +
    expand_limits(y = 0) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0, vjust = 1, size = 16, face = "bold"),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(color = "black", size = 14),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  final_slope_plot <- (p_combined + (p_x50 / p_slope)) +
    plot_layout(widths = c(2, 1), guides = "collect") &
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 14),
      legend.key = element_rect(fill = "white", color = NA),
      legend.background = element_blank(),
      legend.box.background = element_blank()
    )
  
  ggsave(filename = save_slope_fig, plot = final_slope_plot, width = 10, height = 8, dpi = 300)
  
  # Coefficient plot
  coeff_stats_list <- lapply(levels(data_clean$species), function(sp) {
    mod <- nls_models[[sp]]
    s <- summary(mod)$coefficients
    df <- as.data.frame(s)
    df$Coefficient <- rownames(df)
    df$species <- sp
    df
  })
  coeff_stats <- bind_rows(coeff_stats_list)
  
  coeffs_long <- coef_df %>%
    select(species, a, b) %>%
    pivot_longer(cols = c("a", "b"), names_to = "Coefficient", values_to = "Value")
  
  coeffs_long <- left_join(coeffs_long, coeff_stats %>% select(species, Coefficient, `Std. Error`, `Pr(>|t|)`),
                           by = c("species", "Coefficient")) %>%
    mutate(label = if_else(`Pr(>|t|)` < 0.05, "*", sprintf("%.2f", `Pr(>|t|)`)))
  
  p_coeffs <- ggplot(coeffs_long, aes(x = species, y = Value, fill = species)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_errorbar(aes(ymin = Value - `Std. Error`, ymax = Value + `Std. Error`), position = position_dodge(width = 0.9), width = 0.2) +
    geom_text(aes(label = label, y = Value/2), position = position_dodge(width = 0.9), color = "black", size = 5) +
    facet_wrap(~ Coefficient, scales = "free_y") +
    scale_fill_manual(values = cb_palette, name = "") +
    labs(title = "NDVI quantiles ~ soil water potential",
         subtitle = expression(NDVI == a * log(x) + b + epsilon),
         x = "", y = "coefficient value") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 14),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(color = "black", size = 14),
      legend.position = "bottom",
      legend.text = element_text(size = 14),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    )
  
  ggsave(filename = save_coeff_fig, plot = p_coeffs, width = 10, height = 8, dpi = 300)
}

plot_NDVI_PDM_PSIbin_exp_slope <- function(data, save_coeff_fig, save_slope_fig) {
  
  # Process using your custom NDVI_PSIbin function and remove missing rows
  data <- NDVI_PSIbin(data)
  data <- na.omit(data)
  
  # Load required libraries
  library(ggplot2)
  library(nlme)
  library(dplyr)
  library(tibble)
  library(patchwork)
  library(purrr)
  library(car)      # for deltaMethod
  library(tidyr)    # for pivot_longer
  
  # Identify the value column and set species order + color palette
  value_col <- "avg_value"
  data$species <- factor(data$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  cb_palette <- c("Oak"   = "#E69F00",  # Orange
                  "Beech" = "#0072B2",  # Deep blue
                  "Spruce"= "#009E73",  # Bluish-green
                  "Pine"  = "#F0E442")  # Yellow
  
  # Convert soil water potential to positive values (x = -bin_median)
  data <- data %>% mutate(x = -bin_median)
  
  # Clean data: remove rows with missing or non-finite values
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  
  # Set threshold manually to 1 instead of using the median of avg_value
  threshold <- 1
  
  #### NONLINEAR MODELING PER SPECIES ####
  # Set starting parameters and control parameters for the nls models
  start_list <- list(a = 0.9, b = 0.03, c = 0.001)
  control_params <- nls.control(maxiter = 1000, minFactor = 1e-5)
  
  nls_models <- nlsList(avg_value ~ a + b * exp(-c * x) | species,
                        data = data_clean,
                        start = start_list,
                        control = control_params)
  
  # Extract coefficients into a data frame and compute x50 and slope at x50
  coef_df <- as.data.frame(coef(nls_models), optional = TRUE) %>% 
    rownames_to_column(var = "species") %>%
    filter(!is.na(a))
  coef_df$species <- factor(coef_df$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  
  coef_df <- coef_df %>%
    mutate(x50 = ifelse((threshold - a) > 0 & b > 0,
                        -log((threshold - a) / b) / c,
                        NA),
           slope50 = -c * (threshold - a))
  
  #### SLOPE FIGURE (Panels A, B, and C) ####
  # Panel A: Observed data and fitted curves
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
    geom_line(data = pred_all, aes(x = x, y = pred, color = species), linewidth = 1) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "black", linewidth = 1) +
    annotate("text", x = 2000, y = threshold, label = "median", 
             hjust = -0.1, vjust = -0.3, fontface = "italic", size = 4) +
    scale_color_manual(values = cb_palette) +
    scale_x_continuous(trans = "reverse", labels = function(x) -x) +
    labs(x = "Soil Water Potential", y = "NDVI Proportions") +
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
  p_x50 <- ggplot(coef_df, aes(x = species, y = -x50, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = cb_palette) +
    labs(y = "Soil Water Potential") +
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
  
  # Panel C: Bar plot of absolute slope at x50 with error bars and annotations
  stats_list <- lapply(levels(data_clean$species), function(sp) {
    mod <- nls_models[[sp]]
    coefs <- coef(mod)
    slope50_est <- -coefs["c"] * (threshold - coefs["a"])
    
    # Compute standard error for c*(a - threshold) using deltaMethod
    dm_result <- deltaMethod(mod, paste0("c*(a - ", threshold, ")"),
                             parameterNames = c("a", "b", "c"))
    se <- dm_result$SE
    df_mod <- summary(mod)$df[2]
    t_val <- slope50_est / se
    p_val <- 2 * (1 - pt(abs(t_val), df_mod))
    
    # Compute R² for the model for this species
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
  
  # Create annotation labels for Panel C:
  # If p < 0.05, display R² with a star; else, show R² and the p-value
  stats_df <- stats_df %>%
    mutate(label = ifelse(p_val < 0.05,
                          sprintf("%.2f*", r_squared),
                          sprintf("%.2f\np=%.3f", r_squared, p_val)))
  
  p_slope <- ggplot(stats_df, aes(x = species, y = slope_abs, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = cb_palette) +
    labs(y = "Absolute Slope") +
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
          plot.caption.position = "plot") +
    geom_text(aes(label = label, y = slope_abs/2), color = "black", size = 4)
  
  # Combine Panels A, B, and C into the slope figure
  final_plot <- p_combined + (p_x50 / p_slope) + plot_layout(widths = c(2, 1))
  print(final_plot)
  
  # Save the slope figure to file
  dir.create(dirname(save_slope_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = save_slope_fig, plot = final_plot, width = 10, height = 8, dpi = 300)
  
  #### COEFFICIENT FIGURE (Separate) ####
  # Get coefficient summary statistics for each species from the nls models.
  coeff_stats_list <- lapply(levels(data_clean$species), function(sp) {
    mod <- nls_models[[sp]]
    s <- summary(mod)$coefficients  # matrix with Estimate, Std. Error, t value, Pr(>|t|)
    df <- as.data.frame(s)
    df$Coefficient <- rownames(df)
    df$species <- sp
    df
  })
  coeff_stats <- bind_rows(coeff_stats_list)
  
  # Create a label: if p < 0.05 then "*", else "p=<pvalue>"
  coeff_stats <- coeff_stats %>%
    mutate(label = ifelse(`Pr(>|t|)` < 0.05,
                          "*",
                          sprintf("p=%.3f", `Pr(>|t|)`)))
  
  # Pivot the coefficient estimates from coef_df into long format and join with the labels
  coeffs_long <- coef_df %>%
    select(species, a, b, c) %>%
    pivot_longer(cols = c("a", "b", "c"),
                 names_to = "Coefficient",
                 values_to = "Value")
  
  coeffs_long <- left_join(coeffs_long, coeff_stats %>% select(species, Coefficient, label),
                           by = c("species", "Coefficient"))
  
  coeffs_long$species <- factor(coeffs_long$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
  
  p_coeffs <- ggplot(coeffs_long, aes(x = species, y = Value, fill = species)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = label, y = Value/2),
              color = "black", size = 4,
              position = position_dodge(width = 0.9)) +
    scale_fill_manual(values = cb_palette) +
    facet_wrap(~ Coefficient, scales = "free_y") +
    labs(title = "Coefficients by Species for NDVI Proportions",
         subtitle = expression(italic(NDVI) == a + b * e^(-c * italic(x))),
         x = "",
         y = "Coefficient Value") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      plot.subtitle = element_text(hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    )
  print(p_coeffs)
  
  # Save the coefficient figure to file
  dir.create(dirname(save_coeff_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = save_coeff_fig, plot = p_coeffs, device = "png", width = 10, height = 8, dpi = 300)
}

plot_TDiff_PSIbin_poly_3_slope <- function(data, coef_output, figure_output) {
  # Load required libraries (if not already loaded)
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
             hjust = -0.1, vjust = -0.3, fontface = "italic", size = 5) +
    geom_hline(yintercept = line_median, linetype = "dashed", color = "black", linewidth = 1) +
    annotate("text", x = min(TDiff_PSIbin_df$bin_median, na.rm = TRUE), 
             y = line_median, 
             label = paste0("median: ", round(line_median, 2)),
             hjust = -0.1, vjust = -0.3, fontface = "italic", size = 5) +
    scale_color_manual(values = cb_palette) +
    labs(x = "Soil Water Potential (bin_median)",
         y = "Average Transpiration Deficit") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
      axis.title = element_text(face = "bold", size = 20),
      axis.text = element_text(color = "black", size = 16),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 24, face = "bold", color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.text = element_text(size = 16),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 16)
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
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
      axis.title = element_text(face = "bold", size = 20),
      axis.text = element_text(color = "black", size = 16),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
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
  
  # Format the label for Panel C: if p < 0.05, append an asterisk to the R².
  local_slope_data <- local_slope_data %>%
    mutate(slope_abs = abs(slope),
           label_text = ifelse(p_value < 0.05,
                               sprintf("%.2f*", r2),
                               sprintf("p = %.2f\nR² = %.2f", p_value, r2)))
  
  p_bar <- ggplot(local_slope_data, aes(x = species, y = slope_abs, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = slope_abs - slope_se, ymax = slope_abs + slope_se),
                  width = 0.2, color = "black") +
    geom_text(aes(y = slope_abs/2, label = label_text), size = 5, color = "black") +
    scale_fill_manual(values = cb_palette) +
    labs(x = "", y = "Absolute Slope at Mean PSI") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
      axis.title = element_text(face = "bold", size = 20),
      axis.text = element_text(color = "black", size = 16),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  ## Combine Panels: Panel A on the left, Panels B (top) and C (bottom) on the right.
  right_panel <- p_bar_x / p_bar
  final_plot <- plot_mixed + right_panel + plot_layout(widths = c(2, 1))
  print(final_plot)
  
  # Save the combined plot with width = 10 and height = 8.
  ggsave(figure_output, plot = final_plot, width = 10, height = 8, dpi = 300)
  
  ## New Figure: Grouped Bar Plot for All Model Coefficients per Species
  # Extract the species-specific (conditional) coefficients.
  species_coef <- coef(model)$species
  
  # Convert rownames to a proper column and pivot the data to long format.
  coeff_data <- species_coef %>%
    tibble::rownames_to_column("species") %>%
    pivot_longer(cols = -species, names_to = "term", values_to = "value") %>%
    mutate(species = factor(species, levels = species_levels),
           term = dplyr::recode(term,
                                "(Intercept)" = "a",
                                "poly(bin_median, 3)1" = "b",
                                "poly(bin_median, 3)2" = "c",
                                "poly(bin_median, 3)3" = "d"))
  
  # Set the order of the term factor.
  coeff_data$term <- factor(coeff_data$term, levels = c("a", "b", "c", "d"))
  
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
    mutate(term = dplyr::recode(term,
                                "(Intercept)" = "a",
                                "poly(bin_median, 3)1" = "b",
                                "poly(bin_median, 3)2" = "c",
                                "poly(bin_median, 3)3" = "d")) %>%
    dplyr::select(species, term, p_value = `Pr(>|t|)`)
  
  # Join the p-values with the coefficient data.
  coeff_data <- left_join(coeff_data, coeff_stats, by = c("species", "term"))
  
  # Create a label: if p < 0.05 then "*" else the p_value with 2 decimals.
  coeff_data <- coeff_data %>%
    mutate(label_text = ifelse(p_value < 0.05, "*", sprintf("%.2f", p_value)))
  
  coeff_data$species <- factor(coeff_data$species, levels = species_levels)
  
  # Create the grouped bar plot with p-value labels centered in each bar.
  plot_coeff <- ggplot(coeff_data, aes(x = species, y = value, fill = species)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = label_text, y = value/2), 
              color = "black", size = 5, 
              position = position_dodge(width = 0.8)) +
    facet_wrap(~ term, scales = "free_y") +
    scale_fill_manual(values = cb_palette) +
    labs(title = "Coefficients of Transpiration Deficit (y) with Soil Water Potential (x)", 
         subtitle = expression(hat(Y) == a + b*x + c*x^2 + d*x^3),
         x = "Coefficient Term", 
         y = "Coefficient Value") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 18, face = "italic"),
      axis.title = element_text(face = "bold", size = 20),
      axis.text = element_text(color = "black", size = 16),
      panel.grid = element_blank(),
      legend.position = "top",
      legend.text = element_text(size = 16),
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(face = "bold", size = 16)
    )
  
  print(plot_coeff)
  
  # Save the coefficient plot with width = 10 and height = 8.
  ggsave(coef_output, plot = plot_coeff, width = 10, height = 8, dpi = 300)
  
  print(coeff_data)
  # Return both plots as a list (if desired)
  return(list(combined_plot = final_plot, coeff_plot = plot_coeff))
}

plot_NDVI_PSIbin_line <- function(data, figure_output = NULL) {
  # Process the data using your NDVI_PSIbin function
  df <- NDVI_PSIbin(data)
  
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
    labs(title = "NDVI vs PSI Bin",
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

plot_NDVI_PSIbin_linear <- function(data, save_slope_fig) {
  
  # Process data with your custom NDVI_PSIbin function and remove missing values
  data <- NDVI_PSIbin(data)
  data <- na.omit(data)
  
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
    labs(x = "Soil Water Potential", y = "NDVI Quantiles") +
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
    labs(x = "Species", y = "Soil Water Potential") +
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

plot_NDVI_PSIbin_slope_each <- function(data, save_coeff_fig, save_slope_fig) {
  
  # Process data with your custom NDVI_PSIbin function and remove missing values
  data <- NDVI_PSIbin(data)
  data <- na.omit(data)
  
  # Load required libraries
  library(ggplot2)
  library(nlme)
  library(dplyr)
  library(tibble)
  library(patchwork)
  library(purrr)
  library(car)      # (kept in case of future needs)
  library(tidyr)    # for pivot_longer
  library(ggpattern)  # for patterned bars in Panel C
  
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
  
  #### Compute species-specific medians ####
  medians_df <- data_clean %>% 
    group_by(species) %>% 
    summarise(median_value = median(.data[[value_col]]))
  
  #### Compute maximum and minimum soil water potential for each species ####
  max_x_df <- data_clean %>%
    group_by(species) %>%
    summarise(x_max = max(x))
  min_x_df <- data_clean %>%
    group_by(species) %>%
    summarise(x_min = min(x))
  
  #### NONLINEAR MODELING PER SPECIES ####
  start_list <- list(a = 5, b = 3, c = 0.001)
  control_params <- nls.control(maxiter = 200, minFactor = 1e-4)
  
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
  
  # Reverse the x-axis so that soil water potential is shown as positive
  x_scale <- scale_x_continuous(trans = "reverse", labels = function(x) -x)
  
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
    x_scale +
    labs(x = "Soil Water Potential", y = "NDVI Quantiles") +
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
  p_x50 <- ggplot(coef_df, aes(x = species, y = -x50, fill = species)) +
    geom_col(width = 0.7) +
    labs(x = "Species", y = "Soil Water Potential") +
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
    # geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
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
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0),
          plot.caption.position = "plot")
  
  #### Combine Panels ####
  # Arrange right side: Panel B on top and the combined Panel C below.
  right_side <- p_x50 / p_slope_combined
  final_slope_plot <- p_combined + right_side + plot_layout(widths = c(2, 1))
  print(final_slope_plot)
  
  # Save the slope figure to file
  dir.create(dirname(save_slope_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = save_slope_fig, plot = final_slope_plot, width = 10, height = 8, dpi = 300)
  
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
  dir.create(dirname(save_coeff_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = save_coeff_fig, plot = p_coeffs, device = "png", width = 10, height = 8, dpi = 300)
}

plot_NDVI_PDM_PSIbin_slope_each <- function(data, save_coeff_fig, save_slope_fig) {
  
  # Process data with your custom NDVI_PSIbin function and remove missing values
  data <- NDVI_PSIbin(data)
  data <- na.omit(data)
  
  # Load required libraries
  library(ggplot2)
  library(nlme)
  library(dplyr)
  library(tibble)
  library(patchwork)
  library(purrr)
  library(car)      # (kept in case of future needs)
  library(tidyr)    # for pivot_longer
  library(ggpattern)  # for patterned bars in Panel C
  
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
  
  #### Compute species-specific medians ####
  medians_df <- data_clean %>% 
    group_by(species) %>% 
    summarise(median_value = median(.data[[value_col]]))
  
  #### Compute maximum and minimum soil water potential for each species ####
  max_x_df <- data_clean %>%
    group_by(species) %>%
    summarise(x_max = max(x))
  min_x_df <- data_clean %>%
    group_by(species) %>%
    summarise(x_min = min(x))
  
  #### NONLINEAR MODELING PER SPECIES ####
  start_list <- list(a = 1, b = 0.01, c = 0.001)
  control_params <- nls.control(maxiter = 200, minFactor = 1e-4)
  
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
  
  # Reverse the x-axis so that soil water potential is shown as positive
  x_scale <- scale_x_continuous(trans = "reverse", labels = function(x) -x)
  
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
    x_scale +
    labs(x = "Soil Water Potential", y = "NDVI Quantiles") +
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
  p_x50 <- ggplot(coef_df, aes(x = species, y = -x50, fill = species)) +
    geom_col(width = 0.7) +
    labs(x = "Species", y = "Soil Water Potential") +
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
    # geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
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
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0),
          plot.caption.position = "plot")
  
  #### Combine Panels ####
  # Arrange right side: Panel B on top and the combined Panel C below.
  right_side <- p_x50 / p_slope_combined
  final_slope_plot <- p_combined + right_side + plot_layout(widths = c(2, 1))
  print(final_slope_plot)
  
  # Save the slope figure to file
  dir.create(dirname(save_slope_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = save_slope_fig, plot = final_slope_plot, width = 10, height = 8, dpi = 300)
  
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
  dir.create(dirname(save_coeff_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = save_coeff_fig, plot = p_coeffs, device = "png", width = 10, height = 8, dpi = 300)
}

plot_TDiff_PSIbin_slope_each <- function(NDVI_PSI_TDiff_species, save_coeff_fig, save_slope_fig) {
  # Load required libraries
  library(dplyr)      # Data manipulation
  library(ggplot2)    # Plotting
  library(tidyr)      # Data reshaping
  library(patchwork)  # Combining plots
  library(ggpattern)  # For patterned bars
  library(purrr)      # For mapping
  
  # Process the data using your custom function and remove missing values.
  TDiff_PSIbin_df <- TDiff_PSIbin(NDVI_PSI_TDiff_species)
  TDiff_PSIbin_df <- na.omit(TDiff_PSIbin_df)
  
  # Define a custom color palette and species order.
  species_levels <- c("Oak", "Beech", "Spruce", "Pine")
  cb_palette <- c("Oak"   = "#E69F00",
                  "Beech" = "#0072B2",
                  "Spruce"= "#009E73",
                  "Pine"  = "#F0E442")
  TDiff_PSIbin_df$species <- factor(TDiff_PSIbin_df$species, levels = species_levels)
  
  # Create a positive soil water potential variable (x = -bin_median).
  TDiff_PSIbin_df <- TDiff_PSIbin_df %>% mutate(x = -bin_median)
  
  #### Compute Species-specific Summary Statistics ####
  medians_df <- TDiff_PSIbin_df %>% 
    group_by(species) %>% 
    summarise(median_deficit = median(avg_transpiration_deficit))
  
  max_x_df <- TDiff_PSIbin_df %>% 
    group_by(species) %>% 
    summarise(x_max = max(x))
  min_x_df <- TDiff_PSIbin_df %>% 
    group_by(species) %>% 
    summarise(x_min = min(x))
  
  #### Fit Separate Quadratic Models for Each Species ####
  # Use a quadratic model: response ~ bin_median + I(bin_median^2)
  models_df <- TDiff_PSIbin_df %>% 
    group_by(species) %>% 
    do(model = lm(avg_transpiration_deficit ~ bin_median + I(bin_median^2), data = .))
  
  #### Create Prediction Data ####
  pred_data <- TDiff_PSIbin_df %>%
    group_by(species) %>%
    summarise(min_bin = min(bin_median), max_bin = max(bin_median)) %>%
    group_by(species) %>%
    do(data.frame(bin_median = seq(.$min_bin, .$max_bin, length.out = 100))) %>%
    ungroup() %>%
    mutate(x = -bin_median)
  pred_data$species <- factor(pred_data$species, levels = species_levels)
  
  # Get predicted values by joining the species-specific model for each row.
  pred_data <- pred_data %>%
    left_join(models_df, by = "species") %>%
    rowwise() %>%
    mutate(predicted = predict(model, newdata = data.frame(bin_median = bin_median))) %>%
    ungroup() %>%
    select(-model)
  
  #### Identify the Median Point (x50) ####
  median_points <- pred_data %>% 
    left_join(medians_df, by = "species") %>%
    group_by(species) %>%
    slice(which.min(abs(predicted - median_deficit))) %>%
    ungroup() %>%
    rename(x50 = x) %>%  
    select(species, x50)
  
  # Combine summary statistics.
  coef_df <- medians_df %>% 
    left_join(max_x_df, by = "species") %>%
    left_join(min_x_df, by = "species") %>%
    left_join(median_points, by = "species")
  
  #### Determine Endpoints from the Prediction Data ####
  endpoints_max <- pred_data %>% 
    left_join(max_x_df, by = "species") %>%
    group_by(species) %>%
    slice(which.min(abs(x - x_max))) %>%
    ungroup() %>%
    rename(x_upper = x, y_upper = predicted) %>%
    select(species, x_upper, y_upper)
  
  endpoints_min <- pred_data %>% 
    left_join(min_x_df, by = "species") %>%
    group_by(species) %>%
    slice(which.min(abs(x - x_min))) %>%
    ungroup() %>%
    rename(x_lower = x, y_lower = predicted) %>%
    select(species, x_lower, y_lower)
  
  coef_df <- coef_df %>% 
    left_join(endpoints_max, by = "species") %>%
    left_join(endpoints_min, by = "species")
  
  #### Create Dashed Lines for Linear Segments ####
  linear_lines_max <- coef_df %>%
    pivot_longer(cols = c(x50, x_upper), names_to = "point", values_to = "x") %>%
    mutate(y = ifelse(point == "x50", median_deficit, y_upper))
  
  linear_lines_min <- coef_df %>%
    pivot_longer(cols = c(x50, x_lower), names_to = "point", values_to = "x") %>%
    mutate(y = ifelse(point == "x50", median_deficit, y_lower))
  
  #### Panel A: Plot Observed Data, Fitted Curves, Dashed Regression Lines, and Emphasized Median ####
  p_combined <- ggplot() +
    geom_point(data = TDiff_PSIbin_df, aes(x = x, y = avg_transpiration_deficit, color = species)) +
    geom_line(data = pred_data, aes(x = x, y = predicted, color = species), linewidth = 1) +
    geom_line(data = linear_lines_max, aes(x = x, y = y, group = species), 
              color = "black", linetype = "dashed", size = 0.5) +
    geom_line(data = linear_lines_min, aes(x = x, y = y, group = species), 
              color = "black", linetype = "dashed", size = 0.5) +
    geom_point(data = coef_df, aes(x = x50, y = median_deficit), 
               shape = 10, size = 5, color = "black") +
    scale_color_manual(values = cb_palette) +
    scale_x_continuous(trans = "reverse", labels = function(x) -x) +
    labs(x = "Soil Water Potential", y = "Average Transpiration Deficit") +
    theme_minimal() +
    labs(caption = "(a)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top",
          plot.background = element_rect(fill = "white", color = "white"),
          panel.background = element_rect(fill = "white"))
  
  #### Panel B: Bar Plot of X-axis Values (x50) ####
  p_x50 <- ggplot(coef_df, aes(x = species, y = -x50, fill = species)) +
    geom_col(width = 0.7) +
    labs(x = "Species", y = "Soil Water Potential") +
    scale_fill_manual(values = cb_palette) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          axis.title.x = element_blank(),
          axis.text = element_text(color = "black"),
          plot.background = element_rect(fill = "white", color = "white"),
          panel.background = element_rect(fill = "white"),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(caption = "(b)") +
    theme(plot.caption = element_text(face = "bold", size = 14, hjust = 0),
          plot.caption.position = "plot")
  
  #### Panel C: Bar Plot of Local Slopes (Two Regressions) with Patterns ####
  stats_list_max <- lapply(levels(TDiff_PSIbin_df$species), function(sp) {
    sp_coef <- coef_df %>% filter(species == sp)
    if(is.na(sp_coef$x50) | is.na(sp_coef$x_upper)) {
      return(tibble(species = sp, slope = NA, se = NA, p_val = NA))
    }
    df_sp <- pred_data %>% filter(species == sp)
    lower_bound <- min(sp_coef$x50, sp_coef$x_upper)
    upper_bound <- max(sp_coef$x50, sp_coef$x_upper)
    df_subset <- df_sp %>% filter(x >= lower_bound, x <= upper_bound)
    if(nrow(df_subset) < 2) {
      x_points <- c(sp_coef$x50, sp_coef$x_upper)
      y_points <- c(sp_coef$median_deficit, sp_coef$y_upper)
      lm_fit <- lm(y_points ~ x_points)
    } else {
      lm_fit <- lm(predicted ~ x, data = df_subset)
    }
    slope_val <- coef(lm_fit)["x"]
    se <- summary(lm_fit)$coefficients["x", "Std. Error"]
    p_val <- summary(lm_fit)$coefficients["x", "Pr(>|t|)"]
    tibble(species = sp, slope = slope_val, se = se, p_val = p_val, RegType = "Max")
  })
  stats_df_max <- bind_rows(stats_list_max)
  
  stats_list_min <- lapply(levels(TDiff_PSIbin_df$species), function(sp) {
    sp_coef <- coef_df %>% filter(species == sp)
    if(is.na(sp_coef$x50) | is.na(sp_coef$x_lower)) {
      return(tibble(species = sp, slope = NA, se = NA, p_val = NA))
    }
    df_sp <- pred_data %>% filter(species == sp)
    lower_bound <- min(sp_coef$x50, sp_coef$x_lower)
    upper_bound <- max(sp_coef$x50, sp_coef$x_lower)
    df_subset <- df_sp %>% filter(x >= lower_bound, x <= upper_bound)
    if(nrow(df_subset) < 2) {
      x_points <- c(sp_coef$x_lower, sp_coef$x50)
      y_points <- c(sp_coef$y_lower, sp_coef$median_deficit)
      lm_fit <- lm(y_points ~ x_points)
    } else {
      lm_fit <- lm(predicted ~ x, data = df_subset)
    }
    slope_val <- coef(lm_fit)["x"]
    se <- summary(lm_fit)$coefficients["x", "Std. Error"]
    p_val <- summary(lm_fit)$coefficients["x", "Pr(>|t|)"]
    tibble(species = sp, slope = slope_val, se = se, p_val = p_val, RegType = "Min")
  })
  stats_df_min <- bind_rows(stats_list_min)
  
  stats_df_combined <- bind_rows(stats_df_max, stats_df_min) %>%
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
          axis.text = element_text(color = "black"),
          plot.background = element_rect(fill = "white", color = "white"),
          panel.background = element_rect(fill = "white"),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(caption = "(c)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0),
          plot.caption.position = "plot")
  
  #### Combine Panels ####
  final_slope_plot <- p_combined + (p_x50 / p_slope_combined) + plot_layout(widths = c(2, 1))
  print(final_slope_plot)
  
  dir.create(dirname(save_slope_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = save_slope_fig, plot = final_slope_plot, width = 10, height = 8, dpi = 300)
  
  #### Coefficient Plot (Separate Figure) ####
  # Extract coefficients and p-values from each species-specific quadratic model.
  species_coefs <- models_df %>%
    rowwise() %>%
    mutate(coefs = list(summary(model)$coefficients)) %>%
    ungroup() %>%
    mutate(species = as.character(species)) %>%
    select(species, coefs)
  
  # Rename the coefficient columns so they can be pivoted.
  species_coefs_long <- species_coefs %>%
    mutate(a_est = map_dbl(coefs, ~ .x["(Intercept)", "Estimate"]),
           a_pval = map_dbl(coefs, ~ .x["(Intercept)", "Pr(>|t|)"]),
           b_est = map_dbl(coefs, ~ .x["bin_median", "Estimate"]),
           b_pval = map_dbl(coefs, ~ .x["bin_median", "Pr(>|t|)"]),
           c_est = map_dbl(coefs, ~ .x["I(bin_median^2)", "Estimate"]),
           c_pval = map_dbl(coefs, ~ .x["I(bin_median^2)", "Pr(>|t|)"])) %>%
    select(species, a_est, a_pval, b_est, b_pval, c_est, c_pval) %>%
    pivot_longer(cols = -species, names_to = c("Coefficient", ".value"), names_sep = "_")
  
  # Ensure the species order is the same.
  species_coefs_long$species <- factor(species_coefs_long$species, levels = species_levels)
  
  p_coeffs <- ggplot(species_coefs_long, aes(x = species, y = est, fill = species)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = ifelse(pval < 0.05, "*", sprintf("%.3f", pval)), y = est/2), 
              color = "black", size = 3.5, position = position_dodge(width = 0.9)) +
    scale_fill_manual(values = cb_palette) +
    facet_wrap(~ Coefficient, scales = "free_y") +
    labs(title = "Coefficients by Species for Transpiration Deficit",
         x = "", y = "Coefficient Value") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),  # Center the title
      axis.text.x = element_text(hjust = 1, size = 10),
      axis.text = element_text(color = "black", size = 10),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.position = "top",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  print(p_coeffs)
  
  dir.create(dirname(save_coeff_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = save_coeff_fig, plot = p_coeffs, device = "png", width = 10, height = 8, dpi = 300)
  
  return(list(slope_plot = final_slope_plot, coeff_plot = p_coeffs))
}

plot_TDiff_PSIbin_poly_2_slope <- function(data, coef_output, figure_output) {
  # Load required libraries
  library(lme4)      # For mixed-effects modeling
  library(dplyr)     # For data manipulation
  library(ggplot2)   # For plotting
  library(patchwork) # For combining ggplots
  library(tidyr)     # For pivoting data
  library(car)       # For deltaMethod
  library(nlme)      # For nlsList (if needed)
  library(tibble)    # For rownames_to_column
  
  # Process the data and remove NAs
  TDiff_PSIbin_df <- TDiff_PSIbin(data)
  TDiff_PSIbin_df <- na.omit(TDiff_PSIbin_df)
  
  # Define color palette and species order
  cb_palette <- c("Oak"   = "#E69F00",   # Orange
                  "Beech" = "#0072B2",   # Deep blue
                  "Spruce"= "#009E73",   # Bluish-green
                  "Pine"  = "#F0E442")   # Yellow
  species_levels <- c("Oak", "Beech", "Spruce", "Pine")
  TDiff_PSIbin_df$species <- factor(TDiff_PSIbin_df$species, levels = species_levels)
  
  ## Panel A: Mixed-Effects Model Plot
  model <- lmer(avg_transpiration_deficit ~ poly(bin_median, 2, raw = TRUE) + 
                  (poly(bin_median, 2, raw = TRUE) | species),
                data = TDiff_PSIbin_df)
  
  # Create prediction data
  pred_data <- TDiff_PSIbin_df %>%
    group_by(species) %>%
    summarise(min_bin = min(bin_median),
              max_bin = max(bin_median)) %>%
    group_by(species) %>%
    do(data.frame(bin_median = seq(.$min_bin, .$max_bin, length.out = 100))) %>%
    ungroup()
  pred_data$species <- factor(pred_data$species, levels = species_levels)
  pred_data$predicted <- predict(model, newdata = pred_data, re.form = NULL)
  
  # Keep only the segment from maximum fitted value
  pred_data <- pred_data %>%
    group_by(species) %>%
    filter(bin_median >= bin_median[which.max(predicted)]) %>%
    ungroup()
  
  # Calculate mean and median transpiration deficit
  line_val <- mean(data$transpiration_deficit)
  line_median <- median(data$transpiration_deficit)
  
  # Create the mixed-effects plot
  plot_mixed <- ggplot(TDiff_PSIbin_df, aes(x = bin_median, y = avg_transpiration_deficit, color = species)) +
    geom_point() +
    geom_line(data = pred_data, aes(x = bin_median, y = predicted, color = species), linewidth = 1) +
    geom_hline(yintercept = line_val, linetype = "dashed", color = "black", linewidth = 1) +
    annotate("text", x = min(TDiff_PSIbin_df$bin_median, na.rm = TRUE), 
             y = line_val, label = paste0("mean: ", round(line_val, 2)),
             hjust = -0.1, vjust = -0.3, fontface = "italic", size = 5) +
    geom_hline(yintercept = line_median, linetype = "dashed", color = "black", linewidth = 1) +
    annotate("text", x = min(TDiff_PSIbin_df$bin_median, na.rm = TRUE), 
             y = line_median, label = paste0("median: ", round(line_median, 2)),
             hjust = -0.1, vjust = -0.3, fontface = "italic", size = 5) +
    scale_color_manual(values = cb_palette, name = "") +
    labs(x = "soil water potential (kPa)", y = "transpiration deficit (mm)") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      plot.caption = element_text(face = "bold", size = 16, hjust = 0),
      plot.caption.position = "plot",
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(color = "black", size = 14),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      legend.text = element_text(size = 14) 
    ) +
    labs(caption = "(a)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0.5)) +
    coord_cartesian(clip = "off")
  
  ## Analytical Calculation of x50 and Slope (using delta method)
  species_list <- species_levels
  stats_list <- lapply(species_list, function(sp) {
    df_sp <- TDiff_PSIbin_df %>% filter(species == sp)
    mod_sp <- lm(avg_transpiration_deficit ~ poly(bin_median, 2, raw = TRUE), data = df_sp)
    coefs <- coef(mod_sp)
    a <- coefs[1]; b <- coefs[2]; c <- coefs[3]
    disc <- b^2 - 4 * c * (a - line_val)
    x50 <- ifelse(disc >= 0, (-b + sqrt(disc)) / (2 * c), NA_real_)
    if(!is.na(x50)) {
      slope_expr <- paste0("b + 2 * c * ", x50)
      dm_result <- deltaMethod(mod_sp, slope_expr, parameterNames = c("a", "b", "c"))
      t_val <- (b + 2 * c * x50) / dm_result$SE
      df_resid <- df.residual(mod_sp)
      p_val <- 2 * (1 - pt(abs(t_val), df_resid))
      fitted_vals <- predict(mod_sp)
      r_squared <- 1 - sum((df_sp$avg_transpiration_deficit - fitted_vals)^2) /
        sum((df_sp$avg_transpiration_deficit - mean(df_sp$avg_transpiration_deficit))^2)
      tibble(
        species = sp,
        x50 = x50,
        slope50 = b + 2 * c * x50,
        slope_abs = abs(b + 2 * c * x50),
        se = dm_result$SE,
        p_val = p_val,
        r_squared = r_squared
      )
    } else {
      tibble(
        species = sp,
        x50 = NA_real_,
        slope50 = NA_real_,
        slope_abs = NA_real_,
        se = NA_real_,
        p_val = NA_real_,
        r_squared = NA_real_
      )
    }
  })
  
  stats_df <- bind_rows(stats_list)
  stats_df$species <- factor(stats_df$species, levels = species_levels)
  
  # Create label text with R² and significance
  stats_df <- stats_df %>%
    mutate(label_text = ifelse(p_val < 0.05,
                               sprintf("%.2f*", r_squared),
                               sprintf("%.2f\np = %.2f", r_squared, p_val)))
  
  ## Panel B: Bar Plot of x50 Values
  p_bar_x <- ggplot(stats_df, aes(x = species, y = x50, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = cb_palette) +
    labs(x = "", y = "soil water potential (kPa)") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(color = "black", size = 14),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(caption = "(b)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0.5),
          plot.caption.position = "plot")
  
  ## Panel C: Bar Plot of Absolute Slopes with Error Bars (using delta method SE)
  p_bar <- ggplot(stats_df, aes(x = species, y = slope_abs, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, slope_abs - se), 
                      ymax = slope_abs + se), 
                  width = 0.2) +
    geom_text(aes(label = label_text, y = slope_abs/2), 
              color = "black", size = 5) +
    scale_fill_manual(values = cb_palette) +
    labs(x = "", y = "absolute slope") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(color = "black", size = 14),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(caption = "(c)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0.5),
          plot.caption.position = "plot")
  
  # Combine panels (Panel A + Panel B/C)
  final_plot <- (plot_mixed + (p_bar_x / p_bar)) +
    plot_layout(widths = c(2, 1), guides = "collect") &
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 14),
      legend.key = element_rect(fill = "white", color = NA),
      legend.background = element_blank(),
      legend.box.background = element_blank()
    )
  
  print(final_plot)
  ggsave(figure_output, plot = final_plot, width = 10, height = 8, dpi = 300)
  
  ## Coefficient Plot: Extract and visualize species-specific (conditional) coefficients
  # Extract species-specific coefficients from the mixed model
  species_coef <- coef(model)$species
  
  # Convert rownames to a proper column and pivot the data to long format.
  coeff_data <- species_coef %>%
    tibble::rownames_to_column("species") %>%
    pivot_longer(cols = -species, names_to = "term", values_to = "value") %>%
    mutate(species = factor(species, levels = species_levels),
           term = dplyr::recode(term,
                                "(Intercept)" = "a",
                                "poly(bin_median, 2, raw = TRUE)1" = "b",
                                "poly(bin_median, 2, raw = TRUE)2" = "c"))
  
  # Set the order of the term factor.
  coeff_data$term <- factor(coeff_data$term, levels = c("a", "b", "c"))
  
  ## Compute p-values for coefficients by fitting separate linear models for each species.
  coeff_stats_list <- lapply(species_levels, function(sp) {
    subdata <- subset(TDiff_PSIbin_df, species == sp)
    mod_sp <- lm(avg_transpiration_deficit ~ poly(bin_median, 2, raw = TRUE), data = subdata)
    summ <- summary(mod_sp)$coefficients
    df <- as.data.frame(summ)
    df$term <- rownames(df)
    df$species <- sp
    df
  })
  
  coeff_stats <- do.call(rbind, coeff_stats_list)
  coeff_stats <- coeff_stats %>%
    mutate(term = dplyr::recode(term,
                                "(Intercept)" = "a",
                                "poly(bin_median, 2, raw = TRUE)1" = "b",
                                "poly(bin_median, 2, raw = TRUE)2" = "c")) %>%
    dplyr::select(species, term, p_value = `Pr(>|t|)`)
  
  # Join the p-values with the coefficient data.
  coeff_data <- left_join(coeff_data, coeff_stats, by = c("species", "term"))
  
  # Create a label: if p < 0.05 then "*" else the p-value with 2 decimals.
  coeff_data <- coeff_data %>%
    mutate(label_text = ifelse(p_value < 0.05, "*", sprintf("%.2f", p_value)))
  
  coeff_data$species <- factor(coeff_data$species, levels = species_levels)
  
  # Create the grouped bar plot with p-value labels centered in each bar.
  plot_coeff <- ggplot(coeff_data, aes(x = species, y = value, fill = species)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = label_text, y = value/2), 
              color = "black", size = 5, 
              position = position_dodge(width = 0.8)) +
    facet_wrap(~ term, scales = "free_y") +
    scale_fill_manual(values = cb_palette, name = "") +
    labs(title = "transpiration deficit ~ soil water potential", 
         subtitle = expression(NDVI == a + b*x + c*x^2 + epsilon),
         x = "", 
         y = "coefficient value") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      plot.subtitle = element_text(hjust = 0.5, size =14),
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
  
  print(plot_coeff)
  
  # Save the coefficient plot with width = 10 and height = 8.
  ggsave(coef_output, plot = plot_coeff, width = 10, height = 8, dpi = 300)
  
  return(list(combined_plot = final_plot, stats_df = stats_df, coeff_plot = plot_coeff))
}
