plot_time_series_Quantiles_PSI_TDiff_year_month_species <- function(data, figure_output) {
  # Load required libraries
  library(ggplot2)
  library(dplyr)
  library(lubridate)
  library(gridExtra)
  library(grid)  # For textGrob if needed
  
  # Define species order and update the data accordingly
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  data$species <- factor(data$species, levels = species_order)
  
  # Define the color-blind palette with the specified order
  cb_palette <- c("Oak"   = "#E69F00",   # Orange
                  "Beech" = "#0072B2",   # Deep blue
                  "Spruce"= "#009E73",   # Bluish-green
                  "Pine"  = "#F0E442")   # Yellow
  
  # Filter data for years 2003-2024 and months April to August
  data_filtered <- data %>%
    filter(year >= 2003, year <= 2024,
           month %in% c("April", "May", "June", "July", "August"))
  
  # Aggregate data: take the mean value for each month in each year for each species
  data_summary <- data_filtered %>%
    group_by(year, month, species) %>%
    summarise(
      Quantiles = mean(Quantiles, na.rm = TRUE),
      soil_water_potential = mean(soil_water_potential, na.rm = TRUE),
      transpiration_deficit = mean(transpiration_deficit, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    # Create a date variable (assuming the 15th day of each month)
    mutate(date = as.Date(paste(year, month, "15", sep = "-"), format = "%Y-%B-%d"))
  
  # Define a common x-axis scale with breaks for every year
  common_x_scale <- scale_x_date(
    breaks = seq(as.Date("2003-01-01"), as.Date("2024-12-31"), by = "year"),
    date_labels = "%Y"
  )
  
  # Define a common custom theme for all plots
  common_theme <- theme(
    plot.background = element_rect(fill = "white", color = "white"),
    panel.background = element_rect(fill = "white"),
    legend.background = element_rect(fill = "white", color = "white"),
    axis.title = element_text(face = "bold", size = 20),  # Larger axis titles
    axis.text = element_text(color = "black", size = 16), # Larger axis text
    legend.text = element_text(size = 18),  # Larger legend text
    legend.title = element_text(size = 20, face = "bold"), # Larger legend title
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.5),
    panel.grid.minor.y = element_blank()
  )
  
  # Create Panel (a): Time series for Mean Quantiles with thicker lines
  plot_a <- ggplot(data_summary, aes(x = date, y = Quantiles, color = species)) +
    geom_line(linewidth = 1.2) +  # Thicker lines
    labs(x = "", y = "Quantiles") +
    scale_color_manual(values = cb_palette) +
    common_x_scale +
    common_theme +
    theme(legend.position = "top")
  
  # Create Panel (b): Time series for Mean Soil Water Potential with thicker lines
  plot_b <- ggplot(data_summary, aes(x = date, y = soil_water_potential, color = species)) +
    geom_line(linewidth = 1.2) +
    labs(x = "", y = "Soil Water Potential") +
    scale_color_manual(values = cb_palette) +
    common_x_scale +
    common_theme +
    theme(legend.position = "none")
  
  # Create Panel (c): Time series for Mean Transpiration Deficit with thicker lines
  plot_c <- ggplot(data_summary, aes(x = date, y = transpiration_deficit, color = species)) +
    geom_line(linewidth = 1.2) +
    labs(x = "Year", y = "Transpiration Deficit") +
    scale_color_manual(values = cb_palette) +
    common_x_scale +
    common_theme +
    theme(legend.position = "none",
          axis.text.x = element_text(hjust = 1, size = 16))  # Larger x-axis text
  
  # Arrange the three plots vertically (each showing full axes)
  final_plot <- grid.arrange(plot_a, plot_b, plot_c, ncol = 1)
  
  # Save the final plot to the specified output file with high resolution
  ggsave(filename = figure_output, plot = final_plot, width = 20, height = 15, dpi = 300)
  
  return(final_plot)
}

NDVI_PSIbin_month <- function(df, bin_width = 100) {
  library(dplyr)
  
  # Identify the correct column dynamically (here we use "Quantiles")
  value_column <- if ("Quantiles" %in% names(df)) "Quantiles" else "Proportions"
  
  # Define bin breaks dynamically based on soil water potential range
  psi_min <- floor(min(df$soil_water_potential, na.rm = TRUE))
  psi_max <- ceiling(max(df$soil_water_potential, na.rm = TRUE))
  bin_breaks <- seq(psi_min, psi_max, by = bin_width)
  
  df <- df %>%
    mutate(PSI_bin = cut(soil_water_potential, breaks = bin_breaks, include.lowest = TRUE, right = FALSE))
  
  # Compute mean for the chosen value column per species, per month, and per PSI_bin
  meanNDVI_PSIbin_species <- df %>%
    group_by(species, month, PSI_bin) %>%
    summarise(
      avg_value = mean(.data[[value_column]], na.rm = TRUE),
      count = n(),
      .groups = 'drop'
    ) %>%
    filter(count >= 2000) %>%  # Only keep bins with at least 2000 observations
    mutate(bin_median = sapply(as.character(PSI_bin), function(bin_label) {
      nums <- as.numeric(strsplit(gsub("\\[|\\]|\\(|\\)", "", bin_label), ",")[[1]])
      mean(nums)
    })) %>%
    select(species, month, PSI_bin, bin_median, avg_value)
  
  return(meanNDVI_PSIbin_species)
}

plot_Quantiles_PSIbin_month_species <- function(data, bin_width = 100, figure_output = NULL) {
  library(ggplot2)
  library(dplyr)
  library(viridis)  # For color-blind-friendly colors
  
  # Define the order of species and months
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  month_order <- c("April", "May", "June", "July", "August")
  
  # Convert species and month to factors with defined order
  data$species <- factor(data$species, levels = species_order)
  data$month <- factor(data$month, levels = month_order)
  
  # Aggregate the data
  aggregated_data <- NDVI_PSIbin_month(data, bin_width)
  
  # Convert month to factor for correct ordering in the plot
  aggregated_data$month <- factor(aggregated_data$month, levels = month_order)
  
  # Create the plot using a color-blind-friendly palette
  p <- ggplot(aggregated_data, aes(x = bin_median, y = avg_value, color = month, group = month)) +
    geom_line(linewidth = 1) +  # Make lines more visible
    geom_point(size = 2) +
    facet_wrap(~ species, scales = "free_x") +
    scale_color_viridis_d(option = "C", begin = 0.2, end = 0.8) +  # Color-blind friendly
    labs(
      x = "Soil Water Potential (Bin Median)",
      y = "Mean Quantiles",
      color = "Month",
      title = "Mean Quantiles vs Soil Water Potential by Month"
    ) +
    theme_minimal() +
    theme(
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
  
  # Save the plot if an output filename is provided
  if (!is.null(figure_output)) {
    ggsave(filename = figure_output, plot = p, width = 10, height = 8)
  }
  
  return(p)
}

NDVI_TDiffbin_month <- function(df, bin_width = 5) {
  library(dplyr)
  
  # Identify the correct column dynamically (either "Quantiles" or "Proportions")
  value_column <- if ("Quantiles" %in% names(df)) "Quantiles" else "Proportions"
  
  # Define bin breaks dynamically based on transpiration deficit range
  tdiff_min <- floor(min(df$transpiration_deficit, na.rm = TRUE))
  tdiff_max <- ceiling(max(df$transpiration_deficit, na.rm = TRUE))
  bin_breaks <- seq(tdiff_min, tdiff_max, by = bin_width)
  
  df <- df %>%
    mutate(TDiff_bin = cut(transpiration_deficit, breaks = bin_breaks, include.lowest = TRUE, right = FALSE))
  
  # Compute mean for the chosen value column per species, per month, and per TDiff_bin
  meanNDVI_TDiff_species <- df %>%
    group_by(species, month, TDiff_bin) %>%
    summarise(
      avg_value = mean(.data[[value_column]], na.rm = TRUE),
      count = n(),
      .groups = 'drop'
    ) %>%
    filter(count >= 2000) %>%  # Only keep bins with at least 2000 observations
    mutate(bin_median = sapply(as.character(TDiff_bin), function(bin_label) {
      nums <- as.numeric(strsplit(gsub("\\[|\\]|\\(|\\)", "", bin_label), ",")[[1]])
      mean(nums)
    })) %>%
    select(species, month, TDiff_bin, bin_median, avg_value)
  
  return(meanNDVI_TDiff_species)
}

plot_Quantiles_TDiffbin_month_species <- function(data, bin_width = 5, figure_output = NULL) {
  library(ggplot2)
  library(dplyr)
  library(viridis)  # For color-blind-friendly colors
  
  # Define the order of species and months
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  month_order <- c("April", "May", "June", "July", "August")
  
  # Convert species and month to factors with defined order
  data$species <- factor(data$species, levels = species_order)
  data$month <- factor(data$month, levels = month_order)
  
  # Aggregate the data using NDVI_TDiffbin_month function
  aggregated_data <- NDVI_TDiffbin_month(data, bin_width)
  
  # Convert month to factor for correct ordering in the plot
  aggregated_data$month <- factor(aggregated_data$month, levels = month_order)
  
  # Create the plot using a color-blind-friendly palette
  p <- ggplot(aggregated_data, aes(x = bin_median, y = avg_value, color = month, group = month)) +
    geom_line(linewidth = 1.2) +  # Thicker lines for better visibility
    geom_point(size = 2) +
    facet_wrap(~ species, scales = "free_x") +
    scale_color_viridis_d(option = "C", begin = 0.2, end = 0.8) +  # Color-blind friendly
    labs(
      x = "Transpiration Deficit (Bin Median)",
      y = "Mean Quantiles",
      color = "Month",
      title = "Mean Quantiles vs Transpiration Deficit by Month"
    ) +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(color = "black", size = 14),
      panel.border = element_rect(color = NA, fill = NA, linewidth = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    )
  
  # Save the plot if an output filename is provided
  if (!is.null(figure_output)) {
    ggsave(filename = figure_output, plot = p, width = 10, height = 8, dpi = 300)
  }
  
  return(p)
}

TDiff_PSIbin_month <- function(df, bin_width = 100) {
  library(dplyr)
  
  # Identify the correct column dynamically (here we use "transpiration_deficit")
  value_column <- "transpiration_deficit"
  
  # Define bin breaks dynamically based on soil water potential range
  psi_min <- floor(min(df$soil_water_potential, na.rm = TRUE))
  psi_max <- ceiling(max(df$soil_water_potential, na.rm = TRUE))
  bin_breaks <- seq(psi_min, psi_max, by = bin_width)
  
  df <- df %>%
    mutate(PSI_bin = cut(soil_water_potential, breaks = bin_breaks, include.lowest = TRUE, right = FALSE))
  
  # Compute mean for transpiration deficit per species, per month, and per PSI_bin
  meanTDiff_PSIbin_species <- df %>%
    group_by(species, month, PSI_bin) %>%
    summarise(
      avg_value = mean(.data[[value_column]], na.rm = TRUE),
      count = n(),
      .groups = 'drop'
    ) %>%
    filter(count >= 2000) %>%  # Only keep bins with at least 2000 observations
    rowwise() %>%
    mutate(bin_median = ifelse(!is.na(PSI_bin),
                               mean(as.numeric(unlist(strsplit(gsub("\\[|\\]|\\(|\\)", "", as.character(PSI_bin)), ",")))),
                               NA_real_)) %>% # Ensure numeric conversion
    ungroup() %>%
    filter(!is.na(bin_median)) %>%  # Remove rows where bin_median couldn't be calculated
    select(species, month, PSI_bin, bin_median, avg_value)
  
  return(meanTDiff_PSIbin_species)
}

plot_TDiff_PSIbin_month_species <- function(data, bin_width = 100, figure_output = NULL) {
  library(ggplot2)
  library(dplyr)
  library(viridis)  # For color-blind-friendly colors
  
  # Define the order of species and months
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  month_order <- c("April", "May", "June", "July", "August")
  
  # Convert species and month to factors with defined order
  data$species <- factor(data$species, levels = species_order)
  data$month <- factor(data$month, levels = month_order)
  
  # Aggregate the data using TDiff_PSIbin_month function
  aggregated_data <- TDiff_PSIbin_month(data, bin_width)
  
  # Convert month to factor for correct ordering in the plot
  aggregated_data$month <- factor(aggregated_data$month, levels = month_order)
  
  # Ensure bin_median is numeric and exclude missing values
  aggregated_data <- aggregated_data %>%
    filter(!is.na(bin_median)) %>%
    mutate(bin_median = as.numeric(bin_median))
  
  # Create the plot using a color-blind-friendly palette
  p <- ggplot(aggregated_data, aes(x = bin_median, y = avg_value, color = month, group = interaction(month, species))) +
    geom_line(linewidth = 1.2) +  # Thicker lines for better visibility
    geom_point(size = 2) +
    facet_wrap(~ species, scales = "free_x") +
    scale_color_viridis_d(option = "C", begin = 0.2, end = 0.8) +  # Color-blind friendly
    labs(
      x = "Soil Water Potential (Bin Median)",
      y = "Mean Transpiration Deficit",
      color = "Month",
      title = "Mean Transpiration Deficit vs Soil Water Potential by Month"
    ) +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(color = "black", size = 14),
      panel.border = element_rect(color = NA, fill = NA, linewidth = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    )
  
  # Save the plot if an output filename is provided
  if (!is.null(figure_output)) {
    ggsave(filename = figure_output, plot = p, width = 10, height = 8, dpi = 300)
  }
  
  return(p)
}

plot_Quantiles_PSIbin_month_species_linear_agg <- function(data, bin_width = 100,
                                                           figure_output = NULL,
                                                           figure_output2 = NULL) {
  # Load required libraries
  library(ggplot2)
  library(dplyr)
  library(viridis)  # For color-blind-friendly colors
  library(tidyr)
  
  # Define the order of species and months
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  month_order <- c("April", "May", "June", "July", "August")
  
  # Convert species and month in the raw data to factors with the defined order
  data$species <- factor(data$species, levels = species_order)
  data$month <- factor(data$month, levels = month_order)
  
  # Aggregate the data using the helper function (returns bin_median and avg_value)
  aggregated_data <- NDVI_PSIbin_month(data, bin_width)
  
  # Ensure correct month ordering in the aggregated data
  aggregated_data$month <- factor(aggregated_data$month, levels = month_order)
  
  ## Figure 1: Scatter plot with fitted dashed linear regression lines (using aggregated data)
  
  p1 <- ggplot(aggregated_data, aes(x = bin_median, y = avg_value, color = month)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed", aes(group = month)) +
    facet_wrap(~ species, scales = "free_x") +
    scale_color_viridis_d(option = "C", begin = 0.2, end = 0.8) +
    labs(
      x = "Soil Water Potential (Bin Median)",
      y = "Mean Quantiles",
      color = "Month",
      title = "Mean Quantiles vs Soil Water Potential\nwith Fitted Linear Regression"
    ) +
    theme_minimal() +
    theme(
      plot.background    = element_rect(fill = "white", color = "white"),
      panel.background   = element_rect(fill = "white"),
      legend.background  = element_rect(fill = "white", color = "white"),
      plot.title         = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title         = element_text(face = "bold"),
      axis.text          = element_text(color = "black"),
      strip.background   = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text         = element_text(face = "bold", size = 12),
      legend.position    = "top"
    )
  
  ## Figure 2: Bar plot of regression parameters (intercept and slope) for each species and month
  
  # Compute regression parameters for each species and month by filtering the aggregated data
  params <- aggregated_data %>%
    group_by(species, month) %>%
    group_modify(~ {
      # .x is the data for each species-month group (filtered automatically)
      mod <- lm(avg_value ~ bin_median, data = .x)
      tibble(intercept = coef(mod)[1], slope = coef(mod)[2])
    }) %>% 
    ungroup()
  
  # Pivot the parameters into a long format so both intercept and slope can be plotted
  params_long <- params %>%
    pivot_longer(cols = c(intercept, slope), names_to = "parameter", values_to = "value")
  
  p2 <- ggplot(params_long, aes(x = month, y = value, fill = month)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    facet_grid(parameter ~ species, scales = "free_y") +
    scale_fill_viridis_d(option = "C", begin = 0.2, end = 0.8) +
    labs(
      x = "Month",
      y = "Parameter Value",
      fill = "Month",
      title = "Regression Parameters for each Species and Month\n(Y = a + bx)"
    ) +
    theme_minimal() +
    theme(
      plot.background    = element_rect(fill = "white", color = "white"),
      panel.background   = element_rect(fill = "white"),
      legend.background  = element_rect(fill = "white", color = "white"),
      plot.title         = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title         = element_text(face = "bold"),
      axis.text          = element_text(color = "black"),
      strip.background   = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text         = element_text(face = "bold", size = 12),
      legend.position    = "top"
    )
  
  # Save the plots if output filenames are provided
  if (!is.null(figure_output)) {
    ggsave(filename = figure_output, plot = p1, width = 10, height = 8)
  }
  if (!is.null(figure_output2)) {
    ggsave(filename = figure_output2, plot = p2, width = 10, height = 8)
  }
  
  # Return the plots as a list
  return(list(linear_fit_plot = p1, regression_parameters_barplot = p2))
}

plot_Quantiles_PSIbin_month_species_linear_orig <- function(data,
                                                            figure_output = NULL,
                                                            figure_output2 = NULL) {
  # Load required libraries
  library(ggplot2)
  library(dplyr)
  library(viridis)  # For color-blind-friendly colors
  library(tidyr)
  
  # Define the order of species and months
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  month_order   <- c("April", "May", "June", "July", "August")
  
  # Convert species and month in the raw data to factors with the defined order
  data$species <- factor(data$species, levels = species_order)
  data$month   <- factor(data$month, levels = month_order)
  
  ## Figure 1: Scatter plot with fitted dashed linear regression lines using original data
  
  p1 <- ggplot(data, aes(x = soil_water_potential, y = Quantiles, color = month)) +
    # geom_point() +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed", aes(group = month)) +
    facet_wrap(~ species, scales = "free_x") +
    scale_color_viridis_d(option = "C", begin = 0.2, end = 0.8) +
    labs(
      x = "Soil Water Potential",
      y = "Quantiles",
      color = "Month",
      title = "Quantiles vs Soil Water Potential\nwith Fitted Linear Regression"
    ) +
    theme_minimal() +
    theme(
      plot.background   = element_rect(fill = "white", color = "white"),
      panel.background  = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title        = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title        = element_text(face = "bold"),
      axis.text         = element_text(color = "black"),
      strip.background  = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text        = element_text(face = "bold", size = 12),
      legend.position   = "top"
    )
  
  ## Figure 2: Bar plot of regression parameters (intercept and slope) for each species and month
  
  # Compute regression parameters for each species and month from the original data
  params <- data %>%
    group_by(species, month) %>%
    group_modify(~ {
      mod <- lm(Quantiles ~ soil_water_potential, data = .x)
      tibble(intercept = coef(mod)[1], slope = coef(mod)[2])
    }) %>%
    ungroup()
  
  # Pivot the parameters into a long format so both intercept and slope can be plotted
  params_long <- params %>%
    pivot_longer(cols = c(intercept, slope), names_to = "parameter", values_to = "value")
  
  p2 <- ggplot(params_long, aes(x = month, y = value, fill = month)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    facet_grid(parameter ~ species, scales = "free_y") +
    scale_fill_viridis_d(option = "C", begin = 0.2, end = 0.8) +
    labs(
      x = "Month",
      y = "Parameter Value",
      fill = "Month",
      title = "Regression Parameters for each Species and Month\n(Y = a + bx)"
    ) +
    theme_minimal() +
    theme(
      plot.background   = element_rect(fill = "white", color = "white"),
      panel.background  = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title        = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title        = element_text(face = "bold"),
      axis.text         = element_text(color = "black"),
      strip.background  = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text        = element_text(face = "bold", size = 12),
      legend.position   = "top"
    )
  
  # Save the plots if output filenames are provided
  if (!is.null(figure_output)) {
    ggsave(filename = figure_output, plot = p1, width = 10, height = 8)
  }
  if (!is.null(figure_output2)) {
    ggsave(filename = figure_output2, plot = p2, width = 10, height = 8)
  }
  
  # Return the plots as a list
  return(list(linear_fit_plot = p1, regression_parameters_barplot = p2))
}

plot_TDiff_PSIbin_month_species_original <- function(data, bin_width = 100, figure_output = NULL) {
  library(ggplot2)
  library(dplyr)
  library(viridis)  # For color-blind-friendly colors
  
  # Define the order of species and months
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  month_order <- c("April", "May", "June", "July", "August")
  
  # Convert species and month to factors with defined order
  data$species <- factor(data$species, levels = species_order)
  data$month <- factor(data$month, levels = month_order)
  
  # Calculate binned means using your function (assumes TDiff_PSIbin_month is defined)
  summary_data <- TDiff_PSIbin_month(data, bin_width = bin_width)
  
  # Create the plot:
  # - Plot original data points with geom_point
  # - Overlay a line using the mean values from each PSI bin for each species and month
  p <- ggplot() +
    geom_point(data = data, 
               aes(x = soil_water_potential, y = transpiration_deficit, color = month),
               size = 2, alpha = 0.6) +
    geom_line(data = summary_data, 
              aes(x = bin_median, y = avg_value, color = month, group = interaction(species, month)),
              linewidth = 1.2) +
    facet_wrap(~ species, scales = "free_x") +
    scale_color_viridis_d(option = "C", begin = 0.2, end = 0.8) +
    labs(
      x = "Soil Water Potential",
      y = "Transpiration Deficit",
      color = "Month",
      title = "Transpiration Deficit vs Soil Water Potential by Month and Species"
    ) +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(color = "black", size = 14),
      panel.border = element_rect(color = NA, fill = NA, linewidth = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    )
  
  # Save the plot if an output filename is provided
  if (!is.null(figure_output)) {
    ggsave(filename = figure_output, plot = p, width = 10, height = 8, dpi = 300)
  }
  
  return(p)
}

