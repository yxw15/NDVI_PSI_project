setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

load("results/Data/All_Species_Quantiles_PSI_TDiff.RData")

plot_species_distribution <- function(data, output_path) {
  # Load required libraries
  library(dplyr)
  library(ggplot2)
  library(maps)
  
  # Define species order and color palette
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  cb_palette <- c("Oak"   = "#E69F00",  # Orange
                  "Beech" = "#0072B2",  # Deep blue
                  "Spruce"= "#009E73",  # Bluish-green
                  "Pine"  = "#F0E442")  # Yellow
  
  # Extract unique x, y coordinates for each species
  unique_points <- data %>% 
    distinct(x, y, species)
  
  # Ensure species factor levels are in the desired order
  unique_points$species <- factor(unique_points$species, levels = species_order)
  
  # Load simplified map data for Germany
  germany_map <- map_data("world", region = "Germany")
  
  # Create the plot
  p <- ggplot() +
    # Draw the map
    geom_polygon(data = germany_map, 
                 aes(x = long, y = lat, group = group),
                 fill = "white", 
                 color = "black", 
                 linewidth = 0.5) +
    # Add species points
    geom_point(data = unique_points, 
               aes(x = x, y = y, color = species), 
               size = 1) +
    # Apply color palette
    scale_color_manual(values = cb_palette, name = "") +
    # Set horizontal facet layout (1 row, 4 columns)
    facet_wrap(~ species, ncol = 4) +
    # Use coord_quickmap() to accurately project map data 
    coord_quickmap() +
    # Plot title and axis labels
    labs(title = "Distribution of Tree Species in Germany",
         x = "longitude",
         y = "latitude") +
    # Bigger legend points
    guides(color = guide_legend(override.aes = list(size = 4))) +
    # Theme settings
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
  
  # Save and return the plot
  ggsave(filename = output_path, plot = p, width = 14, height = 6, dpi = 300)
  print(p)
  return(p)
}

plot_moran_I_distance <- function(input_path, output_figure) {
  # Load necessary libraries
  library(ggplot2)
  library(dplyr)
  
  # Define your custom color palette and species ordering
  cb_palette <- c("Oak"   = "#E69F00",    # Orange
                  "Beech" = "#0072B2",    # Deep blue
                  "Spruce"= "#009E73",    # Bluish-green
                  "Pine"  = "#F0E442")    # Yellow
  species_levels <- c("Oak", "Beech", "Spruce", "Pine")
  
  # Custom theme function
  create_custom_theme <- function() {
    theme_minimal() +
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
  }
  
  # Get the list of CSV files from the folder
  files <- list.files(input_path, full.names = TRUE)
  
  # Helper function to read file and add a "Parameter" column based on filename
  read_moran_data <- function(file) {
    df <- read.csv(file)
    
    # Determine parameter based on file name contents:
    parameter <- if (grepl("Quantiles", file)) {
      "Quantiles"
    } else if (grepl("soil_water_potential", file)) {
      "Soil Water Potential"
    } else if (grepl("transpiration_deficit", file)) {
      "Transpiration Deficit"
    } else {
      NA
    }
    
    df$Parameter <- parameter
    return(df)
  }
  
  # Combine all CSV files into one data frame
  all_data <- do.call(rbind, lapply(files, read_moran_data))
  
  # Factorize Species and recode Parameter for desired labels in the plots
  all_data$Species <- factor(all_data$Species, levels = species_levels)
  all_data$Parameter <- factor(all_data$Parameter, 
                               levels = c("Quantiles", "Soil Water Potential", "Transpiration Deficit"),
                               labels = c("NDVI quantiles", "soil water potential", "transpiration deficit"))
  
  # Create the plot
  p <- ggplot(all_data, aes(x = Distance, y = Moran_I, color = Species)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    # Add a geom_text layer: if Significance is "*" then display "*" at the corresponding point
    geom_text(aes(label = ifelse(Significance == "*", "*", "")), 
              vjust = 2, 
              size = 6,
              show.legend = FALSE) +
    facet_wrap(~ Parameter, scales = "free_y") +
    scale_color_manual(values = cb_palette, name = "") +
    labs(title = "Moran's I across different parameters and species",
         x = "distance (m)",
         y = "Moran's I") +
    create_custom_theme()
  
  # Print the plot to the current device
  print(p)
  
  # Save the plot to the specified output file
  ggsave(filename = output_figure, plot = p, width = 10, height = 6)
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
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(color = "black", size = 14),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  # Create an x-axis scale that shows breaks every 2 years from 2003 to 2024
  x_scale <- scale_x_continuous(breaks = seq(2003, 2024, by = 2), limits = c(2003, 2024))
  
  ### TIME SERIES PLOTS ###
  p1 <- ggplot(df_species, aes(x = as.numeric(year), y = .data[[ndvi_column]], color = species)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_line() +
    labs(title = "", x = "", y = expression(atop(bold("NDVI quantiles"), bold("(rank)"))), color = "Species") +
    scale_color_manual(values = cb_palette, name = "", guide = guide_legend(nrow = 1)) +
    custom_theme +
    theme(legend.position = "bottom", legend.direction = "horizontal") +
    x_scale +
    annotate("text", x = -Inf, y = Inf, label = "(a)", hjust = -0.2, vjust = 1.5, size = 6, fontface = "bold")
  
  p2 <- ggplot(df_species, aes(x = as.numeric(year), y = avg_psi, color = species)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_line() +
    labs(title = "", x = "", y = expression(atop(bold("soil water potential"), bold("(kPa)"))), color = "Species") +
    scale_color_manual(values = cb_palette, name = "") +
    custom_theme +
    guides(color = "none") +
    x_scale +
    annotate("text", x = -Inf, y = Inf, label = "(b)", hjust = -0.2, vjust = 1.5, size = 6, fontface = "bold")
  
  p3 <- ggplot(df_species, aes(x = as.numeric(year), y = avg_tdiff, color = species)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_line() +
    labs(title = "", x = "", y = expression(atop(bold("transpiration deficit"), bold("(mm)"))), color = "Species") +
    scale_color_manual(values = cb_palette, name = "") +
    custom_theme +
    guides(color = "none") +
    x_scale +
    annotate("text", x = -Inf, y = Inf, label = "(c)", hjust = -0.2, vjust = 1.5, size = 6, fontface = "bold")
  
  ts_plot <- p1 / p2 / p3
  
  ### CORRELATION PLOTS ###
  annotations_d <- df_species %>%
    group_by(species) %>%
    summarize(corr = round(cor(avg_psi, .data[[ndvi_column]], use = "complete.obs"), 2))
  
  print(annotations_d)
  
  annotations_e <- df_species %>%
    group_by(species) %>%
    summarize(corr = round(cor(avg_tdiff, .data[[ndvi_column]], use = "complete.obs"), 2))
  print(annotations_e)
  
  # text_d <- paste0(species_order, ": r = ", 
  #                  annotations_d$corr[match(species_order, annotations_d$species)],
  #                  collapse = "\n")
  # 
  # text_e <- paste0(species_order, ": r = ", 
  #                  annotations_e$corr[match(species_order, annotations_e$species)],
  #                  collapse = "\n")
  
  p_d <- ggplot(df_species, aes(x = avg_psi, y = .data[[ndvi_column]], color = species)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(x = "soil water potential (kPa)", y = expression(atop(bold("NDVI quantiles"), bold("(rank)"))), color = "Species") +
    scale_color_manual(values = cb_palette, name = "") +
    custom_theme +
    guides(color = "none") +
    annotate("text", x = -Inf, y = Inf, label = "(d)", hjust = -0.2, vjust = 1.5, size = 6, fontface = "bold") 
    # annotate("text", x = -Inf, y = Inf, label = text_d, hjust = 0, vjust = 4.5, size = 5, color = "black")
  
  p_e <- ggplot(df_species, aes(x = avg_tdiff, y = .data[[ndvi_column]], color = species)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(x = "transpiration deficit (mm)", y = "", color = "Species") +
    scale_color_manual(values = cb_palette, name = "") +
    custom_theme +
    guides(color = "none") +
    annotate("text", x = Inf, y = Inf, label = "(e)", hjust = 1.1, vjust = 1.5, size = 6, fontface = "bold") 
    # annotate("text", x = Inf, y = Inf, label = text_e, hjust = 1, vjust = 1.5, size = 5, color = "black")
  
  corr_plot <- p_d | p_e
  
  combined_main <- ts_plot / corr_plot
  combined_main <- combined_main & theme(legend.position = "top")
  
  legend <- suppressWarnings(cowplot::get_legend(p1))
  
  final_plot <- cowplot::plot_grid(combined_main, legend, ncol = 1, rel_heights = c(1, 0.1))
  
  print(final_plot)
  ggsave(filename = output_path, plot = final_plot, width = 12, height = 14, dpi = 300)
}

plot_yearly_mean_linear_coeffs <- function(df_all, output_path) {
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  
  # Get species-level averaged data (assuming df_average_year_species is available)
  df_species <- df_average_year_species(df_all)
  
  # Set species order and assign factor levels
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  df_species$species <- factor(df_species$species, levels = species_order)
  
  # Determine NDVI column and label
  ndvi_column <- if ("avg_quantile" %in% names(df_species)) {
    "avg_quantile"
  } else if ("avg_proportion" %in% names(df_species)) {
    "avg_proportion"
  } else {
    stop("Neither avg_quantile nor avg_proportion column found in species data.")
  }
  
  ndvi_label <- if (ndvi_column == "avg_quantile") "NDVI quantiles" else "NDVI proportions"
  
  # Define the color palette for the species
  cb_palette <- c("Oak"   = "#E69F00", 
                  "Beech" = "#0072B2",
                  "Spruce"= "#009E73", 
                  "Pine"  = "#F0E442")
  
  # Define a custom theme similar to your provided theme.
  custom_theme <- theme(
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
  
  ### NDVI ~ soil water potential Regression Coefficients (Panel a) ###
  coeffs_psi <- do.call(rbind, lapply(species_order, function(sp) {
    df_sp <- df_species %>% filter(species == sp)
    mod <- lm(as.formula(paste0(ndvi_column, " ~ avg_psi")), data = df_sp)
    coef_vals <- coef(mod)
    p_vals <- summary(mod)$coefficients[,4]
    data.frame(species = sp, 
               term = names(coef_vals), 
               estimate = coef_vals, 
               p.value = p_vals,
               stringsAsFactors = FALSE)
  }))
  # Rename coefficient terms: intercept -> a, slope -> b
  coeffs_psi$term <- dplyr::recode(coeffs_psi$term, "(Intercept)" = "a", "avg_psi" = "b")
  # Ensure species factor order is maintained
  coeffs_psi$species <- factor(coeffs_psi$species, levels = species_order)
  coeffs_psi$label <- ifelse(coeffs_psi$p.value < 0.05, "*", sprintf("%.2f", coeffs_psi$p.value))
  
  p_coeff_psi <- ggplot(coeffs_psi, aes(x = species, y = estimate, fill = species)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = label, y = estimate/2),
              position = position_dodge(width = 0.9), color = "black", size = 4) +
    scale_fill_manual(values = cb_palette, name = "") +
    facet_wrap(~ term, scales = "free_y") +
    labs(title = "NDVI quantiles ~ soil water potential",
         subtitle = expression(NDVI == a + b * x + epsilon),
         x = "", y = "coefficient value") +
    custom_theme +
    labs(caption = "(a)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0.5))
  
  ### NDVI ~ transpiration deficit Regression Coefficients (Panel b) ###
  coeffs_tdiff <- do.call(rbind, lapply(species_order, function(sp) {
    df_sp <- df_species %>% filter(species == sp)
    mod <- lm(as.formula(paste0(ndvi_column, " ~ avg_tdiff")), data = df_sp)
    coef_vals <- coef(mod)
    p_vals <- summary(mod)$coefficients[,4]
    data.frame(species = sp, 
               term = names(coef_vals), 
               estimate = coef_vals, 
               p.value = p_vals,
               stringsAsFactors = FALSE)
  }))
  # Rename coefficient terms: intercept -> a, slope -> b
  coeffs_tdiff$term <- dplyr::recode(coeffs_tdiff$term, "(Intercept)" = "a", "avg_tdiff" = "b")
  # Ensure species factor order is maintained
  coeffs_tdiff$species <- factor(coeffs_tdiff$species, levels = species_order)
  coeffs_tdiff$label <- ifelse(coeffs_tdiff$p.value < 0.05, "*", sprintf("%.2f", coeffs_tdiff$p.value))
  
  p_coeff_tdiff <- ggplot(coeffs_tdiff, aes(x = species, y = estimate, fill = species)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = label, y = estimate/2),
              position = position_dodge(width = 0.9), color = "black", size = 4) +
    scale_fill_manual(values = cb_palette, name = "") +
    facet_wrap(~ term, scales = "free_y") +
    labs(title = "NDVI quantile ~ transpiration deficit",
         subtitle = expression(NDVI == a + b * x + epsilon),
         x = "", y = "") +
    custom_theme +
    labs(caption = "(b)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0.5))
  
  ### Combine the Two Panels Side-by-Side with a Shared Legend and Overall Caption ###
  combined_plot <- p_coeff_psi + p_coeff_tdiff + 
    plot_layout(guides = "collect") &
    theme(legend.position = "top")
  
  print(combined_plot)
  ggsave(filename = output_path, plot = combined_plot, width = 14, height = 8, dpi = 300)
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

NDVI_TDiffbin <- function(df, bin_width = 3) {
  
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
    mutate(percentage = count / total_pixels) %>%
    filter(percentage >= 0.0001) %>%
    select(species, TDiff_bin, bin_median, avg_value, count, total_pixels, percentage)
  
  return(meanNDVI_TDiffbin_species)
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

plot_combined_AIC_R2 <- function(data, save_combined_fig) {
  
  # Load required libraries
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(patchwork)
  
  # Define common species order
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  
  #### Panel A: NDVI ~ PSIbin ####
  data_a <- NDVI_PSIbin(data)
  data_a <- na.omit(data_a)
  value_col_a <- "avg_value"
  data_a$species <- factor(data_a$species, levels = species_order)
  data_a <- data_a %>% mutate(x = -bin_median)
  
  start_list_a <- list(a = 5, b = 3, c = 0.001)
  control_params_a <- nls.control(maxiter = 200, minFactor = 1e-4)
  
  aic_results_a <- list()
  for (sp in levels(data_a$species)) {
    sp_data <- data_a %>% filter(species == sp)
    
    # Linear model
    lm_linear <- lm(avg_value ~ x, data = sp_data)
    aic_linear <- AIC(lm_linear)
    r2_linear <- summary(lm_linear)$r.squared
    
    # Exponential model
    aic_exp <- NA
    r2_exp <- NA
    nls_exp <- tryCatch({
      nls(avg_value ~ a + b * exp(-c * x),
          data = sp_data,
          start = start_list_a,
          control = control_params_a)
    }, error = function(e) NULL)
    if (!is.null(nls_exp)) {
      aic_exp <- AIC(nls_exp)
      res <- resid(nls_exp)
      ss_res <- sum(res^2)
      ss_tot <- sum((sp_data[[value_col_a]] - mean(sp_data[[value_col_a]]))^2)
      r2_exp <- 1 - ss_res / ss_tot
    }
    
    sp_results <- data.frame(
      species = sp,
      Model = c("linear", "exponential"),
      AIC = c(aic_linear, aic_exp),
      R2 = c(r2_linear, r2_exp)
    )
    sp_results$y_label_pos <- sp_results$AIC / 2
    aic_results_a[[sp]] <- sp_results
  }
  aic_df_a <- do.call(rbind, aic_results_a)
  aic_df_a$species <- factor(aic_df_a$species, levels = species_order)
  
  # Shared color palette for Panels A and B
  model_palette_shared <- c("linear" = "orange",
                            "exponential" = "dodgerblue")
  
  p_a <- ggplot(aic_df_a, aes(x = species, y = AIC, fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(aes(label = round(R2, 2), y = y_label_pos),
              position = position_dodge(width = 0.9),
              vjust = 0.5, size = 8, color = "black") +
    labs(x = "", y = "AIC", title = "NDVI quantiles ~ soil water potential") +
    scale_fill_manual(values = model_palette_shared, name = "") +
    theme_minimal() +
    theme(
      axis.text.x    = element_text(angle = 0, hjust = 0.5, size = 14),
      plot.background= element_rect(fill = "white", color = "white"),
      panel.background=element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title     = element_text(hjust = 0.5, size = 22, face = "bold", color = "black"),
      axis.title     = element_text(face = "bold", size = 16),
      axis.text      = element_text(color = "black", size = 14),
      panel.border   = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.text    = element_text(size = 14)
    )+
    labs(caption = "(a)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0.5))
  
  #### Panel B: NDVI ~ TDiffbin ####
  data_b <- NDVI_TDiffbin(data)
  data_b <- na.omit(data_b)
  value_col_b <- "avg_value"
  data_b$species <- factor(data_b$species, levels = species_order)
  data_b <- data_b %>% mutate(x = bin_median)
  
  start_list_b <- list(a = 5, b = 7, c = 0.04)
  control_params_b <- nls.control(maxiter = 1200, minFactor = 1e-09)
  
  aic_results_b <- list()
  for (sp in levels(data_b$species)) {
    sp_data <- data_b %>% filter(species == sp)
    
    lm_linear <- lm(avg_value ~ x, data = sp_data)
    aic_linear <- AIC(lm_linear)
    r2_linear <- summary(lm_linear)$r.squared
    
    aic_exp <- NA
    r2_exp <- NA
    nls_exp <- tryCatch({
      nls(avg_value ~ a + b * exp(-c * x),
          data = sp_data,
          start = start_list_b,
          control = control_params_b)
    }, error = function(e) NULL)
    if (!is.null(nls_exp)) {
      aic_exp <- AIC(nls_exp)
      res <- resid(nls_exp)
      ss_res <- sum(res^2)
      ss_tot <- sum((sp_data[[value_col_b]] - mean(sp_data[[value_col_b]]))^2)
      r2_exp <- 1 - ss_res / ss_tot
    }
    
    sp_results <- data.frame(
      species = sp,
      Model = c("linear", "exponential"),
      AIC = c(aic_linear, aic_exp),
      R2 = c(r2_linear, r2_exp)
    )
    sp_results$y_label_pos <- sp_results$AIC / 2
    aic_results_b[[sp]] <- sp_results
  }
  aic_df_b <- do.call(rbind, aic_results_b)
  aic_df_b$species <- factor(aic_df_b$species, levels = species_order)
  
  p_b <- ggplot(aic_df_b, aes(x = species, y = AIC, fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(aes(label = round(R2, 2), y = y_label_pos),
              position = position_dodge(width = 0.9),
              vjust = 0.5, size = 8, color = "black") +
    labs(x = "", y = "AIC", title = "NDVI quantiles ~ transpiration deficit") +
    scale_fill_manual(values = model_palette_shared, name = "") +
    theme_minimal() +
    theme(
      axis.text.x    = element_text(angle = 0, hjust = 0.5, size = 14),
      plot.background= element_rect(fill = "white", color = "white"),
      panel.background=element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title     = element_text(hjust = 0.5, size = 22, face = "bold", color = "black"),
      axis.title     = element_text(face = "bold", size = 16),
      axis.text      = element_text(color = "black", size = 14),
      panel.border   = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.text    = element_text(size = 14)
    )+
    labs(caption = "(b)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0.5))
  
  #### Panel C: TDiff ~ PSIbin ####
  data_c <- TDiff_PSIbin(data)
  data_c <- na.omit(data_c)
  value_col_c <- "avg_transpiration_deficit"
  data_c$species <- factor(data_c$species, levels = species_order)
  data_c <- data_c %>% mutate(x = bin_median)
  
  aic_results_c <- list()
  for (sp in levels(data_c$species)) {
    sp_data <- data_c %>% filter(species == sp)
    
    lm_linear <- lm(avg_transpiration_deficit ~ x, data = sp_data)
    aic_linear <- AIC(lm_linear)
    r2_linear <- summary(lm_linear)$r.squared
    
    lm_poly2 <- lm(avg_transpiration_deficit ~ x + I(x^2), data = sp_data)
    aic_poly2 <- AIC(lm_poly2)
    r2_poly2 <- summary(lm_poly2)$r.squared
    
    lm_poly3 <- lm(avg_transpiration_deficit ~ x + I(x^2) + I(x^3), data = sp_data)
    aic_poly3 <- AIC(lm_poly3)
    r2_poly3 <- summary(lm_poly3)$r.squared
    
    sp_results <- data.frame(
      species = sp,
      Model = c("linear", "poly2", "poly3"),
      AIC = c(aic_linear, aic_poly2, aic_poly3),
      R2 = c(r2_linear, r2_poly2, r2_poly3)
    )
    sp_results$y_label_pos <- sp_results$AIC / 2
    aic_results_c[[sp]] <- sp_results
  }
  aic_df_c <- do.call(rbind, aic_results_c)
  aic_df_c$species <- factor(aic_df_c$species, levels = species_order)
  
  model_palette_c <- c("linear" = "orange",
                       "poly2"  = "dodgerblue",
                       "poly3"  = "green4")
  
  p_c <- ggplot(aic_df_c, aes(x = species, y = AIC, fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(aes(label = round(R2, 2), y = y_label_pos),
              position = position_dodge(width = 0.9),
              vjust = 0.5, size = 8, color = "black") +
    labs(x = "", y = "AIC", title = "transpiration deficit ~ soil water potential") +
    scale_fill_manual(values = model_palette_c, name = "") +
    theme_minimal() +
    theme(
      axis.text.x    = element_text(angle = 0, hjust = 0.5, size = 14),
      plot.background= element_rect(fill = "white", color = "white"),
      panel.background=element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title     = element_text(hjust = 0.5, size = 22, face = "bold", color = "black"),
      axis.title     = element_text(face = "bold", size = 16),
      axis.text      = element_text(color = "black", size = 14),
      panel.border   = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.text    = element_text(size = 14)
    )+
    labs(caption = "(c)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0.5))
  
  # Combine Panels A and B with a shared legend placed below (centered)
  combined_top <- (p_a | p_b) + 
    plot_layout(guides = "collect") & 
    theme(legend.position = "top", legend.box = "horizontal")
  
  # Combine with Panel C below
  combined_plot <- combined_top / p_c
  
  # Save the combined figure to file
  dir.create(dirname(save_combined_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = save_combined_fig, plot = combined_plot, width = 14, height = 12, dpi = 300)
  
  print(combined_plot)
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
  cb_palette <- c("Oak"   = "#E69F00",
                  "Beech" = "#0072B2",
                  "Spruce"= "#009E73",
                  "Pine"  = "#F0E442")
  
  data <- data %>% mutate(x = bin_median)
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  threshold <- 11.5
  
  #### NONLINEAR MODELING PER SPECIES ####
  start_list <- list(a = 5, b = 3, c = 0.001)
  control_params <- nls.control(maxiter = 200, minFactor = 1e-4)
  nls_models <- nlsList(avg_value ~ a + b * exp(c * x) | species,
                        data = data_clean,
                        start = start_list,
                        control = control_params)
  print(summary(nls_models))
  
  coef_df <- as.data.frame(coef(nls_models), optional = TRUE) %>% 
    rownames_to_column(var = "species") %>%
    filter(!is.na(a))
  coef_df$species <- factor(coef_df$species, levels = species_order)
  
  # Compute x50, slope50, and plateau point
  coef_df <- coef_df %>%
    mutate(x50 = ifelse((threshold - a) > 0 & b > 0,
                        log((threshold - a)/b) / c,
                        NA),
           slope50 = c * (threshold - a),
           x_plateau = log(0.05) / c)  # new line for plateau point
  
  # Print plateau values for all species
  plateau_df <- coef_df %>%
    select(species, x_plateau)
  
  print("Estimated soil water potential (x_plateau) where NDVI is ~95% of asymptotic minimum (plateau):")
  print(plateau_df)
  
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
  
  x_scale <- scale_x_continuous()
  p_combined <- ggplot() +
    geom_point(data = data_clean, aes(x = x, y = avg_value, color = species)) +
    geom_line(data = pred_all, aes(x = x, y = pred, color = species), linewidth = 1) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "black", linewidth = 1) +
    annotate("text", x = min(data_clean$x), y = threshold, label = "median", 
             hjust = -0.1, vjust = -0.3, fontface = "italic", size = 6) +
    scale_color_manual(values = cb_palette, name = "") +
    x_scale +
    labs(x = "soil water potential (kPa)", y = "NDVI quantiles (rank)") +
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
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0.5))+
    coord_cartesian(clip = "off")
  
  #### PANEL B: Bar Plot of x50 Values ####
  p_x50 <- ggplot(coef_df, aes(x = species, y = x50, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    labs(x = "", y = "soil water potential (kPa)") +
    scale_fill_manual(values = cb_palette) +
    expand_limits(y = 0) +
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
  
  #### PANEL C: Absolute Slope at x50 ####
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
  
  print(stats_df)
  
  p_slope <- ggplot(stats_df, aes(x = species, y = slope_abs, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, slope_abs - se), ymax = slope_abs + se), width = 0.2) +
    geom_text(aes(y = slope_abs/2, label = label_text), color = "black", size = 5) +
    labs(x = "", y = "absolute slope") +
    scale_fill_manual(values = cb_palette) +
    expand_limits(y = 0) +
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
  
  #### Combine Panels and Save Plot ####
  final_slope_plot <- (p_combined + (p_x50 / p_slope)) +
    plot_layout(widths = c(2, 1), guides = "collect") &
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 14),
      legend.key = element_rect(fill = "white", color = NA),
      legend.background = element_blank(),
      legend.box.background = element_blank()
    )
  
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
    mutate(label = ifelse(`Pr(>|t|)` < 0.05, "*", sprintf("%.2f", `Pr(>|t|)`)))
  
  coeffs_long <- coef_df %>%
    select(species, a, b, c) %>%
    pivot_longer(cols = c("a", "b", "c"),
                 names_to = "Coefficient",
                 values_to = "Value")
  
  coeffs_long <- left_join(coeffs_long, coeff_stats %>% select(species, Coefficient, label),
                           by = c("species", "Coefficient"))
  coeffs_long$species <- factor(coeffs_long$species, levels = species_order)
  
  print(coeffs_long)
  
  p_coeffs <- ggplot(coeffs_long, aes(x = species, y = Value, fill = species)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = label, y = Value/2),
              color = "black", size = 4,
              position = position_dodge(width = 0.9)) +
    scale_fill_manual(values = cb_palette, name = "") +
    facet_wrap(~ Coefficient, scales = "free_y") +
    labs(title = "NDVI quantiles ~ soil water potential",
         subtitle = expression(NDVI == a + b * e^{c * italic(x)} + epsilon),
         x = "",
         y = "coefficient value") +
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
  
  print(p_coeffs)
  dir.create(dirname(save_coeff_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = save_coeff_fig, plot = p_coeffs, device = "png", width = 10, height = 8, dpi = 300)
}

plot_Quantiles_TDiff_exp_linear_slope_coeff <- function(data, combined_coef_fig, output_figure) {
  # Process data using the custom NDVI_TDiffbin function and remove missing rows
  data <- NDVI_TDiffbin(data)
  data <- na.omit(data)
  
  # Load required libraries
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(patchwork)
  library(purrr)
  library(nlme)    # for consistency (even though lm is used for linear)
  library(car)     # for deltaMethod
  library(broom)
  
  # Define the value column and species groups.
  # Use linear regression for Oak and Beech,
  # and the exponential model for Spruce and Pine.
  value_col <- "avg_value"
  linear_species <- c("Oak", "Beech")
  exp_species <- c("Spruce", "Pine")
  data$species <- factor(data$species, levels = c(linear_species, exp_species))
  
  # Define the color palette (same for all species)
  cb_palette <- c("Oak"   = "#E69F00",  # Orange
                  "Beech" = "#0072B2",  # Deep blue
                  "Spruce"= "#009E73",  # Bluish-green
                  "Pine"  = "#F0E442")  # Yellow
  
  # For transpiration deficit, set x = bin_median
  data <- data %>% mutate(x = bin_median)
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  
  # Set threshold (using a fixed value as before)
  threshold <- 11.5
  
  ##########################
  # Model Fitting by Group #
  ##########################
  
  # Linear models for Oak and Beech
  models_linear <- list()
  for(sp in linear_species) {
    sp_data <- data_clean %>% filter(species == sp)
    lm_model <- lm(avg_value ~ x, data = sp_data)
    models_linear[[sp]] <- lm_model
  }
  
  # Exponential models for Spruce and Pine
  start_list <- list(a = 5, b = 7, c = 0.04)
  control_params <- nls.control(maxiter = 1200, minFactor = 1e-09)
  models_exp <- list()
  for(sp in exp_species) {
    sp_data <- data_clean %>% filter(species == sp)
    exp_model <- tryCatch({
      nls(avg_value ~ a + b * exp(-c * x), 
          data = sp_data, 
          start = start_list, 
          control = control_params)
    }, error = function(e) {
      NULL
    })
    models_exp[[sp]] <- exp_model
  }
  
  ##############################################
  # Coefficient Extraction and Plotting - Linear
  ##############################################
  
  lin_coef_list <- lapply(linear_species, function(sp) {
    mod <- models_linear[[sp]]
    if(is.null(mod)) return(NULL)
    summ <- summary(mod)$coefficients
    df <- data.frame(
      species = sp,
      Coefficient = rownames(summ),
      Value = summ[, "Estimate"],
      pvalue = summ[, "Pr(>|t|)"],
      stringsAsFactors = FALSE
    )
    # Rename coefficients for clarity using dplyr::recode
    df$Coefficient <- dplyr::recode(df$Coefficient, "(Intercept)" = "a", "x" = "b")
    return(df)
  })
  lin_coef_df <- do.call(rbind, lin_coef_list)
  lin_coef_df$species <- factor(lin_coef_df$species, levels = linear_species)
  lin_coef_df <- lin_coef_df %>% 
    mutate(label = if_else(pvalue < 0.05, "*", sprintf("%.2f", pvalue)))
  
  p_coeff_linear <- ggplot(lin_coef_df, aes(x = species, y = Value, fill = species)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(aes(label = label, y = Value/2),
              position = position_dodge(width = 0.9),
              vjust = 0.5, color = "black", size = 4) +
    scale_fill_manual(values = cb_palette, name = "") +
    facet_wrap(~ Coefficient, scales = "free_y") +
    labs(title = "linear models",
         # subtitle = expression(italic(NDVI[Quantiles]) == a + b * x + epsilon),
         subtitle = expression(NDVI == a + b * x + epsilon),
         x = "",
         y = "coefficient value") +
    theme_minimal() +
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
    ) +
    labs(caption = "(a)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0.5))
  
  ##############################################
  # Coefficient Extraction and Plotting - Exponential
  ##############################################
  
  exp_coef_list <- lapply(exp_species, function(sp) {
    mod <- models_exp[[sp]]
    if(is.null(mod)) return(NULL)
    summ <- summary(mod)$coefficients
    df <- data.frame(
      species = sp,
      Coefficient = rownames(summ),
      Value = summ[, "Estimate"],
      pvalue = summ[, "Pr(>|t|)"],
      stringsAsFactors = FALSE
    )
    # Convert coefficient names to lowercase (a, b, c)
    df$Coefficient <- tolower(df$Coefficient)
    return(df)
  })
  exp_coef_df <- do.call(rbind, exp_coef_list)
  exp_coef_df$species <- factor(exp_coef_df$species, levels = exp_species)
  exp_coef_df <- exp_coef_df %>% 
    mutate(label = if_else(pvalue < 0.05, "*", sprintf("%.2f", pvalue)))
  
  p_coeff_exp <- ggplot(exp_coef_df, aes(x = species, y = Value, fill = species)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(aes(label = label, y = Value/2),
              position = position_dodge(width = 0.9),
              vjust = 0.5, color = "black", size = 4) +
    scale_fill_manual(values = cb_palette, name = "") +
    facet_wrap(~ Coefficient, scales = "free_y") +
    labs(title = "exponential models",
         # subtitle = expression(italic(NDVI[Quantiles]) == a + b * e^{-c * italic(x)} + epsilon),
         subtitle = expression(NDVI == a + b * e^{-c * italic(x)} + epsilon),
         x = "",
         y = "") +
    theme_minimal() +
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
    ) +
    labs(caption = "(b)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0.5))
  
  # Combine the coefficient plots into one figure with panel labels (a) and (b)
  combined_coeff <- p_coeff_linear + p_coeff_exp 
  print(combined_coeff)
  
  # Save the combined coefficient plot
  ggsave(filename = combined_coef_fig, plot = combined_coeff, device = "png", width = 10, height = 8, dpi = 300)
  
  #################################################
  # Combined Final Plot with Model Predictions (Panel A)
  #################################################
  
  # Generate predictions for linear species (Oak & Beech)
  pred_linear <- lapply(linear_species, function(sp) {
    sp_data <- data_clean %>% filter(species == sp)
    x_seq <- seq(min(sp_data$x, na.rm = TRUE), max(sp_data$x, na.rm = TRUE), length.out = 100)
    lm_mod <- models_linear[[sp]]
    pred <- predict(lm_mod, newdata = data.frame(x = x_seq))
    data.frame(species = sp, x = x_seq, pred = pred)
  })
  pred_linear_df <- do.call(rbind, pred_linear)
  
  # Generate predictions for exponential species (Spruce & Pine)
  pred_exp <- lapply(exp_species, function(sp) {
    sp_data <- data_clean %>% filter(species == sp)
    x_seq <- seq(min(sp_data$x, na.rm = TRUE), max(sp_data$x, na.rm = TRUE), length.out = 100)
    nls_mod <- models_exp[[sp]]
    if(is.null(nls_mod)) return(NULL)
    pred <- predict(nls_mod, newdata = data.frame(x = x_seq))
    data.frame(species = sp, x = x_seq, pred = pred)
  })
  pred_exp_df <- do.call(rbind, pred_exp)
  
  # Combine all predictions for Panel A
  pred_all <- bind_rows(pred_linear_df, pred_exp_df)
  
  p_combined <- ggplot() +
    geom_point(data = data_clean, aes(x = x, y = avg_value, color = species)) +
    geom_line(data = pred_all, aes(x = x, y = pred, color = species), linewidth = 1) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "black", linewidth = 1) +
    annotate("text", x = min(data_clean$x, na.rm = TRUE), y = threshold,
             label = "median", hjust = -5.2, vjust = -0.3, fontface = "italic", size = 6) +
    scale_color_manual(values = cb_palette, name = "") +
    labs(x = "transpiration deficit (mm)", y = "NDVI quantiles (rank)") +
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
  
  #################################################
  # Additional Analysis for Panel B & C: x50 and Slope for All Species
  #################################################
  
  # For linear models: calculate x50 and absolute slope (and p-value from slope)
  lin_stats <- lapply(linear_species, function(sp) {
    mod <- models_linear[[sp]]
    if(is.null(mod)) return(NULL)
    summ <- summary(mod)$coefficients
    intercept <- summ["(Intercept)", "Estimate"]
    slope <- summ["x", "Estimate"]
    se_slope <- summ["x", "Std. Error"]
    p_val_lin <- summ["x", "Pr(>|t|)"]
    x50_lin <- (threshold - intercept) / slope
    abs_slope_lin <- abs(slope)
    r2_lin <- summary(mod)$r.squared
    data.frame(species = sp, model = "linear", x50 = x50_lin, abs_slope = abs_slope_lin,
               r_squared = r2_lin, se = se_slope, p_val = p_val_lin)
  })
  
  # For exponential models: calculate x50 and absolute slope with deltaMethod and compute p-value
  exp_stats <- lapply(exp_species, function(sp) {
    mod <- models_exp[[sp]]
    if(is.null(mod)) return(NULL)
    coefs <- coef(mod)
    a_val <- coefs["a"]
    b_val <- coefs["b"]
    c_val <- coefs["c"]
    x50_exp <- ifelse((threshold - a_val) > 0 & b_val > 0,
                      -log((threshold - a_val)/b_val) / c_val,
                      NA)
    slope50_exp <- -c_val * (threshold - a_val)
    abs_slope_exp <- abs(slope50_exp)
    sp_data <- data_clean %>% filter(species == sp)
    fitted_vals <- predict(mod, newdata = sp_data)
    r2_exp <- 1 - sum((sp_data[[value_col]] - fitted_vals)^2) / 
      sum((sp_data[[value_col]] - mean(sp_data[[value_col]]))^2)
    dm_result <- deltaMethod(mod, paste0("c*(a - ", threshold, ")"),
                             parameterNames = c("a", "b", "c"))
    se_exp <- dm_result$SE
    df_exp <- summary(mod)$df[2]
    t_val_exp <- slope50_exp / se_exp
    p_val_exp <- 2 * (1 - pt(abs(t_val_exp), df_exp))
    data.frame(species = sp, model = "exponential", x50 = x50_exp, abs_slope = abs_slope_exp,
               r_squared = r2_exp, se = se_exp, p_val = p_val_exp)
  })
  
  stats_all <- rbind(do.call(rbind, lin_stats), do.call(rbind, exp_stats))
  stats_all$species <- factor(stats_all$species, levels = c(linear_species, exp_species))
  
  # Panel B: Bar plot of x50 (transpiration deficit) for all species
  p_x50 <- ggplot(stats_all, aes(x = species, y = x50, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    # geom_text(aes(label = sprintf("%.2f", x50), y = x50/2),
    #          color = "black", size = 4) +
    scale_fill_manual(values = cb_palette) +
    labs(x = "", y = "transpiration deficit (mm)") +
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
    labs(caption = "(b)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0.5))
  
  # Panel C: Bar plot of absolute slope for all species with error bars and R² annotations.
  # The annotation includes an asterisk if p < 0.05.
  p_slope <- ggplot(stats_all, aes(x = species, y = abs_slope, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, abs_slope - se), ymax = abs_slope + se), width = 0.2) +
    geom_text(aes(label = if_else(p_val < 0.05,
                                  sprintf("%.2f*", r_squared),
                                  sprintf("%.2f", r_squared)),
                  y = abs_slope/2),
              color = "black", size = 5) +
    scale_fill_manual(values = cb_palette) +
    labs(x = "", y = "absolute slope") +
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
    labs(caption = "(c)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0.5))
  
  # Combine panels: Panel (a) on the left; Panels (b) and (c) in the right column
  final_plot <- (p_combined + (p_x50 / p_slope)) +
    plot_layout(widths = c(2, 1), guides = "collect") &
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 14),
      legend.key = element_rect(fill = "white", color = NA),
      legend.background = element_blank(),         # <- No background
      legend.box.background = element_blank()      # <- No box border
    )
  print(final_plot)
  
  #### Plateau point calculation for exponential models (Spruce & Pine only)
  plateau_df <- lapply(exp_species, function(sp) {
    mod <- models_exp[[sp]]
    if (is.null(mod)) return(NULL)
    coefs <- coef(mod)
    c_val <- coefs["c"]
    x_plateau <- -log(0.05) / c_val
    data.frame(species = sp, x_plateau = x_plateau)
  }) %>% bind_rows()
  
  print("Estimated transpiration deficit (x_plateau) where NDVI is ~95% of its asymptotic minimum:")
  print(plateau_df)
  
  # Save the combined final plot
  ggsave(filename = output_figure, plot = final_plot, device = "png", width = 10, height = 8, dpi = 300)
}

plot_Quantiles_TDiff_exp_slope_coeff <- function(data, coef_fig, output_figure) {
  # Process data using the custom NDVI_TDiffbin function and remove missing values
  data <- NDVI_TDiffbin(data)
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
    mutate(label = if_else(round(pvalue, 2) <= 0.05, "*", sprintf("%.2f", pvalue)))
  
  
  # Create the faceted coefficient bar plot with p-value labels
  p_coeff <- ggplot(coeffs_long, aes(x = species, y = Value, fill = species)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(aes(label = label, y = Value/2),
              position = position_dodge(width = 0.9),
              vjust = 0.5, color = "black", size = 4) +
    scale_fill_manual(values = cb_palette, name = "") +
    facet_wrap(~ Coefficient, scales = "free_y") +
    labs(title = "NDVI quantiles ~ transpiration deficit",
         subtitle = expression(NDVI == a + b * e^{-~c *x} + epsilon),
         x = "",
         y = "Coefficient Value") +
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
  print(p_coeff)
  
  # Save the coefficient plot
  ggsave(filename = coef_fig, plot = p_coeff, device = "png", width = 10, height = 8, dpi = 300)
  
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
             label = "median", hjust = -5.2, vjust = -0.3, fontface = "italic", size = 6) +
    scale_color_manual(values = cb_palette) +
    labs(x = "transpiration deficit (mm)", y = "NDVI quantiles (rank)") +
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
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0),
          plot.caption.position = "plot") +
    coord_cartesian(clip = "off")
  
  # Panel B: Bar plot of x50 values
  p_x50 <- ggplot(coef_df, aes(x = species, y = x50, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    labs(x = "", y = "transpiration deficit") +
    scale_fill_manual(values = cb_palette) +
    expand_limits(y = 0) +
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
    labs(x = "", y = "absolue slope") +
    scale_fill_manual(values = cb_palette) +
    expand_limits(y = 0) +
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
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0),
          plot.caption.position = "plot")
  
  # Combine panels: Panel (a) on the left; Panel (b) (top) and Panel (c) (bottom) in right column
  final_plot <- (p_combined + (p_x50 / p_slope)) +
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
  # Save the combined final plot
  ggsave(filename = output_figure, plot = final_plot, device = "png", width = 10, height = 8, dpi = 300)
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
  
  ## Panel A: Mixed-Effects Model Plot (same as before)
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
  
  # Create the plot
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
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0.5))+
    coord_cartesian(clip = "off")
  
  ## Analytical Calculation of x50 and Slope (modified with delta method)
  
  # Extract species-specific coefficients
  species_coef <- coef(model)$species %>%
    tibble::rownames_to_column("species") %>%
    rename(a = `(Intercept)`, 
           b = `poly(bin_median, 2, raw = TRUE)1`, 
           c = `poly(bin_median, 2, raw = TRUE)2`)
  
  # For each species, fit a separate model to use with deltaMethod
  species_list <- species_levels
  stats_list <- lapply(species_list, function(sp) {
    # Fit model for this species only
    df_sp <- TDiff_PSIbin_df %>% filter(species == sp)
    mod_sp <- lm(avg_transpiration_deficit ~ poly(bin_median, 2, raw = TRUE), data = df_sp)
    
    # Get coefficients for this species
    coefs <- coef(mod_sp)
    a <- coefs[1]; b <- coefs[2]; c <- coefs[3]
    
    # Calculate x50 (where curve crosses mean transpiration deficit)
    disc <- b^2 - 4 * c * (a - line_val)
    x50 <- ifelse(disc >= 0, (-b + sqrt(disc)) / (2 * c), NA_real_)
    
    # Use delta method to get SE for slope at x50 (derivative: b + 2*c*x)
    if(!is.na(x50)) {
      # Expression for slope at x50
      slope_expr <- paste0("b + 2 * c * ", x50)
      dm_result <- deltaMethod(mod_sp, slope_expr, parameterNames = c("a", "b", "c"))
      
      # Calculate p-value
      t_val <- (b + 2 * c * x50) / dm_result$SE
      df_resid <- df.residual(mod_sp)
      p_val <- 2 * (1 - pt(abs(t_val), df_resid))
      
      # Calculate R-squared
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
  
  # Combine plots
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
  
  ## Coefficient Plot (same as before)
  # Extract the species-specific (conditional) coefficients.
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
  species_list <- species_levels
  coeff_stats_list <- lapply(species_list, function(sp) {
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
  
  return(list(combined_plot = final_plot, stats_df = stats_df))
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
    scale_color_manual(values = cb_palette, name = "") +
    labs(x = "soil water potential",
         y = "transpiration deficit") +
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
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0.5))+
    coord_cartesian(clip = "off")
  
  ## Panel B: Bar Plot of x-values at Mean Transpiration Deficit
  
  x_at_mean <- pred_data %>%
    group_by(species) %>%
    summarise(x_at_mean = bin_median[which.min(abs(predicted - line_val))])
  x_at_mean$species <- factor(x_at_mean$species, levels = species_levels)
  
  p_bar_x <- ggplot(x_at_mean, aes(x = species, y = x_at_mean, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = cb_palette, name = "") +
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
    scale_fill_manual(values = cb_palette, name = "") +
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
  
  ## Combine Panels: Panel A on the left, Panels B (top) and C (bottom) on the right.
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
    scale_fill_manual(values = cb_palette, name = "") +
    labs(title = "transpiration deficit ~ soil water potential", 
         subtitle = expression(NDVI == a + b*x + c*x^2 + d*x^3),
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
  
  print(coeff_data)
  # Return both plots as a list (if desired)
  return(list(combined_plot = final_plot, coeff_plot = plot_coeff))
}

plot_density_Quantiles_PSI_linear <- function(data, swp_bin_width = 50, output_path) {
  library(dplyr)
  library(tidyr)      # Added to provide drop_na() function
  library(ggplot2)
  
  # Define species order
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  
  # Filter data for the four species and set factor levels
  species_df <- data %>% 
    filter(species %in% species_order) %>%
    mutate(species = factor(species, levels = species_order))
  
  # Bin Quantiles (1-22) and Soil Water Potential (using swp_bin_width)
  species_binned <- species_df %>%
    mutate(
      quantile_bin = cut(Quantiles, 
                         breaks = seq(0, 22, 1), 
                         include.lowest = TRUE, 
                         right = FALSE),
      swp_bin = cut(soil_water_potential, 
                    breaks = seq(floor(min(soil_water_potential) / swp_bin_width) * swp_bin_width,
                                 ceiling(max(soil_water_potential) / swp_bin_width) * swp_bin_width,
                                 swp_bin_width),
                    include.lowest = TRUE, 
                    right = FALSE)
    ) %>%
    drop_na(quantile_bin, swp_bin)
  
  # Calculate density percentages for each bin combination, normalizing by each swp_bin (column)
  density_df <- species_binned %>%
    group_by(species, quantile_bin, swp_bin) %>%
    summarize(count = n(), .groups = "drop") %>%
    group_by(species, swp_bin) %>%  # Group by species and soil water potential bin
    mutate(density_percent = (count / sum(count)) * 100) %>%
    ungroup()
  
  # Join densities back to binned data
  species_density <- species_binned %>%
    left_join(density_df, by = c("species", "quantile_bin", "swp_bin"))
  
  # Calculate maximum density and create five equally spaced breaks for the legend
  max_density <- ceiling(max(species_density$density_percent, na.rm = TRUE))
  density_breaks <- seq(0, max_density, length.out = 5)
  
  # Create the plot with fitted linear regression lines and facets for each species
  p <- ggplot(species_density, aes(x = soil_water_potential, y = Quantiles, color = density_percent)) +
    geom_point(size = 2.5, alpha = 0.9) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1) +
    scale_color_viridis_c(
      name = "density (%)", 
      option = "inferno", 
      direction = -1, 
      trans = "sqrt",
      breaks = density_breaks
    ) +
    labs(
      title = "NDVI quantiles ~ soil water potential",
      x = "soil water potential",
      y = "NDVI quantiles"
    ) +
    facet_wrap(~species, ncol = 2) +
    theme_minimal() +
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
      legend.title = element_text(size = 14),  # Match legend text size
      legend.key.width = unit(1.8, "cm"),      # Widen legend keys slightly
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    )
  
  # Print and save the plot to the specified output path
  print(p)
  ggsave(output_path, plot = p, width = 12, height = 8, dpi = 300)
}

plot_density_Quantiles_TDiff_linear <- function(data, tdiff_bin_width = 3, output_path) {
  
  library(dplyr)
  library(tidyr)      # Added to provide drop_na() function
  library(ggplot2)
  
  # Define species order
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  
  # Filter data for the four species and set factor levels
  species_df <- data %>% 
    filter(species %in% species_order) %>%
    mutate(species = factor(species, levels = species_order))
  
  # Bin Quantiles (1-22) and Transpiration Deficit (using tdiff_bin_width)
  species_binned <- species_df %>%
    mutate(
      quantile_bin = cut(Quantiles, breaks = seq(0, 22, 1), include.lowest = TRUE, right = FALSE),
      tdiff_bin = cut(transpiration_deficit, 
                      breaks = seq(floor(min(transpiration_deficit) / tdiff_bin_width) * tdiff_bin_width,
                                   ceiling(max(transpiration_deficit) / tdiff_bin_width) * tdiff_bin_width,
                                   tdiff_bin_width),
                      include.lowest = TRUE, right = FALSE)
    ) %>%
    drop_na(quantile_bin, tdiff_bin)
  
  # Calculate density percentages for each bin combination,
  # normalizing by each transpiration deficit bin (tdiff_bin) per species
  density_df <- species_binned %>%
    group_by(species, quantile_bin, tdiff_bin) %>%
    summarize(count = n(), .groups = "drop") %>%
    group_by(species, tdiff_bin) %>%  # Group by species and tdiff_bin
    mutate(density_percent = (count / sum(count)) * 100) %>%
    ungroup()
  
  # Join densities back to binned data
  species_density <- species_binned %>%
    left_join(density_df, by = c("species", "quantile_bin", "tdiff_bin"))
  
  # Calculate maximum density and create five equally spaced breaks for the legend
  max_density <- ceiling(max(species_density$density_percent, na.rm = TRUE))
  density_breaks <- seq(0, max_density, length.out = 5)
  
  # Create the plot with fitted linear regression lines and facets for each species
  p <- ggplot(species_density, aes(x = transpiration_deficit, y = Quantiles, color = density_percent)) +
    geom_point(size = 2.5, alpha = 0.9) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1) +
    scale_color_viridis_c(
      name = "density (%)", 
      option = "inferno", 
      direction = -1, 
      trans = "sqrt",
      breaks = density_breaks
    ) +
    labs(
      title = "NDVI quantiles ~ transpiration deficit",
      x = "transpiration deficit",
      y = "NDVI quantiles"
    ) +
    facet_wrap(~species, ncol = 2) +
    theme_minimal() +
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
      legend.title = element_text(size = 14),  # Match legend text size
      legend.key.width = unit(1.8, "cm"),      # Widen legend keys slightly
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    )
  
  # Print and save the plot to the specified output path
  print(p)
  ggsave(output_path, plot = p, width = 12, height = 8, dpi = 300)
}

plot_species_distribution(all_results_df, "results/key_displays/species_map.png")

plot_moran_I_distance("results_moran/data", "results/key_displays/moran_distance.png")

plot_combined_AIC_R2(all_results_df,  "results/key_displays/AIC.png")

plot_time_series_and_correlation_combined(all_results_df,
                                          "results/key_displays/time_series_Quantiles_PSI_TDiff_species.png")

plot_yearly_mean_linear_coeffs(all_results_df, "results/key_displays/yearly_mean_linear_coeffs.png")

plot_NDVI_Q_PSIbin_exp_slope_negPSI(all_results_df, 
                                    "results/key_displays/NDVI_Q_PSIbin_exp_coeff_negPSI.png",
                                    "results/key_displays/NDVI_Q_PSIbin_exp_slope_negPSI.png")

plot_Quantiles_TDiff_exp_linear_slope_coeff(all_results_df,
                                            "results/key_displays/NDVI_Q_TDiffbin_exp_linear_coeff.png",
                                            "results/key_displays/NDVI_Q_TDiffbin_exp_linear_slope.png")

plot_Quantiles_TDiff_exp_slope_coeff(all_results_df, 
                                     "results/key_displays/NDVI_Q_TDiffbin_exp_coeff.png",
                                     "results/key_displays/NDVI_Q_TDiffbin_exp_slope.png")

plot_TDiff_PSIbin_poly_2_slope(all_results_df, 
                               "results/key_displays/TDiff_PSIbin_poly_2_coeff.png",
                               "results/key_displays/TDiff_PSIbin_poly_2_slope.png")

plot_TDiff_PSIbin_poly_3_slope(all_results_df, 
                               "results/key_displays/TDiff_PSIbin_poly_3_coeff.png",
                               "results/key_displays/TDiff_PSIbin_poly_3_slope.png")

plot_density_Quantiles_PSI_linear(all_results_df, 
                                  output_path = "results/key_displays/Quantiles_PSI_density_0728.png")

plot_density_Quantiles_TDiff_linear(all_results_df, 
                                    output_path = "results/key_displays/Quantiles_TDiff_density_0728.png")

load("results/Data/All_species_month_year_Quantiles_PSI_TDiff.RData")

NDVI_PSIbin_month <- function(df, bin_width = 50) {
  library(dplyr)
  
  # Identify the correct column dynamically (here we use "Quantiles")
  value_column <- if ("Quantiles" %in% names(df)) "Quantiles" else "Proportions"
  
  # Define bin breaks dynamically based on soil water potential range
  psi_min <- floor(min(df$soil_water_potential, na.rm = TRUE))
  psi_max <- ceiling(max(df$soil_water_potential, na.rm = TRUE))
  bin_breaks <- seq(psi_min, psi_max, by = bin_width)
  
  df <- df %>%
    mutate(PSI_bin = cut(soil_water_potential, breaks = bin_breaks, include.lowest = TRUE, right = FALSE))
  
  # Count total pixels per species
  species_totals <- df %>%
    group_by(species) %>%
    summarise(total_count = n(), .groups = 'drop')
  
  # Compute mean for the chosen value column per species, per month, and per PSI_bin
  meanNDVI_PSIbin_species <- df %>%
    group_by(species, month, PSI_bin) %>%
    summarise(
      avg_value = mean(.data[[value_column]], na.rm = TRUE),
      count = n(),
      .groups = 'drop'
    ) %>%
    left_join(species_totals, by = "species") %>%
    mutate(
      percent_of_total = (count / total_count) * 100
    ) %>%
    filter(percent_of_total >= 0.1) %>%  # Only keep bins >= 0.1% of total pixels per species
    mutate(
      bin_median = sapply(as.character(PSI_bin), function(bin_label) {
        nums <- as.numeric(strsplit(gsub("\\[|\\]|\\(|\\)", "", bin_label), ",")[[1]])
        mean(nums)
      })
    ) %>%
    select(species, month, PSI_bin, bin_median, avg_value)
  
  return(meanNDVI_PSIbin_species)
}

NDVI_TDiffbin_month <- function(df, bin_width = 3) {
  library(dplyr)
  
  # Identify the correct column dynamically (either "Quantiles" or "Proportions")
  value_column <- if ("Quantiles" %in% names(df)) "Quantiles" else "Proportions"
  
  # Define bin breaks dynamically based on transpiration deficit range
  tdiff_min <- floor(min(df$transpiration_deficit, na.rm = TRUE))
  tdiff_max <- ceiling(max(df$transpiration_deficit, na.rm = TRUE))
  bin_breaks <- seq(tdiff_min, tdiff_max, by = bin_width)
  
  df <- df %>%
    mutate(TDiff_bin = cut(transpiration_deficit, breaks = bin_breaks, include.lowest = TRUE, right = FALSE))
  
  # Count total pixels per species
  species_totals <- df %>%
    group_by(species) %>%
    summarise(total_count = n(), .groups = 'drop')
  
  # Compute mean for the chosen value column per species, per month, and per TDiff_bin
  meanNDVI_TDiff_species <- df %>%
    group_by(species, month, TDiff_bin) %>%
    summarise(
      avg_value = mean(.data[[value_column]], na.rm = TRUE),
      count = n(),
      .groups = 'drop'
    ) %>%
    left_join(species_totals, by = "species") %>%
    mutate(
      percent_of_total = (count / total_count) * 100
    ) %>%
    filter(percent_of_total >= 0.01) %>%  # Keep bins >= 0.01% of species' total pixels
    mutate(
      bin_median = sapply(as.character(TDiff_bin), function(bin_label) {
        nums <- as.numeric(strsplit(gsub("\\[|\\]|\\(|\\)", "", bin_label), ",")[[1]])
        mean(nums)
      })
    ) %>%
    select(species, month, TDiff_bin, bin_median, avg_value)
  
  return(meanNDVI_TDiff_species)
}

plot_Quantiles_PSI_month_species_linear_orig <- function(data, figure_output = NULL, figure_output2 = NULL) {
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
    # geom_smooth(method = "lm", se = FALSE, linetype = "dashed", aes(group = month)) +
    geom_smooth(method = "lm", se = FALSE, aes(group = month)) +
    facet_wrap(~ species, scales = "free_x") +
    scale_color_viridis_d(option = "C", begin = 0.2, end = 0.8) +
    labs(
      x = "Soil Water Potential",
      y = "NDVI quantiles",
      color = "",
      title = "NDVI quantiles ~ soil water potential",
      subtitle = expression(NDVI == a + bx + epsilon),
    ) +
    theme_minimal() +
    theme(
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
    pivot_longer(cols = c(intercept, slope), names_to = "parameter", values_to = "value") %>%
    mutate(parameter = dplyr::recode(parameter, intercept = "a", slope = "b"))
  
  p2 <- ggplot(params_long, aes(x = month, y = value, fill = month)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    facet_grid(parameter ~ species, scales = "free_y") +
    scale_fill_viridis_d(option = "C", begin = 0.2, end = 0.8) +
    labs(
      x = "",
      y = "coefficient value",
      fill = "",
      title = "NDVI quantile ~ soil water potential",
      subtitle = expression(NDVI == a + bx + epsilon),
    ) +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      plot.subtitle = element_text(hjust = 0.5, size =14),
      axis.title = element_text(face = "bold", size = 16),
      axis.text.y = element_text(color = "black", size = 14),
      axis.text.x = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.text = element_text(size = 14),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
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

plot_Quantiles_TDiff_month_species_linear_orig <- function(data, figure_output = NULL, figure_output2 = NULL) {
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
  
  p1 <- ggplot(data, aes(x = transpiration_deficit, y = Quantiles, color = month)) +
    # geom_point() +
    # geom_smooth(method = "lm", se = FALSE, linetype = "dashed", aes(group = month)) +
    geom_smooth(method = "lm", se = FALSE, aes(group = month)) +
    facet_wrap(~ species, scales = "free_x") +
    scale_color_viridis_d(option = "C", begin = 0.2, end = 0.8) +
    labs(
      x = "transpiration deficit",
      y = "NDVI quantiles",
      color = "",
      title = "NDVI quantiles ~ transpiration deficit",
      subtitle = expression(NDVI == a + bx + epsilon),
    ) +
    theme_minimal() +
    theme(
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
  
  ## Figure 2: Bar plot of regression parameters (intercept and slope) for each species and month
  
  # Compute regression parameters for each species and month from the original data
  params <- data %>%
    group_by(species, month) %>%
    group_modify(~ {
      mod <- lm(Quantiles ~ transpiration_deficit, data = .x)
      tibble(intercept = coef(mod)[1], slope = coef(mod)[2])
    }) %>%
    ungroup()
  
  # Pivot the parameters into a long format so both intercept and slope can be plotted
  params_long <- params %>%
    pivot_longer(cols = c(intercept, slope), names_to = "parameter", values_to = "value") %>%
    mutate(parameter = dplyr::recode(parameter, intercept = "a", slope = "b"))
  
  p2 <- ggplot(params_long, aes(x = month, y = value, fill = month)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    facet_grid(parameter ~ species, scales = "free_y") +
    scale_fill_viridis_d(option = "C", begin = 0.2, end = 0.8) +
    labs(
      x = "",
      y = "coefficient value",
      fill = "",
      title = "NDVI quantile ~ transpiration deficit",
      subtitle = expression(NDVI == a + bx + epsilon),
    ) +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      plot.subtitle = element_text(hjust = 0.5, size =14),
      axis.title = element_text(face = "bold", size = 16),
      axis.text.y = element_text(color = "black", size = 14),
      axis.text.x = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.text = element_text(size = 14),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
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

# depth 50
plot_Quantiles_PSI_month_species_linear_orig(final_df, 
                                             figure_output = "results/key_displays/Quantiles_PSI_linear_slope_monthly.png",
                                             figure_output2 = "results/key_displays/Quantiles_PSI_linear_coeff_monthly.png")

plot_Quantiles_TDiff_month_species_linear_orig(final_df, 
                                               figure_output = "results/key_displays/Quantiles_TDiff_linear_slope_monthly.png",
                                               figure_output2 = "results/key_displays/Quantiles_TDiff_linear_coeff_monthly.png")

# depth 100
load("results/Data/final_df_depth100.RData")

plot_Quantiles_PSI_month_species_linear_orig(final_df_depth100, 
                                             figure_output = "results/key_displays/Quantiles_PSI_linear_slope_monthly_depth100.png",
                                             figure_output2 = "results/key_displays/Quantiles_PSI_linear_coeff_monthly_depth100.png")

plot_Quantiles_TDiff_month_species_linear_orig(final_df_depth100, 
                                               figure_output = "results/key_displays/Quantiles_TDiff_linear_slope_monthly_depth100.png",
                                               figure_output2 = "results/key_displays/Quantiles_TDiff_linear_coeff_monthly_depth100.png")

load("results/Data/final_df_depth150.RData")
# depth 150
plot_Quantiles_PSI_month_species_linear_orig(final_df_depth150, 
                                             figure_output = "results/key_displays/Quantiles_PSI_linear_slope_monthly_depth150.png",
                                             figure_output2 = "results/key_displays/Quantiles_PSI_linear_coeff_monthly_depth150.png")

plot_Quantiles_TDiff_month_species_linear_orig(final_df_depth150, 
                                               figure_output = "results/key_displays/Quantiles_TDiff_linear_slope_monthly_depth150.png",
                                               figure_output2 = "results/key_displays/Quantiles_TDiff_linear_coeff_monthly_depth150.png")

### End of July and August ###
library(dplyr)

plot_Quantiles_PSI_exp_linear_slope_coeff <- function(data, combined_coef_fig, output_figure) {
  # Process data using the custom NDVI_TDiffbin function and remove missing rows
  data <- NDVI_PSIbin(data)
  data <- na.omit(data)
  
  # Load required libraries
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(patchwork)
  library(purrr)
  library(nlme)    # for consistency (even though lm is used for linear)
  library(car)     # for deltaMethod
  library(broom)
  
  # Define the value column and species groups.
  # Use linear regression for Oak and Beech,
  # and the exponential model for Spruce and Pine.
  value_col <- "avg_value"
  linear_species <- c("Oak", "Beech")
  exp_species <- c("Spruce", "Pine")
  data$species <- factor(data$species, levels = c(linear_species, exp_species))
  
  # Define the color palette (same for all species)
  cb_palette <- c("Oak"   = "#E69F00",  # Orange
                  "Beech" = "#0072B2",  # Deep blue
                  "Spruce"= "#009E73",  # Bluish-green
                  "Pine"  = "#F0E442")  # Yellow
  
  # For transpiration deficit, set x = bin_median
  data <- data %>% mutate(x = bin_median)
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  
  # Set threshold (using a fixed value as before)
  threshold <- 11.5
  
  ##########################
  # Model Fitting by Group #
  ##########################
  
  # Linear models for Oak and Beech
  models_linear <- list()
  for(sp in linear_species) {
    sp_data <- data_clean %>% filter(species == sp)
    lm_model <- lm(avg_value ~ x, data = sp_data)
    models_linear[[sp]] <- lm_model
  }
  
  # Exponential models for Spruce and Pine
  start_list <- list(a = 5, b = 3, c = 0.001)
  control_params <- nls.control(maxiter = 1200, minFactor = 1e-09)
  models_exp <- list()
  for(sp in exp_species) {
    sp_data <- data_clean %>% filter(species == sp)
    
    # Check that sp_data has enough data
    if(nrow(sp_data) < 5) {
      message("Not enough data for species ", sp)
      models_exp[[sp]] <- NULL
      next
    }
    
    exp_model <- tryCatch({
      nls(avg_value ~ a + b * exp(c * x), 
          data = sp_data, 
          start = start_list, 
          control = control_params)
    }, error = function(e) {
      message("Error for species ", sp, ": ", e$message)
      NULL
    })
    
    models_exp[[sp]] <- exp_model
  }
  
  ##############################################
  # Coefficient Extraction and Plotting - Linear
  ##############################################
  
  lin_coef_list <- lapply(linear_species, function(sp) {
    mod <- models_linear[[sp]]
    if(is.null(mod)) return(NULL)
    summ <- summary(mod)$coefficients
    df <- data.frame(
      species = sp,
      Coefficient = rownames(summ),
      Value = summ[, "Estimate"],
      pvalue = summ[, "Pr(>|t|)"],
      stringsAsFactors = FALSE
    )
    # Rename coefficients for clarity using dplyr::recode
    df$Coefficient <- dplyr::recode(df$Coefficient, "(Intercept)" = "a", "x" = "b")
    return(df)
  })
  lin_coef_df <- do.call(rbind, lin_coef_list)
  lin_coef_df$species <- factor(lin_coef_df$species, levels = linear_species)
  lin_coef_df <- lin_coef_df %>% 
    mutate(label = if_else(pvalue < 0.05, "*", sprintf("%.2f", pvalue)))
  
  p_coeff_linear <- ggplot(lin_coef_df, aes(x = species, y = Value, fill = species)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(aes(label = label, y = Value/2),
              position = position_dodge(width = 0.9),
              vjust = 0.5, color = "black", size = 4) +
    scale_fill_manual(values = cb_palette, name = "") +
    facet_wrap(~ Coefficient, scales = "free_y") +
    labs(title = "linear models",
         # subtitle = expression(italic(NDVI[Quantiles]) == a + b * x + epsilon),
         subtitle = expression(NDVI == a + b * x + epsilon),
         x = "",
         y = "coefficient value") +
    theme_minimal() +
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
    ) +
    labs(caption = "(a)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0.5))
  
  ##############################################
  # Coefficient Extraction and Plotting - Exponential
  ##############################################
  
  exp_coef_list <- lapply(exp_species, function(sp) {
    mod <- models_exp[[sp]]
    if (is.null(mod)) return(NULL)
    summ <- summary(mod)$coefficients
    df <- data.frame(
      species = sp,
      Coefficient = rownames(summ),
      Value = summ[, "Estimate"],
      pvalue = summ[, "Pr(>|t|)"],
      stringsAsFactors = FALSE
    )
    # Convert coefficient names to lowercase (a, b, c)
    df$Coefficient <- tolower(df$Coefficient)
    return(df)
  })
  # Remove any NULL entries:
  exp_coef_list <- exp_coef_list[!sapply(exp_coef_list, is.null)]
  exp_coef_df <- do.call(rbind, exp_coef_list)
  exp_coef_df$species <- factor(exp_coef_df$species, levels = exp_species)
  exp_coef_df <- exp_coef_df %>% 
    mutate(label = if_else(pvalue < 0.05, "*", sprintf("%.2f", pvalue)))
  
  p_coeff_exp <- ggplot(exp_coef_df, aes(x = species, y = Value, fill = species)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(aes(label = label, y = Value/2),
              position = position_dodge(width = 0.9),
              vjust = 0.5, color = "black", size = 4) +
    scale_fill_manual(values = cb_palette, name = "") +
    facet_wrap(~ Coefficient, scales = "free_y") +
    labs(title = "exponential models",
         # subtitle = expression(italic(NDVI[Quantiles]) == a + b * e^{-c * italic(x)} + epsilon),
         subtitle = expression(NDVI == a + b * e^{c * italic(x)} + epsilon),
         x = "",
         y = "") +
    theme_minimal() +
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
    ) +
    labs(caption = "(b)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0.5))
  
  # Combine the coefficient plots into one figure with panel labels (a) and (b)
  combined_coeff <- p_coeff_linear + p_coeff_exp 
  print(combined_coeff)
  
  # Save the combined coefficient plot
  ggsave(filename = combined_coef_fig, plot = combined_coeff, device = "png", width = 10, height = 8, dpi = 300)
  
  #################################################
  # Combined Final Plot with Model Predictions (Panel A)
  #################################################
  
  # Generate predictions for linear species (Oak & Beech)
  pred_linear <- lapply(linear_species, function(sp) {
    sp_data <- data_clean %>% filter(species == sp)
    x_seq <- seq(min(sp_data$x, na.rm = TRUE), max(sp_data$x, na.rm = TRUE), length.out = 100)
    lm_mod <- models_linear[[sp]]
    pred <- predict(lm_mod, newdata = data.frame(x = x_seq))
    data.frame(species = sp, x = x_seq, pred = pred)
  })
  pred_linear_df <- do.call(rbind, pred_linear)
  
  # Generate predictions for exponential species (Spruce & Pine)
  pred_exp <- lapply(exp_species, function(sp) {
    sp_data <- data_clean %>% filter(species == sp)
    x_seq <- seq(min(sp_data$x, na.rm = TRUE), max(sp_data$x, na.rm = TRUE), length.out = 100)
    nls_mod <- models_exp[[sp]]
    if(is.null(nls_mod)) return(NULL)
    pred <- predict(nls_mod, newdata = data.frame(x = x_seq))
    data.frame(species = sp, x = x_seq, pred = pred)
  })
  pred_exp_df <- do.call(rbind, pred_exp)
  
  # Combine all predictions for Panel A
  pred_all <- bind_rows(pred_linear_df, pred_exp_df)
  
  p_combined <- ggplot() +
    geom_point(data = data_clean, aes(x = x, y = avg_value, color = species)) +
    geom_line(data = pred_all, aes(x = x, y = pred, color = species), linewidth = 1) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "black", linewidth = 1) +
    annotate("text", x = min(data_clean$x, na.rm = TRUE), y = threshold,
             label = "median", hjust = -0.1, vjust = -0.3, fontface = "italic", size = 6) +
    scale_color_manual(values = cb_palette, name = "") +
    labs(x = "soil water potential (kPa)", y = "NDVI quantiles (rank)") +
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
  
  #################################################
  # Additional Analysis for Panel B & C: x50 and Slope for All Species
  #################################################
  
  # For linear models: calculate x50 and absolute slope (and p-value from slope)
  lin_stats <- lapply(linear_species, function(sp) {
    mod <- models_linear[[sp]]
    if(is.null(mod)) return(NULL)
    summ <- summary(mod)$coefficients
    intercept <- summ["(Intercept)", "Estimate"]
    slope <- summ["x", "Estimate"]
    se_slope <- summ["x", "Std. Error"]
    p_val_lin <- summ["x", "Pr(>|t|)"]
    x50_lin <- (threshold - intercept) / slope
    abs_slope_lin <- abs(slope)
    r2_lin <- summary(mod)$r.squared
    data.frame(species = sp, model = "linear", x50 = x50_lin, abs_slope = abs_slope_lin,
               r_squared = r2_lin, se = se_slope, p_val = p_val_lin)
  })
  
  # For exponential models: calculate x50 and absolute slope with deltaMethod and compute p-value
  exp_stats <- lapply(exp_species, function(sp) {
    mod <- models_exp[[sp]]
    if(is.null(mod)) return(NULL)
    coefs <- coef(mod)
    a_val <- coefs["a"]
    b_val <- coefs["b"]
    c_val <- coefs["c"]
    x50_exp <- ifelse((threshold - a_val) > 0 & b_val > 0,
                      log((threshold - a_val)/b_val) / c_val,
                      NA)
    slope50_exp <- c_val * (threshold - a_val)
    abs_slope_exp <- abs(slope50_exp)
    sp_data <- data_clean %>% filter(species == sp)
    fitted_vals <- predict(mod, newdata = sp_data)
    r2_exp <- 1 - sum((sp_data[[value_col]] - fitted_vals)^2) / 
      sum((sp_data[[value_col]] - mean(sp_data[[value_col]]))^2)
    dm_result <- deltaMethod(mod, paste0("c*(a - ", threshold, ")"),
                             parameterNames = c("a", "b", "c"))
    se_exp <- dm_result$SE
    df_exp <- summary(mod)$df[2]
    t_val_exp <- slope50_exp / se_exp
    p_val_exp <- 2 * (1 - pt(abs(t_val_exp), df_exp))
    data.frame(species = sp, model = "exponential", x50 = x50_exp, abs_slope = abs_slope_exp,
               r_squared = r2_exp, se = se_exp, p_val = p_val_exp)
  })
  
  stats_all <- rbind(do.call(rbind, lin_stats), do.call(rbind, exp_stats))
  stats_all$species <- factor(stats_all$species, levels = c(linear_species, exp_species))
  
  # Panel B: Bar plot of x50 (transpiration deficit) for all species
  p_x50 <- ggplot(stats_all, aes(x = species, y = x50, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    # geom_text(aes(label = sprintf("%.2f", x50), y = x50/2),
    #          color = "black", size = 4) +
    scale_fill_manual(values = cb_palette) +
    labs(x = "", y = "soil water potential (kPa)") +
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
    labs(caption = "(b)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0.5))
  
  # Panel C: Bar plot of absolute slope for all species with error bars and R² annotations.
  # The annotation includes an asterisk if p < 0.05.
  p_slope <- ggplot(stats_all, aes(x = species, y = abs_slope, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, abs_slope - se), ymax = abs_slope + se), width = 0.2) +
    geom_text(aes(label = if_else(p_val < 0.05,
                                  sprintf("%.2f*", r_squared),
                                  sprintf("%.2f", r_squared)),
                  y = abs_slope/2),
              color = "black", size = 5) +
    scale_fill_manual(values = cb_palette) +
    labs(x = "", y = "absolute slope") +
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
    labs(caption = "(c)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0.5))
  
  # Combine panels: Panel (a) on the left; Panels (b) and (c) in the right column
  final_plot <- (p_combined + (p_x50 / p_slope)) +
    plot_layout(widths = c(2, 1), guides = "collect") &
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 14),
      legend.key = element_rect(fill = "white", color = NA),
      legend.background = element_blank(),         # <- No background
      legend.box.background = element_blank()      # <- No box border
    )
  print(final_plot)
  
  #### Plateau point calculation for exponential models (Spruce & Pine only)
  plateau_df <- lapply(exp_species, function(sp) {
    mod <- models_exp[[sp]]
    if (is.null(mod)) return(NULL)
    coefs <- coef(mod)
    c_val <- coefs["c"]
    x_plateau <- -log(0.05) / c_val
    data.frame(species = sp, x_plateau = x_plateau)
  }) %>% bind_rows()
  
  print("Estimated soil water potential (x_plateau) where NDVI is ~95% of its asymptotic minimum:")
  print(plateau_df)
  
  # Save the combined final plot
  ggsave(filename = output_figure, plot = final_plot, device = "png", width = 10, height = 8, dpi = 300)
}

df <- final_df %>% filter(month %in% c("July", "August"))

plot_combined_AIC_R2(df,  "results/key_displays_July_August/AIC.png")

plot_time_series_and_correlation_combined(df,
                                          "results/key_displays_July_August/time_series_Quantiles_PSI_TDiff_species.png")

plot_yearly_mean_linear_coeffs(df, "results/key_displays_July_August/yearly_mean_linear_coeffs.png")

plot_Quantiles_PSI_exp_linear_slope_coeff(df, 
                                          "results/key_displays_July_August/NDVI_Q_PSIbin_exp_linear_coeff_negPSI.png",
                                          "results/key_displays_July_August/NDVI_Q_PSIbin_exp_linear_slope_negPSI.png")

plot_Quantiles_TDiff_exp_slope_coeff(df, 
                                     "results/key_displays_July_August/NDVI_Q_TDiffbin_exp_coeff.png",
                                     "results/key_displays_July_August/NDVI_Q_TDiffbin_exp_slope.png")

plot_TDiff_PSIbin_poly_2_slope(df, 
                               "results/key_displays_July_August/TDiff_PSIbin_poly_2_coeff.png",
                               "results/key_displays_July_August/TDiff_PSIbin_poly_2_slope.png")

plot_TDiff_PSIbin_poly_3_slope(df, 
                               "results/key_displays_July_August/TDiff_PSIbin_poly_3_coeff.png",
                               "results/key_displays_July_August/TDiff_PSIbin_poly_3_slope.png")

plot_density_Quantiles_PSI_linear(df, 
                                  output_path = "results/key_displays_July_August/Quantiles_PSI_density.png")

plot_density_Quantiles_TDiff_linear(df, 
                                    output_path = "results/key_displays_July_August/Quantiles_TDiff_density.png")

setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

library(dplyr)
load("results/Data/All_species_month_year_Quantiles_PSI_TDiff.RData")
df <- final_df %>% filter(month %in% c("July", "August"))

plot_Quantiles_PSI_exp_linear_slope_coeff <- function(data, combined_coef_fig, output_figure) {
  # Process data using the custom NDVI_TDiffbin function and remove missing rows
  data <- NDVI_PSIbin(data)
  data <- na.omit(data)
  
  # Load required libraries
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(patchwork)
  library(purrr)
  library(nlme)    # for consistency (even though lm is used for linear)
  library(car)     # for deltaMethod
  library(broom)
  
  # Define the value column and species groups.
  # Use linear regression for Oak and Beech,
  # and the exponential model for Spruce and Pine.
  value_col <- "avg_value"
  linear_species <- c("Oak", "Beech")
  exp_species <- c("Spruce", "Pine")
  data$species <- factor(data$species, levels = c(linear_species, exp_species))
  
  # Define the color palette (same for all species)
  cb_palette <- c("Oak"   = "#E69F00",  # Orange
                  "Beech" = "#0072B2",  # Deep blue
                  "Spruce"= "#009E73",  # Bluish-green
                  "Pine"  = "#F0E442")  # Yellow
  
  # For transpiration deficit, set x = bin_median
  data <- data %>% mutate(x = bin_median)
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  
  # Set threshold (using a fixed value as before)
  threshold <- 11.5
  
  ##########################
  # Model Fitting by Group #
  ##########################
  
  # Linear models for Oak and Beech
  models_linear <- list()
  for(sp in linear_species) {
    sp_data <- data_clean %>% filter(species == sp)
    lm_model <- lm(avg_value ~ x, data = sp_data)
    models_linear[[sp]] <- lm_model
  }
  
  # Exponential models for Spruce and Pine
  start_list <- list(a = 5, b = 3, c = 0.001)
  control_params <- nls.control(maxiter = 1200, minFactor = 1e-09)
  models_exp <- list()
  for(sp in exp_species) {
    sp_data <- data_clean %>% filter(species == sp)
    
    # Check that sp_data has enough data
    if(nrow(sp_data) < 5) {
      message("Not enough data for species ", sp)
      models_exp[[sp]] <- NULL
      next
    }
    
    exp_model <- tryCatch({
      nls(avg_value ~ a + b * exp(c * x), 
          data = sp_data, 
          start = start_list, 
          control = control_params)
    }, error = function(e) {
      message("Error for species ", sp, ": ", e$message)
      NULL
    })
    
    models_exp[[sp]] <- exp_model
  }
  
  ##############################################
  # Coefficient Extraction and Plotting - Linear
  ##############################################
  
  lin_coef_list <- lapply(linear_species, function(sp) {
    mod <- models_linear[[sp]]
    if(is.null(mod)) return(NULL)
    summ <- summary(mod)$coefficients
    df <- data.frame(
      species = sp,
      Coefficient = rownames(summ),
      Value = summ[, "Estimate"],
      pvalue = summ[, "Pr(>|t|)"],
      stringsAsFactors = FALSE
    )
    # Rename coefficients for clarity using dplyr::recode
    df$Coefficient <- dplyr::recode(df$Coefficient, "(Intercept)" = "a", "x" = "b")
    return(df)
  })
  lin_coef_df <- do.call(rbind, lin_coef_list)
  lin_coef_df$species <- factor(lin_coef_df$species, levels = linear_species)
  lin_coef_df <- lin_coef_df %>% 
    mutate(label = if_else(pvalue < 0.05, "*", sprintf("%.2f", pvalue)))
  
  p_coeff_linear <- ggplot(lin_coef_df, aes(x = species, y = Value, fill = species)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(aes(label = label, y = Value/2),
              position = position_dodge(width = 0.9),
              vjust = 0.5, color = "black", size = 4) +
    scale_fill_manual(values = cb_palette, name = "") +
    facet_wrap(~ Coefficient, scales = "free_y") +
    labs(title = "linear models",
         # subtitle = expression(italic(NDVI[Quantiles]) == a + b * x + epsilon),
         subtitle = expression(NDVI == a + b * x + epsilon),
         x = "",
         y = "coefficient value") +
    theme_minimal() +
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
    ) +
    labs(caption = "(a)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0.5))
  
  ##############################################
  # Coefficient Extraction and Plotting - Exponential
  ##############################################
  
  exp_coef_list <- lapply(exp_species, function(sp) {
    mod <- models_exp[[sp]]
    if (is.null(mod)) return(NULL)
    summ <- summary(mod)$coefficients
    df <- data.frame(
      species = sp,
      Coefficient = rownames(summ),
      Value = summ[, "Estimate"],
      pvalue = summ[, "Pr(>|t|)"],
      stringsAsFactors = FALSE
    )
    # Convert coefficient names to lowercase (a, b, c)
    df$Coefficient <- tolower(df$Coefficient)
    return(df)
  })
  # Remove any NULL entries:
  exp_coef_list <- exp_coef_list[!sapply(exp_coef_list, is.null)]
  exp_coef_df <- do.call(rbind, exp_coef_list)
  exp_coef_df$species <- factor(exp_coef_df$species, levels = exp_species)
  exp_coef_df <- exp_coef_df %>% 
    mutate(label = if_else(pvalue < 0.05, "*", sprintf("%.2f", pvalue)))
  
  p_coeff_exp <- ggplot(exp_coef_df, aes(x = species, y = Value, fill = species)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(aes(label = label, y = Value/2),
              position = position_dodge(width = 0.9),
              vjust = 0.5, color = "black", size = 4) +
    scale_fill_manual(values = cb_palette, name = "") +
    facet_wrap(~ Coefficient, scales = "free_y") +
    labs(title = "exponential models",
         # subtitle = expression(italic(NDVI[Quantiles]) == a + b * e^{-c * italic(x)} + epsilon),
         subtitle = expression(NDVI == a + b * e^{c * italic(x)} + epsilon),
         x = "",
         y = "") +
    theme_minimal() +
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
    ) +
    labs(caption = "(b)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0.5))
  
  # Combine the coefficient plots into one figure with panel labels (a) and (b)
  combined_coeff <- p_coeff_linear + p_coeff_exp 
  print(combined_coeff)
  
  # Save the combined coefficient plot
  ggsave(filename = combined_coef_fig, plot = combined_coeff, device = "png", width = 10, height = 8, dpi = 300)
  
  #################################################
  # Combined Final Plot with Model Predictions (Panel A)
  #################################################
  
  # Generate predictions for linear species (Oak & Beech)
  pred_linear <- lapply(linear_species, function(sp) {
    sp_data <- data_clean %>% filter(species == sp)
    x_seq <- seq(min(sp_data$x, na.rm = TRUE), max(sp_data$x, na.rm = TRUE), length.out = 100)
    lm_mod <- models_linear[[sp]]
    pred <- predict(lm_mod, newdata = data.frame(x = x_seq))
    data.frame(species = sp, x = x_seq, pred = pred)
  })
  pred_linear_df <- do.call(rbind, pred_linear)
  
  # Generate predictions for exponential species (Spruce & Pine)
  pred_exp <- lapply(exp_species, function(sp) {
    sp_data <- data_clean %>% filter(species == sp)
    x_seq <- seq(min(sp_data$x, na.rm = TRUE), max(sp_data$x, na.rm = TRUE), length.out = 100)
    nls_mod <- models_exp[[sp]]
    if(is.null(nls_mod)) return(NULL)
    pred <- predict(nls_mod, newdata = data.frame(x = x_seq))
    data.frame(species = sp, x = x_seq, pred = pred)
  })
  pred_exp_df <- do.call(rbind, pred_exp)
  
  # Combine all predictions for Panel A
  pred_all <- bind_rows(pred_linear_df, pred_exp_df)
  
  p_combined <- ggplot() +
    geom_point(data = data_clean, aes(x = x, y = avg_value, color = species)) +
    geom_line(data = pred_all, aes(x = x, y = pred, color = species), linewidth = 1) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "black", linewidth = 1) +
    annotate("text", x = min(data_clean$x, na.rm = TRUE), y = threshold,
             label = "median", hjust = -0.1, vjust = -0.3, fontface = "italic", size = 6) +
    scale_color_manual(values = cb_palette, name = "") +
    labs(x = "soil water potential (kPa)", y = "NDVI quantiles (rank)") +
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
  
  #################################################
  # Additional Analysis for Panel B & C: x50 and Slope for All Species
  #################################################
  
  # For linear models: calculate x50 and absolute slope (and p-value from slope)
  lin_stats <- lapply(linear_species, function(sp) {
    mod <- models_linear[[sp]]
    if(is.null(mod)) return(NULL)
    summ <- summary(mod)$coefficients
    intercept <- summ["(Intercept)", "Estimate"]
    slope <- summ["x", "Estimate"]
    se_slope <- summ["x", "Std. Error"]
    p_val_lin <- summ["x", "Pr(>|t|)"]
    x50_lin <- (threshold - intercept) / slope
    abs_slope_lin <- abs(slope)
    r2_lin <- summary(mod)$r.squared
    data.frame(species = sp, model = "linear", x50 = x50_lin, abs_slope = abs_slope_lin,
               r_squared = r2_lin, se = se_slope, p_val = p_val_lin)
  })
  
  # For exponential models: calculate x50 and absolute slope with deltaMethod and compute p-value
  exp_stats <- lapply(exp_species, function(sp) {
    mod <- models_exp[[sp]]
    if(is.null(mod)) return(NULL)
    coefs <- coef(mod)
    a_val <- coefs["a"]
    b_val <- coefs["b"]
    c_val <- coefs["c"]
    x50_exp <- ifelse((threshold - a_val) > 0 & b_val > 0,
                      log((threshold - a_val)/b_val) / c_val,
                      NA)
    slope50_exp <- c_val * (threshold - a_val)
    abs_slope_exp <- abs(slope50_exp)
    sp_data <- data_clean %>% filter(species == sp)
    fitted_vals <- predict(mod, newdata = sp_data)
    r2_exp <- 1 - sum((sp_data[[value_col]] - fitted_vals)^2) / 
      sum((sp_data[[value_col]] - mean(sp_data[[value_col]]))^2)
    dm_result <- deltaMethod(mod, paste0("c*(a - ", threshold, ")"),
                             parameterNames = c("a", "b", "c"))
    se_exp <- dm_result$SE
    df_exp <- summary(mod)$df[2]
    t_val_exp <- slope50_exp / se_exp
    p_val_exp <- 2 * (1 - pt(abs(t_val_exp), df_exp))
    data.frame(species = sp, model = "exponential", x50 = x50_exp, abs_slope = abs_slope_exp,
               r_squared = r2_exp, se = se_exp, p_val = p_val_exp)
  })
  
  stats_all <- rbind(do.call(rbind, lin_stats), do.call(rbind, exp_stats))
  stats_all$species <- factor(stats_all$species, levels = c(linear_species, exp_species))
  
  # Panel B: Bar plot of x50 (transpiration deficit) for all species
  p_x50 <- ggplot(stats_all, aes(x = species, y = x50, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    # geom_text(aes(label = sprintf("%.2f", x50), y = x50/2),
    #          color = "black", size = 4) +
    scale_fill_manual(values = cb_palette) +
    labs(x = "", y = "soil water potential (kPa)") +
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
    labs(caption = "(b)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0.5))
  
  # Panel C: Bar plot of absolute slope for all species with error bars and R² annotations.
  # The annotation includes an asterisk if p < 0.05.
  p_slope <- ggplot(stats_all, aes(x = species, y = abs_slope, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, abs_slope - se), ymax = abs_slope + se), width = 0.2) +
    geom_text(aes(label = if_else(p_val < 0.05,
                                  sprintf("%.2f*", r_squared),
                                  sprintf("%.2f", r_squared)),
                  y = abs_slope/2),
              color = "black", size = 5) +
    scale_fill_manual(values = cb_palette) +
    labs(x = "", y = "absolute slope") +
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
    labs(caption = "(c)") +
    theme(plot.caption = element_text(face = "bold", size = 16, hjust = 0.5))
  
  # Combine panels: Panel (a) on the left; Panels (b) and (c) in the right column
  final_plot <- (p_combined + (p_x50 / p_slope)) +
    plot_layout(widths = c(2, 1), guides = "collect") &
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 14),
      legend.key = element_rect(fill = "white", color = NA),
      legend.background = element_blank(),         # <- No background
      legend.box.background = element_blank()      # <- No box border
    )
  print(final_plot)
  
  #### Plateau point calculation for exponential models (Spruce & Pine only)
  plateau_df <- lapply(exp_species, function(sp) {
    mod <- models_exp[[sp]]
    if (is.null(mod)) return(NULL)
    coefs <- coef(mod)
    c_val <- coefs["c"]
    x_plateau <- -log(0.05) / c_val
    data.frame(species = sp, x_plateau = x_plateau)
  }) %>% bind_rows()
  
  print("Estimated soil water potential (x_plateau) where NDVI is ~95% of its asymptotic minimum:")
  print(plateau_df)
  
  # Save the combined final plot
  ggsave(filename = output_figure, plot = final_plot, device = "png", width = 10, height = 8, dpi = 300)
}

plot_combined_AIC_R2(df,  "results/key_displays_July_August/AIC.png")

plot_time_series_and_correlation_combined(df,
                                          "results/key_displays_July_August/time_series_Quantiles_PSI_TDiff_species.png")

plot_yearly_mean_linear_coeffs(df, "results/key_displays_July_August/yearly_mean_linear_coeffs.png")

plot_Quantiles_PSI_exp_linear_slope_coeff(df, 
                                          "results/key_displays_July_August/NDVI_Q_PSIbin_exp_linear_coeff_negPSI.png",
                                          "results/key_displays_July_August/NDVI_Q_PSIbin_exp_linear_slope_negPSI.png")

plot_Quantiles_TDiff_exp_linear_slope_coeff(df,
                                            "results/key_displays_July_August/NDVI_Q_TDiffbin_exp_linear_coeff.png",
                                            "results/key_displays_July_August/NDVI_Q_TDiffbin_exp_linear_slope.png")

plot_Quantiles_TDiff_exp_slope_coeff(df, 
                                     "results/key_displays_July_August/NDVI_Q_TDiffbin_exp_coeff.png",
                                     "results/key_displays_July_August/NDVI_Q_TDiffbin_exp_slope.png")

plot_TDiff_PSIbin_poly_2_slope(df, 
                               "results/key_displays_July_August/TDiff_PSIbin_poly_2_coeff.png",
                               "results/key_displays_July_August/TDiff_PSIbin_poly_2_slope.png")

plot_TDiff_PSIbin_poly_3_slope(df, 
                               "results/key_displays_July_August/TDiff_PSIbin_poly_3_coeff.png",
                               "results/key_displays_July_August/TDiff_PSIbin_poly_3_slope.png")

plot_density_Quantiles_PSI_linear(df, 
                                  output_path = "results/key_displays_July_August/Quantiles_PSI_density.png")

plot_density_Quantiles_TDiff_linear(df, 
                                    output_path = "results/key_displays_July_August/Quantiles_TDiff_density.png")

##### Change the size of point based on percentage of pixels of each bin

plot_Quantiles_PSI_exp_linear_slope_coeff <- function(data, combined_coef_fig, output_figure) {
  # Process data
  data <- NDVI_PSIbin(data)
  data <- na.omit(data)
  
  # Load libraries
  library(ggplot2); library(dplyr); library(tidyr); library(tibble)
  library(patchwork); library(purrr); library(nlme); library(car); library(broom)
  
  # Define species and palette
  value_col      <- "avg_value"
  linear_species <- c("Oak", "Beech")
  exp_species    <- c("Spruce", "Pine")
  data$species   <- factor(data$species, levels = c(linear_species, exp_species))
  cb_palette     <- c(
    "Oak"    = "#E69F00",
    "Beech"  = "#0072B2",
    "Spruce" = "#009E73",
    "Pine"   = "#F0E442"
  )
  
  # Prepare data
  data <- data %>% mutate(x = bin_median)
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  threshold <- 11.5
  
  ##########################
  # Model Fitting by Group #
  ##########################
  models_linear <- list()
  for (sp in linear_species) {
    sp_data <- data_clean %>% filter(species == sp)
    models_linear[[sp]] <- lm(avg_value ~ x, data = sp_data)
  }
  models_exp <- list()
  start_list <- list(a = 5, b = 3, c = 0.001)
  control_params <- nls.control(maxiter = 1200, minFactor = 1e-9)
  for (sp in exp_species) {
    sp_data <- data_clean %>% filter(species == sp)
    if (nrow(sp_data) < 5) { models_exp[[sp]] <- NULL; next }
    models_exp[[sp]] <- tryCatch({
      nls(avg_value ~ a + b * exp(c * x), data = sp_data,
          start = start_list, control = control_params)
    }, error = function(e) NULL)
  }
  
  ##############################################
  # Predictions for Panel A
  ##############################################
  pred_linear <- bind_rows(lapply(linear_species, function(sp) {
    df <- data_clean %>% filter(species == sp)
    x_seq <- seq(min(df$x), max(df$x), length.out = 100)
    data.frame(species = sp, x = x_seq,
               pred = predict(models_linear[[sp]], newdata = data.frame(x = x_seq)))
  }))
  pred_exp <- bind_rows(lapply(exp_species, function(sp) {
    mod <- models_exp[[sp]]; if (is.null(mod)) return(NULL)
    df <- data_clean %>% filter(species == sp)
    x_seq <- seq(min(df$x), max(df$x), length.out = 100)
    data.frame(species = sp, x = x_seq,
               pred = predict(mod, newdata = data.frame(x = x_seq)))
  }))
  pred_all <- bind_rows(pred_linear, pred_exp)
  
  # Panel A: combined
  p_combined <- ggplot() +
    geom_point(data = data_clean,
               aes(x = x, y = avg_value, color = species, size = percentage), alpha = 0.7) +
    geom_line(data = pred_all,
              aes(x = x, y = pred, color = species), linewidth = 1) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "black", linewidth = 1) +
    annotate("text", x = min(data_clean$x), y = threshold,
             label = "median", fontface = "italic",
             hjust = -0.1, vjust = -0.3, size = 5) +
    scale_color_manual(values = cb_palette, name = "Species") +
    scale_size_continuous(
      name = "% pixels per bin",
      range = c(1, 8)
    ) +
    labs(x = "soil water potential (kPa)",
         y = "NDVI quantiles (rank)") +
    ggtitle("(a)") +
    theme_minimal() +
    theme(
      plot.title        = element_text(hjust = 0, vjust = 1, size = 16, face = "bold"),
      axis.title        = element_text(face = "bold", size = 16),
      axis.text         = element_text(color = "black", size = 14),
      axis.text.x       = element_text(angle = 0, hjust = 0.5),
      plot.background   = element_rect(fill = "white", color = "white"),
      panel.background  = element_rect(fill = "white"),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      legend.background = element_rect(fill = "white", color = "white"),
      legend.text       = element_text(size = 14),
      legend.position   = "bottom"
    )
  
  #################################################
  # Stats for Panels B & C
  #################################################
  lin_stats <- bind_rows(lapply(linear_species, function(sp) {
    mod <- models_linear[[sp]]; summ <- summary(mod)$coefficients
    intercept <- summ["(Intercept)", "Estimate"]; slope <- summ["x", "Estimate"]
    se_slope <- summ["x", "Std. Error"]; p_val <- summ["x", "Pr(>|t|)"]
    x50 <- (threshold - intercept) / slope; abs_sl <- abs(slope)
    r2 <- summary(mod)$r.squared
    data.frame(species = sp, x50 = x50, abs_slope = abs_sl,
               se = se_slope, p_val = p_val, r_squared = r2)
  }))
  exp_stats <- bind_rows(lapply(exp_species, function(sp) {
    mod <- models_exp[[sp]]; if (is.null(mod)) return(NULL)
    co <- coef(mod); a <- co["a"]; b <- co["b"]; c <- co["c"]
    x50 <- ifelse((threshold - a) > 0 & b > 0, log((threshold - a)/b)/c, NA)
    slope50 <- c * (threshold - a); abs_sl <- abs(slope50)
    df_sp <- data_clean %>% filter(species == sp)
    r2 <- 1 - sum((df_sp[[value_col]] - predict(mod, df_sp))^2) /
      sum((df_sp[[value_col]] - mean(df_sp[[value_col]]))^2)
    dm <- deltaMethod(mod, paste0("c*(a - ", threshold, ")"), parameterNames = c("a","b","c"))
    se <- dm$SE; df_res <- summary(mod)$df[2]
    t_val <- slope50 / se; p_val <- 2*(1-pt(abs(t_val), df_res))
    data.frame(species = sp, x50 = x50, abs_slope = abs_sl,
               se = se, p_val = p_val, r_squared = r2)
  }))
  stats_all <- bind_rows(lin_stats, exp_stats)
  stats_all$species <- factor(stats_all$species, levels = c(linear_species, exp_species))
  
  # Panel B: x50
  p_x50 <- ggplot(stats_all, aes(x = species, y = x50, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = cb_palette, guide = FALSE) +
    labs(x = "", y = "soil water potential (kPa)") +
    ggtitle("(b)") +
    theme_minimal() +
    theme(
      plot.title        = element_text(hjust = 0, vjust = 1, size = 16, face = "bold"),
      axis.title        = element_text(face = "bold", size = 16),
      axis.text         = element_text(color = "black", size = 14),
      axis.text.x       = element_text(angle = 0, hjust = 0.5),
      panel.background  = element_rect(fill = "white"),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank()
    )
  
  # Panel C: absolute slope
  p_slope <- ggplot(stats_all, aes(x = species, y = abs_slope, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, abs_slope - se), ymax = abs_slope + se), width = 0.2) +
    geom_text(aes(label = if_else(p_val < 0.05,
                                  sprintf("%.2f*", r_squared),
                                  sprintf("%.2f", r_squared)),
                  y = abs_slope/2), size = 5) +
    scale_fill_manual(values = cb_palette, guide = FALSE) +
    labs(x = "", y = "absolute slope") +
    ggtitle("(c)") +
    theme_minimal() +
    theme(
      plot.title        = element_text(hjust = 0, vjust = 1, size = 16, face = "bold"),
      axis.title        = element_text(face = "bold", size = 16),
      axis.text         = element_text(color = "black", size = 14),
      axis.text.x       = element_text(angle = 0, hjust = 0.5),
      panel.background  = element_rect(fill = "white"),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank()
    )
  
  # Combine panels
  final_plot <- (p_combined + (p_x50 / p_slope)) +
    plot_layout(widths = c(2, 1), guides = "collect") &
    theme(
      legend.position       = "bottom",
      legend.title          = element_blank(),
      legend.text           = element_text(size = 14),
      legend.key            = element_rect(fill = "white", color = NA),
      legend.background     = element_blank(),
      legend.box.background = element_blank()
    )
  
  print(final_plot)
  ggsave(output_figure, final_plot, device = "png", width = 10, height = 8, dpi = 300)
}

plot_Quantiles_PSI_exp_linear_slope_coeff(df, 
                                          "results/key_displays_July_August/NDVI_Q_PSIbin_exp_linear_coeff_negPSI_size.png",
                                          "results/key_displays_July_August/NDVI_Q_PSIbin_exp_linear_slope_negPSI_size.png")

plot_Quantiles_TDiff_exp_linear_slope_coeff <- function(data, combined_coef_fig, output_figure) {
  # Process data
  data <- NDVI_TDiffbin(data)
  data <- na.omit(data)
  
  # Load libraries
  library(ggplot2); library(dplyr); library(tidyr); library(tibble)
  library(patchwork); library(purrr); library(nlme); library(car); library(broom)
  
  # Define species and palette
  value_col      <- "avg_value"
  linear_species <- c("Oak", "Beech")
  exp_species    <- c("Spruce", "Pine")
  data$species   <- factor(data$species, levels = c(linear_species, exp_species))
  cb_palette     <- c(
    "Oak"    = "#E69F00",
    "Beech"  = "#0072B2",
    "Spruce" = "#009E73",
    "Pine"   = "#F0E442"
  )
  
  # Prepare data
  data <- data %>% mutate(x = bin_median)
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  threshold <- 11.5
  
  ##########################
  # Model Fitting by Group #
  ##########################
  models_linear <- list()
  for (sp in linear_species) {
    sp_data <- data_clean %>% filter(species == sp)
    models_linear[[sp]] <- lm(avg_value ~ x, data = sp_data)
  }
  models_exp <- list()
  start_list <- list(a = 5, b = 7, c = 0.04)
  control_params <- nls.control(maxiter = 1200, minFactor = 1e-9)
  for (sp in exp_species) {
    sp_data <- data_clean %>% filter(species == sp)
    if (nrow(sp_data) < 5) { models_exp[[sp]] <- NULL; next }
    models_exp[[sp]] <- tryCatch({
      nls(avg_value ~ a + b * exp(-c * x), data = sp_data,
          start = start_list, control = control_params)
    }, error = function(e) NULL)
  }
  
  ##############################################
  # Predictions for Panel A                    #
  ##############################################
  pred_linear <- bind_rows(lapply(linear_species, function(sp) {
    df   <- data_clean %>% filter(species == sp)
    xseq <- seq(min(df$x), max(df$x), length.out = 100)
    data.frame(species = sp, x = xseq,
               pred = predict(models_linear[[sp]], newdata = data.frame(x = xseq)))
  }))
  pred_exp <- bind_rows(lapply(exp_species, function(sp) {
    mod <- models_exp[[sp]]; if (is.null(mod)) return(NULL)
    df   <- data_clean %>% filter(species == sp)
    xseq <- seq(min(df$x), max(df$x), length.out = 100)
    data.frame(species = sp, x = xseq,
               pred = predict(mod, newdata = data.frame(x = xseq)))
  }))
  pred_all <- bind_rows(pred_linear, pred_exp)
  
  # Panel A: combined with size aesthetic
  p_combined <- ggplot() +
    geom_point(data = data_clean,
               aes(x = x, y = avg_value, color = species, size = percentage), alpha = 0.7) +
    geom_line(data = pred_all,
              aes(x = x, y = pred, color = species), linewidth = 1) +
    geom_hline(yintercept = threshold, linetype = "dashed",
               color = "black", linewidth = 1) +
    annotate("text", x = 35, y = threshold,
             label = "median", fontface = "italic",
             hjust = -0.1, vjust = -0.3, size = 5) +
    scale_color_manual(values = cb_palette, name = "Species") +
    scale_size_continuous(name = "% pixels per bin", range = c(1, 8)) +
    labs(x = "transpiration deficit", y = "NDVI quantiles") +
    ggtitle("(a)") +
    theme_minimal() +
    theme(
      plot.title        = element_text(hjust = 0, vjust = 1, size = 16, face = "bold"),
      axis.title        = element_text(face = "bold", size = 16),
      axis.text         = element_text(color = "black", size = 14),
      axis.text.x       = element_text(angle = 0, hjust = 0.5),
      panel.background  = element_rect(fill = "white"),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      legend.position   = "bottom",
      legend.text       = element_text(size = 14),
      legend.background = element_rect(fill = "white", color = "white")
    )
  
  #################################################
  # Stats for Panels B & C                         #
  #################################################
  lin_stats <- bind_rows(lapply(linear_species, function(sp) {
    mod     <- models_linear[[sp]]; summ <- summary(mod)$coefficients
    a       <- summ["(Intercept)", "Estimate"]; b <- summ["x", "Estimate"]
    se      <- summ["x", "Std. Error"]; pval <- summ["x", "Pr(>|t|)"]
    x50     <- (threshold - a) / b; abs_s  <- abs(b)
    r2      <- summary(mod)$r.squared
    data.frame(species = sp, x50 = x50, abs_slope = abs_s,
               se = se, p_val = pval, r_squared = r2)
  }))
  exp_stats <- bind_rows(lapply(exp_species, function(sp) {
    mod  <- models_exp[[sp]]; if (is.null(mod)) return(NULL)
    co   <- coef(mod); a <- co["a"]; b <- co["b"]; c <- co["c"]
    x50  <- ifelse((threshold - a) > 0 & b > 0,
                   -log((threshold - a)/b) / c, NA)
    slope50 <- -c * (threshold - a); abs_s    <- abs(slope50)
    df_sp    <- data_clean %>% filter(species == sp)
    r2       <- 1 - sum((df_sp[[value_col]] - predict(mod, df_sp))^2) /
      sum((df_sp[[value_col]] - mean(df_sp[[value_col]]))^2)
    dm       <- deltaMethod(mod, paste0("c*(a - ", threshold, ")"),
                            parameterNames = c("a","b","c"))
    se       <- dm$SE; df_res <- summary(mod)$df[2]
    tval     <- slope50 / se; pval <- 2 * (1 - pt(abs(tval), df_res))
    data.frame(species = sp, x50 = x50, abs_slope = abs_s,
               se = se, p_val = pval, r_squared = r2)
  }))
  stats_all <- bind_rows(lin_stats, exp_stats)
  stats_all$species <- factor(stats_all$species, levels = c(linear_species, exp_species))
  
  # Panel B: x50
  p_x50 <- ggplot(stats_all, aes(x = species, y = x50, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = cb_palette, guide = FALSE) +
    labs(x = "", y = "transpiration deficit") +
    ggtitle("(b)") +
    theme_minimal() +
    theme(
      plot.title       = element_text(hjust = 0, vjust = 1, size = 16, face = "bold"),
      axis.title       = element_text(face = "bold", size = 16),
      axis.text        = element_text(color = "black", size = 14),
      axis.text.x      = element_text(angle = 0, hjust = 0.5),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  # Panel C: absolute slope
  p_slope <- ggplot(stats_all, aes(x = species, y = abs_slope, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, abs_slope - se), ymax = abs_slope + se), width = 0.2) +
    geom_text(aes(label = if_else(p_val < 0.05,
                                  sprintf("%.2f*", r_squared),
                                  sprintf("%.2f", r_squared)), y = abs_slope / 2), size = 5) +
    scale_fill_manual(values = cb_palette, guide = FALSE) +
    labs(x = "", y = "absolute slope") +
    ggtitle("(c)") +
    theme_minimal() +
    theme(
      plot.title       = element_text(hjust = 0, vjust = 1, size = 16, face = "bold"),
      axis.title       = element_text(face = "bold", size = 16),
      axis.text        = element_text(color = "black", size = 14),
      axis.text.x      = element_text(angle = 0, hjust = 0.5),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  # Combine panels with unified bottom legends
  final_plot <- (p_combined + (p_x50 / p_slope)) +
    plot_layout(widths = c(2, 1), guides = "collect") &
    theme(
      legend.position       = "bottom",
      legend.title          = element_blank(),
      legend.text           = element_text(size = 14),
      legend.key            = element_rect(fill = "white", color = NA),
      legend.background     = element_blank(),
      legend.box.background = element_blank()
    )
  
  print(final_plot)
  ggsave(output_figure, final_plot, device = "png",
         width = 10, height = 8, dpi = 300)
}

plot_Quantiles_TDiff_exp_linear_slope_coeff(df,
                                            "results/key_displays_July_August/NDVI_Q_TDiffbin_exp_linear_coeff_size.png",
                                            "results/key_displays_July_August/NDVI_Q_TDiffbin_exp_linear_slope_size.png")

plot_TDiff_PSIbin_poly_2_slope <- function(data, coef_output, figure_output) {
  # Load required libraries
  library(lme4)      # For mixed-effects modeling
  library(dplyr)     # For data manipulation
  library(ggplot2)   # For plotting
  library(patchwork) # For combining ggplots
  library(tidyr)     # For pivoting data
  library(car)       # For deltaMethod
  library(nlme)      # For nlsList (if needed)
  
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
  
  ## Panel A: Mixed-Effects Model Plot (same as before)
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
  
  # Create the plot
  plot_mixed <- ggplot(TDiff_PSIbin_df, aes(x = bin_median, y = avg_transpiration_deficit, color = species)) +
    geom_point(data = TDiff_PSIbin_df,
               aes(x = bin_median,
                   y = avg_transpiration_deficit,
                   color = species,
                   size = percentage)) +
    geom_line(data = pred_data,
              aes(x = bin_median,
                  y = predicted,
                  color = species),
                  linewidth = 1) +
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
  
  ## Analytical Calculation of x50 and Slope (modified with delta method)
  
  # Extract species-specific coefficients
  species_coef <- coef(model)$species %>%
    tibble::rownames_to_column("species") %>%
    rename(a = `(Intercept)`, 
           b = `poly(bin_median, 2, raw = TRUE)1`, 
           c = `poly(bin_median, 2, raw = TRUE)2`)
  
  # For each species, fit a separate model to use with deltaMethod
  species_list <- species_levels
  stats_list <- lapply(species_list, function(sp) {
    # Fit model for this species only
    df_sp <- TDiff_PSIbin_df %>% filter(species == sp)
    mod_sp <- lm(avg_transpiration_deficit ~ poly(bin_median, 2, raw = TRUE), data = df_sp)
    
    # Get coefficients for this species
    coefs <- coef(mod_sp)
    a <- coefs[1]; b <- coefs[2]; c <- coefs[3]
    
    # Calculate x50 (where curve crosses mean transpiration deficit)
    disc <- b^2 - 4 * c * (a - line_val)
    x50 <- ifelse(disc >= 0, (-b + sqrt(disc)) / (2 * c), NA_real_)
    
    # Use delta method to get SE for slope at x50 (derivative: b + 2*c*x)
    if(!is.na(x50)) {
      # Expression for slope at x50
      slope_expr <- paste0("b + 2 * c * ", x50)
      dm_result <- deltaMethod(mod_sp, slope_expr, parameterNames = c("a", "b", "c"))
      
      # Calculate p-value
      t_val <- (b + 2 * c * x50) / dm_result$SE
      df_resid <- df.residual(mod_sp)
      p_val <- 2 * (1 - pt(abs(t_val), df_resid))
      
      # Calculate R-squared
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
    scale_fill_manual(values = cb_palette, name ="") +
    guides(fill = "none") +
    labs(x = "", y = "soil water potential (kPa)") +
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
  
  ## Panel C: Bar Plot of Absolute Slopes with Error Bars (using delta method SE)
  p_bar <- ggplot(stats_df, aes(x = species, y = slope_abs, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, slope_abs - se), 
                      ymax = slope_abs + se), 
                  width = 0.2) +
    geom_text(aes(label = label_text, y = slope_abs/2), 
              color = "black", size = 5) +
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
  
  # Combine plots
  final_plot <- (plot_mixed + (p_bar_x / p_bar)) +
    plot_layout(widths = c(2, 1), guides = "collect") &
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 14),
      legend.key = element_rect(fill = "white", color = NA),
      legend.background = element_blank(),
      legend.box.background = element_blank()
    )
  
  print(final_plot)
  ggsave(figure_output, plot = final_plot, width = 10, height = 8, dpi = 300)
  
  ## Coefficient Plot (same as before)
  # Extract the species-specific (conditional) coefficients.
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
  species_list <- species_levels
  coeff_stats_list <- lapply(species_list, function(sp) {
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
  
  return(list(combined_plot = final_plot, stats_df = stats_df))
}

plot_TDiff_PSIbin_poly_2_slope(df, 
                               "results/key_displays_July_August/TDiff_PSIbin_poly_2_coeff_size.png",
                               "results/key_displays_July_August/TDiff_PSIbin_poly_2_slope_size.png")
