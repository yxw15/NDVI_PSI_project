setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")
load("results/Data/AllSpecies_AllMonths_rootzone.RData")

library(tidyverse)
library(patchwork)
library(terra)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)
library(tidyr)
library(grid) 

data_all <- combined %>% filter(Quantiles > 0) %>% filter(month %in% c("May", "June", "July", "August"))
data_all <- na.omit(data_all)
data <- data_all %>% filter(month %in% c("July", "August"))

### S1 original NDVI - PSI ###
plot_NDVI_PSI_month_species_linear_orig <- function(data, figure_output = NULL, figure_output2 = NULL) {
  # Load required libraries
  library(ggplot2)
  library(dplyr)
  library(viridis)  # For color-blind-friendly colors
  library(tidyr)
  library(tidyverse)
  
  data <- data %>% filter(Quantiles > 0)
  
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
      x = "soil water potential (kPa)",
      y = "NDVI quantiles (rank)",
      color = "",
      title = "NDVI quantiles (rank) ~ soil water potential (kPa)",
      subtitle = expression(NDVI == a + bx + epsilon),
    ) +
    ylim(0,22) +
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
      title = "NDVI quantiles (rank) ~ soil water potential (kPa)",
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

plot_NDVI_PSI_month_species_linear_orig(data_all, 
                                        "results_rootzone/Figures/supplementary/NDVI_PSI_linear_slope_monthly_orig.png",
                                        "results_rootzone/Figures/supplementary/NDVI_PSI_linear_coeff_monthly_orig.png")

### S2 original NDVI - TDiff ###
plot_NDVI_TDiff_month_species_linear_orig <- function(data, figure_output = NULL, figure_output2 = NULL) {
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
      x = "transpiration deficit (mm)",
      y = "NDVI quantiles (rank)",
      color = "",
      title = "NDVI quantiles (rank) ~ transpiration deficit (mm)",
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
      title = "NDVI quantile (rank) ~ transpiration deficit (mm)",
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

plot_NDVI_TDiff_month_species_linear_orig(data_all, 
                                          "results_rootzone/Figures/supplementary/NDVI_TDiff_linear_slope_monthly_orig.png",
                                          "results_rootzone/Figures/supplementary/NDVI_TDiff_linear_coeff_monthly_orig.png")

### S1 density NDVI - PSI ###
plot_density_NDVI_PSI_linear <- function(data, swp_bin_width = 50, output_path) {
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
      title = "NDVI quantiles (rank) ~ soil water potential (kPa)",
      x = "soil water potential (kPa)",
      y = "NDVI quantiles (rank)"
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

plot_density_NDVI_PSI_linear(data,
                             swp_bin_width = 50,
                             output_path = "results_rootzone/Figures/supplementary/NDVI_PSI_density.png")

### S2 density NDVI - TDiff ###
plot_density_NDVI_TDiff_linear <- function(data, tdiff_bin_width = 3, output_path) {
  
  library(dplyr)
  library(tidyr)      # Added to provide drop_na() function
  library(ggplot2)
  library(tidyverse)
  
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
      title = "NDVI quantiles (rank) ~ transpiration deficit (mm)",
      x = "transpiration deficit (mm)",
      y = "NDVI quantiles (rank)"
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

plot_density_NDVI_TDiff_linear <- function(data, tdiff_bin_width = 3, output_path) {
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  
  # 1) Clean + restrict to target species
  df <- data %>%
    mutate(
      species = as.character(species),
      Quantiles = suppressWarnings(as.numeric(Quantiles)),
      transpiration_deficit = suppressWarnings(as.numeric(transpiration_deficit))
    ) %>%
    filter(species %in% species_order) %>%
    mutate(species = factor(species, levels = species_order)) %>%
    drop_na(Quantiles, transpiration_deficit)
  
  if (nrow(df) == 0) stop("No rows left after filtering/NA removal.")
  
  # 2) Robust bin edges for transpiration_deficit
  tmin <- min(df$transpiration_deficit, na.rm = TRUE)
  tmax <- max(df$transpiration_deficit, na.rm = TRUE)
  if (!is.finite(tmin) || !is.finite(tmax)) stop("Non-finite transpiration_deficit range.")
  if (tmin == tmax) { tmin <- tmin - tdiff_bin_width/2; tmax <- tmax + tdiff_bin_width/2 }
  
  edges <- seq(
    floor(tmin / tdiff_bin_width) * tdiff_bin_width,
    ceiling(tmax / tdiff_bin_width) * tdiff_bin_width,
    by = tdiff_bin_width
  )
  if (length(edges) < 2) edges <- c(tmin, tmax)
  
  species_binned <- df %>%
    mutate(
      quantile_bin = cut(Quantiles, breaks = seq(0, 22, 1), include.lowest = TRUE, right = FALSE),
      tdiff_bin    = cut(transpiration_deficit, breaks = edges, include.lowest = TRUE, right = FALSE)
    ) %>%
    drop_na(quantile_bin, tdiff_bin)
  
  if (nrow(species_binned) == 0) stop("All rows dropped during binning; check Quantiles range and tdiff_bin_width.")
  
  # 3) Density within each (species, tdiff_bin)
  density_df <- species_binned %>%
    count(species, quantile_bin, tdiff_bin, name = "count") %>%
    group_by(species, tdiff_bin) %>%
    mutate(density_percent = count / sum(count) * 100) %>%
    ungroup()
  
  species_density <- species_binned %>%
    left_join(density_df, by = c("species", "quantile_bin", "tdiff_bin"))
  
  # 4) Legend safety: determine if colorbar should exist
  max_density <- suppressWarnings(max(species_density$density_percent, na.rm = TRUE))
  has_color   <- is.finite(max_density) && max_density > 0
  if (!has_color) max_density <- 1
  density_breaks <- unique(pretty(c(0, max_density), n = 5))
  
  # 5) Only draw smoother where there’s enough variation
  smooth_df <- species_density %>%
    group_by(species) %>%
    filter(n() >= 3,
           dplyr::n_distinct(transpiration_deficit) > 1,
           dplyr::n_distinct(Quantiles) > 1) %>%
    ungroup()
  
  p <- ggplot(
    species_density,
    aes(x = transpiration_deficit, y = Quantiles, color = density_percent)
  ) +
    geom_point(size = 2.5, alpha = 0.9, na.rm = TRUE) +
    # Add smoother only if we have adequate data
    { if (nrow(smooth_df) > 0)
      geom_smooth(
        data = smooth_df, aes(x = transpiration_deficit, y = Quantiles),
        method = "lm", formula = y ~ x, se = TRUE,
        color = "black", linewidth = 1, na.rm = TRUE
      )
      else NULL
    } +
    { if (has_color)
      scale_color_viridis_c(
        name = "density (%)",
        option = "inferno",
        direction = -1,
        trans = "sqrt",
        breaks = density_breaks,
        limits = c(0, max_density),
        na.value = NA
      )
      else
        scale_color_continuous(guide = "none", na.value = NA)
    } +
    labs(
      title = "NDVI quantiles (rank) ~ transpiration deficit (mm)",
      x = "transpiration deficit (mm)",
      y = "NDVI quantiles (rank)"
    ) +
    facet_wrap(~species, ncol = 2, drop = TRUE) +
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
      legend.title = element_text(size = 14),
      legend.key.width = grid::unit(1.8, "cm"),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    )
  
  print(p)
  ggsave(output_path, plot = p, width = 12, height = 8, dpi = 300)
}

plot_density_NDVI_TDiff_linear(data,
                               tdiff_bin_width = 3,
                               output_path = "results_rootzone/Figures/supplementary/NDVI_Tdiff_density.png")

### S3 S4 distribution of mean PSI and TDiff ###
plot_distribution_PSI_TDiff <- function(df,var_name,legend_title, output_file) {
  
  # load libraries once at top
  library(dplyr)
  library(ggplot2)
  library(grid)          # for unit()
  library(sf)            # for spatial features
  library(rnaturalearth)
  library(rnaturalearthdata)
  
  # cache Germany border once
  germany_border <- ne_countries(
    scale       = "medium",
    country     = "Germany",
    returnclass = "sf"
  )
  
  # input checks
  if (!is.data.frame(df)) stop("`df` must be a data frame.")
  req_cols <- c("x","y","species", var_name)
  if (!all(req_cols %in% names(df))) {
    stop("`df` must contain columns: ", paste(req_cols, collapse = ", "))
  }
  
  # compute per‐species mean
  df_mean <- df %>%
    group_by(x, y, species) %>%
    summarise(mean_val = mean(.data[[var_name]], na.rm = TRUE), .groups = "drop")
  
  # enforce species order
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  df_mean$species <- factor(df_mean$species, levels = species_order)
  
  # choose palette
  base_pal <- c("blue","dodgerblue","cyan","yellow","orange","red")
  pal <- if (var_name == "soil_water_potential") rev(base_pal) else base_pal
  
  # build plot
  p <- ggplot(df_mean, aes(x = x, y = y, color = mean_val)) +
    geom_point(size = 0.5) +
    geom_sf(
      data        = germany_border,
      fill        = NA,
      color       = "black",
      inherit.aes = FALSE
    ) +
    scale_color_gradientn(
      colours = pal,
      name    = legend_title
    ) +
    facet_wrap(~ species, nrow = 1) +
    coord_sf(expand = FALSE) +
    labs(x = "longitude", y = "latitude") +
    guides(
      color = guide_colorbar(
        title.position = "top",
        title.hjust    = 0.5,
        barwidth       = unit(8, "cm"),
        barheight      = unit(0.4, "cm"),
        ticks          = TRUE
      )
    ) +
    theme_minimal() +
    theme(
      axis.text.x        = element_text(hjust = 0.5),
      plot.background    = element_rect(fill = "white", color = NA),
      panel.background   = element_rect(fill = "white"),
      legend.background  = element_rect(fill = "white", color = NA),
      legend.position    = "top",
      legend.text        = element_text(size = 14),
      legend.title       = element_text(size = 14, face = "bold"),
      axis.title         = element_text(face = "bold", size = 16),
      axis.text          = element_text(color = "black", size = 14),
      panel.border       = element_rect(color = "black", fill = NA),
      panel.grid         = element_blank(),
      plot.title         = element_text(hjust = 0.5, size = 18, face = "bold"),
      plot.subtitle      = element_text(hjust = 0.5, size = 14),
      strip.background   = element_rect(fill = "white", color = "black"),
      strip.text         = element_text(face = "bold", size = 12)
    )
  
  # save the plot to file
  ggsave(filename = output_file, plot = p, width = 14, height = 6, dpi = 300)
  # print the plot to the current graphics device
  print(p)
  # also return the plot object
  invisible(p)
}

plot_distribution_PSI_TDiff(
  df           = data,
  var_name     = "soil_water_potential",
  legend_title = "soil water potential (kPa)",
  output_file  = "results_rootzone/Figures/supplementary/mean_soil_water_potential.png"
)

plot_distribution_PSI_TDiff(
  df           = data,
  var          = "transpiration_deficit",
  legend_title = "transpiration deficit (mm)",
  output_file  = "results_rootzone/Figures/supplementary/mean_transpiration_deficit.png"
)

### S5 S6 distribution of mean temperature and VPD ###

load("results/Data/mean_Temp_VPD.RData")

plot_mean_climate <- function(data, var, legend_title, output_path) {
  stopifnot(var %in% names(data))   # check var exists
  
  boundary_germany <- ne_countries(
    scale = "medium",
    country = "Germany",
    returnclass = "sf"
  )
  
  
  p <- ggplot(data) +
    geom_point(aes(x = x, y = y, color = .data[[var]]), size = 0.5) +
    geom_sf(data = boundary_germany, fill = NA, color = "black", inherit.aes = FALSE) +
    scale_color_gradientn(
      colours = c("blue", "dodgerblue", "cyan", "yellow", "orange", "red"),
      name    = legend_title
    ) +
    facet_wrap(~species, nrow = 1) +
    coord_sf(expand = FALSE) +
    labs(
      x = "longitude",
      y = "latitude"
    ) +
    guides(
      color = guide_colorbar(
        title.position  = "top",
        title.hjust     = 0.5,
        barwidth        = unit(8, "cm"),
        barheight       = unit(0.4, "cm"),
        ticks           = TRUE
      )
    ) +
    theme_minimal() +
    theme(
      axis.text.x        = element_text(hjust = 0.5),
      plot.background    = element_rect(fill = "white", color = NA),
      panel.background   = element_rect(fill = "white"),
      legend.background  = element_rect(fill = "white", color = NA),
      legend.position    = "top",
      legend.text        = element_text(size = 14),
      legend.title       = element_text(size = 14, face = "bold"),
      axis.title         = element_text(face = "bold", size = 16),
      axis.text          = element_text(color = "black", size = 14),
      panel.border       = element_rect(color = "black", fill = NA),
      panel.grid         = element_blank(),
      plot.title         = element_text(hjust = 0.5, size = 18, face = "bold"),
      plot.subtitle      = element_text(hjust = 0.5, size = 14),
      strip.background   = element_rect(fill = "white", color = "black"),
      strip.text         = element_text(face = "bold", size = 12)
    )
  
  print(p)
  ggsave(filename = output_path, plot = p, width = 14, height = 6, dpi = 300)
}

plot_mean_climate(all_df, var = "temp", legend_title = "temperature (°C)",
                  output_path = "results_rootzone/Figures/supplementary/mean_temperature.png")

plot_mean_climate(all_df, var = "vpd", legend_title = "VPD (kPa)",
                  output_path = "results_rootzone/Figures/supplementary/mean_vpd.png")


### S7 mean bar PSI TDiff ###
plot_bar_mean_PSI_TDiff_rootzone <- function(df, output_file) {
  
  # enforce species order & color palette
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  cb_palette <- c(
    "Oak"    = "#E69F00",
    "Beech"  = "#0072B2",
    "Spruce" = "#009E73",
    "Pine"   = "#F0E442"
  )
  
  # compute mean values per species
  summary_df <- df %>%
    group_by(species) %>%
    summarise(
      Mean_SW_Potential       = mean(soil_water_potential, na.rm = TRUE),
      Mean_Transpiration_Def  = mean(transpiration_deficit, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = c(Mean_SW_Potential, Mean_Transpiration_Def),
      names_to  = "Parameter",
      values_to = "Mean_Value"
    ) %>%
    mutate(
      species   = factor(species, levels = species_order),
      Parameter = recode(Parameter,
                         Mean_SW_Potential       = "soil water potential (kPa)",
                         Mean_Transpiration_Def  = "transpiration deficit (mm)"
      ),
      Parameter = factor(Parameter, levels = c(
        "soil water potential (kPa)",
        "transpiration deficit (mm)"
      ))
    )
  
  # build two-panel bar plot
  p <- ggplot(summary_df, aes(x = species, y = Mean_Value, fill = species)) +
    geom_col(color = NA, width = 0.7) +
    facet_wrap(~ Parameter, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = cb_palette) +
    labs(
      x     = "",
      y     = "mean value",
      title = NULL,
      fill  = NULL
    ) +
    theme_minimal() +
    theme(
      plot.background   = element_rect(fill = "white", color = NA),
      panel.background  = element_rect(fill = "white"),
      panel.border      = element_rect(color = "black", fill = NA, size = 0.5),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      axis.title        = element_text(face = "bold", size = 16),
      axis.text         = element_text(color = "black", size = 14),
      legend.position   = "none",
      strip.background  = element_rect(fill = "white", color = "black", size = 0.5),
      strip.text        = element_text(face = "bold", size = 14)
    ) 
  
  # save & print
  ggsave(filename = output_file, plot = p, width = 10, height = 5, dpi = 300)
  print(p)
  invisible(p)
}

plot_bar_mean_PSI_TDiff_rootzone(data, "results_rootzone/Figures/supplementary/bar_mean_PSI_TDiff_rootzone.png")

### S8 mean bar Temp VPD ###
load("results/Data/mean_Temp_VPD.RData")

plot_bar_mean_Temp_VPD <- function(df, output_file) {
  # packages
  require(dplyr)
  require(tidyr)
  require(ggplot2)
  
  # enforce species order & color palette
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  cb_palette <- c(
    "Oak"    = "#E69F00",
    "Beech"  = "#0072B2",
    "Spruce" = "#009E73",
    "Pine"   = "#F0E442"
  )
  
  # compute mean values per species
  summary_df <- df %>%
    group_by(species) %>%
    summarise(
      Mean_Temp = mean(temp, na.rm = TRUE),
      Mean_VPD  = mean(vpd,  na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = c(Mean_Temp, Mean_VPD),
      names_to  = "Parameter",
      values_to = "Mean_Value"
    ) %>%
    mutate(
      species   = factor(species, levels = species_order),
      Parameter = dplyr::recode(Parameter,
                                Mean_Temp = "temperature (°C)",
                                Mean_VPD  = "VPD (kPa)"),
      Parameter = factor(Parameter, levels = c("temperature (°C)", "VPD (kPa)"))
    )
  
  # build two-panel bar plot
  p <- ggplot(summary_df, aes(x = species, y = Mean_Value, fill = species)) +
    geom_col(color = NA, width = 0.7) +
    facet_wrap(~ Parameter, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = cb_palette) +
    labs(
      x     = "",
      y     = "mean value",
      title = NULL,
      fill  = NULL
    ) +
    theme_minimal() +
    theme(
      plot.background   = element_rect(fill = "white", color = NA),
      panel.background  = element_rect(fill = "white"),
      panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      axis.title        = element_text(face = "bold", size = 16),
      axis.text         = element_text(color = "black", size = 14),
      legend.position   = "none",
      strip.background  = element_rect(fill = "white", color = "black", size = 0.5),
      strip.text        = element_text(face = "bold", size = 14)
    )
  
  # save & print
  ggsave(filename = output_file, plot = p, width = 10, height = 5, dpi = 300)
  print(p)
  invisible(p)
}

plot_bar_mean_Temp_VPD(all_df, "results_rootzone/Figures/supplementary/bar_mean_temp_vpd.png")

# =========================
# S9 time series Temp and VPD (organized)
# =========================

# ====== 0) Libraries ======
suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(ggplot2)
  library(scales)
  library(cowplot)
})

# ====== 1) Helpers ======
es.func  <- function(tmean)  6.11 * exp((2.5e6 / 461) * (1/273 - 1/(273 + tmean))) # hPa
vpd.func <- function(es.hy) { es <- es.hy[1]; hy <- es.hy[2]; ((100 - hy) / 100) * es } # hPa

ensure_crs <- function(r, fallback) {
  if (is.null(r)) return(r)
  cr <- try(crs(r), silent = TRUE)
  if (inherits(cr, "try-error") || is.na(cr) || cr == "") {
    crs(r) <- fallback
  }
  r
}

check_exists <- function(path_vec) {
  missing <- path_vec[!file.exists(path_vec)]
  if (length(missing)) {
    stop(sprintf("Missing required file(s):\n- %s", paste(missing, collapse = "\n- ")), call. = FALSE)
  }
}

# ====== 2) Parameters & I/O ======
species_list <- c("Oak","Beech","Spruce","Pine")
years        <- 2003:2024

paths <- list(
  hu_file      = "../E_OBS_HU/hu_ens_mean_0.1deg_reg_v31.0e.nc",
  mask_pattern = "species_map_MODIS/%s.tif",
  dwd_temp_j   = function(yr) sprintf("../DWD_temp/grids_germany_monthly_air_temp_mean_%d07.asc", yr),
  dwd_temp_a   = function(yr) sprintf("../DWD_temp/grids_germany_monthly_air_temp_mean_%d08.asc", yr),
  rdata_out    = "results/Data/species_Temp_VPD_timeseries.RData",
  csv_out      = "results/Data/species_Temp_VPD_timeseries.csv",
  fig_dir      = "results_rootzone/Figures/supplementary"
)

dwd_crs     <- "epsg:31467"  # DWD ASCII grids (adjust if yours differ)
hu_fallback <- "epsg:4326"   # E-OBS lat/lon, used if missing in file

dir.create("results/Data", showWarnings = FALSE, recursive = TRUE)
dir.create(paths$fig_dir, showWarnings = FALSE, recursive = TRUE)

# Quick upfront existence checks for static inputs
check_exists(paths$hu_file)

# ====== 3) Load humidity cube (once) ======
hu_all   <- rast(paths$hu_file)                 # relative humidity (%)
hu_all   <- ensure_crs(hu_all, hu_fallback)
hu_dates <- as.Date(time(hu_all))

# Precompute RH-layer index/mean per year (Jul 28 & Aug 29)
rh_layers_by_year <- local({
  jul28 <- sprintf("%d-07-28", years)
  aug29 <- sprintf("%d-08-29", years)
  idxJ  <- match(as.Date(jul28), hu_dates)
  idxA  <- match(as.Date(aug29), hu_dates)
  map2(idxJ, idxA, function(iJ, iA) {
    idx <- c(iJ, iA)
    idx <- idx[!is.na(idx)]
    if (!length(idx)) return(NULL)
    if (length(idx) == 1) hu_all[[idx]] else app(hu_all[[idx]], mean)
  }) %>% set_names(as.character(years))
})

# ====== 4) Per-species processor ======
process_species <- function(sp) {
  mask_path <- sprintf(paths$mask_pattern, sp)
  check_exists(mask_path)
  mask_r <- rast(mask_path) %>% ensure_crs(dwd_crs)
  
  recs <- vector("list", length(years))
  
  for (i in seq_along(years)) {
    yr <- years[i]
    
    t_j_path <- paths$dwd_temp_j(yr)
    t_a_path <- paths$dwd_temp_a(yr)
    if (!file.exists(t_j_path) || !file.exists(t_a_path)) {
      warning(sprintf("[%s] Missing temperature grid(s) for %d; skipping year.", sp, yr))
      next
    }
    
    # Temperature (0.1 °C) for Jul & Aug
    t_j <- rast(t_j_path) %>% ensure_crs(dwd_crs)
    t_a <- rast(t_a_path) %>% ensure_crs(dwd_crs)
    tmp_src <- app(c(t_j, t_a), mean) / 10
    tmp_src <- ensure_crs(tmp_src, dwd_crs)
    
    # Project & mask to species area
    tmp <- tmp_src |> project(mask_r, method = "bilinear") |> mask(mask_r)
    
    # Saturation vapour pressure (hPa)
    es <- app(tmp, es.func)
    
    # Relative humidity mean for the two target days
    hy_src <- rh_layers_by_year[[as.character(yr)]]
    if (is.null(hy_src)) {
      warning(sprintf("[%s] No RH data for %d; skipping year.", sp, yr))
      next
    }
    hy_src <- ensure_crs(hy_src, hu_fallback)
    hy     <- hy_src |> project(mask_r, method = "bilinear") |> mask(mask_r)
    
    # VPD (kPa)
    vpd_r <- app(c(es, hy), vpd.func) / 10
    
    # Spatial means
    mean_tmp <- as.numeric(global(tmp,   "mean", na.rm = TRUE))
    mean_vpd <- as.numeric(global(vpd_r, "mean", na.rm = TRUE))
    
    recs[[i]] <- tibble(
      species = sp,
      time    = as.Date(sprintf("%d-08-01", yr)),
      temp    = mean_tmp,   # °C
      vpd     = mean_vpd    # kPa
    )
  }
  
  bind_rows(recs)
}

# ====== 5) Run & save time series ======
ts_df <- map_dfr(species_list, process_species) %>%
  mutate(species = factor(species, levels = species_list)) %>%
  arrange(species, time)

save(ts_df, file = paths$rdata_out)
write_csv(ts_df, paths$csv_out)
print(ts_df)

# ====== 6) Plotting (UNCHANGED theme & saver) ======
# --- Custom Theme (unchanged) ---
theme_rootzone <- function() {
  theme_minimal() +
    theme(
      axis.text.x        = element_text(hjust = 0.5),
      plot.background    = element_rect(fill = "white", color = NA),
      panel.background   = element_rect(fill = "white"),
      panel.border       = element_rect(color = "black", fill = NA),
      panel.grid         = element_blank(),
      legend.position    = "top",
      legend.text        = element_text(size = 14),
      legend.title       = element_text(size = 14, face = "bold"),
      axis.title         = element_text(face = "bold", size = 16),
      axis.text          = element_text(color = "black", size = 14),
      plot.title         = element_text(hjust = 0.5, size = 18, face = "bold"),
      plot.subtitle      = element_text(hjust = 0.5, size = 14),
      strip.background   = element_rect(fill = "white", color = "black"),
      strip.text         = element_text(face = "bold", size = 12)
    )
}

# --- Helper for nice year axis (unchanged) ---
year_breaks_fmt <- function(n = 6) {
  list(
    scale_x_date(
      breaks = pretty_breaks(n),
      date_labels = "%Y",
      expand = expansion(mult = c(0.01, 0.01))
    )
  )
}

# --- Two-row figure saver (unchanged) ---
save_temp_vpd_time_series <- function(df,
                                      out_dir = paths$fig_dir,
                                      cb_palette = c(
                                        "Oak"    = "#E69F00",
                                        "Beech"  = "#0072B2",
                                        "Spruce" = "#009E73",
                                        "Pine"   = "#F0E442"
                                      ),
                                      width = 10, height = 10, dpi = 300,
                                      filename = "TS_temp_VPD_all_species.png",
                                      caption_size = 14) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  df <- df %>%
    mutate(
      time    = as.Date(time),
      species = factor(species, levels = names(cb_palette))
    )
  
  p_temp_base <- ggplot(df, aes(x = time, y = temp, color = species)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.6) +
    scale_color_manual(values = cb_palette, name = "") +
    labs(x = NULL, y = "temperature (°C)") +
    theme_rootzone() +
    year_breaks_fmt()
  
  p_vpd_base <- ggplot(df, aes(x = time, y = vpd, color = species)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.6) +
    scale_color_manual(values = cb_palette, name = "") +
    labs(x = NULL, y = "VPD (kPa)") +
    theme_rootzone() +
    year_breaks_fmt()
  
  legend <- cowplot::get_legend(
    p_temp_base +
      theme(
        legend.position      = "top",
        legend.background    = element_rect(fill = "white", color = NA),
        legend.box.background= element_rect(fill = "white", color = NA),
        legend.key           = element_rect(fill = "white", color = NA)
      )
  )
  
  p_temp <- p_temp_base +
    theme(legend.position = "none",
          plot.margin = margin(t = 5, r = 6, b = 22, l = 6))
  
  p_vpd  <- p_vpd_base +
    theme(legend.position = "none",
          plot.margin = margin(t = 5, r = 6, b = 22, l = 6))
  
  p_temp_labeled <- cowplot::add_sub(
    p_temp, "(a)", x = 0.5, y = 0, hjust = 0.5, vjust = 1,
    fontface = "bold", size = caption_size
  )
  p_vpd_labeled <- cowplot::add_sub(
    p_vpd, "(b)", x = 0.5, y = 0, hjust = 0.5, vjust = 1,
    fontface = "bold", size = caption_size
  )
  
  p_combined <- cowplot::plot_grid(
    legend,
    p_temp_labeled,
    p_vpd_labeled,
    ncol = 1,
    rel_heights = c(0.12, 1, 1)
  )
  
  ggsave(file.path(out_dir, filename), p_combined,
         width = width, height = height, dpi = dpi, bg = "white")
  
  invisible(p_combined)
}

# ====== 7) Save figure ======
cb_palette <- c("Oak"="#E69F00","Beech"="#0072B2","Spruce"="#009E73","Pine"="#F0E442")
save_temp_vpd_time_series(ts_df, cb_palette = cb_palette)

### S10 year mean linear regression coefficients ###
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

plot_mean_linear_coeffs <- function(df_all, output_path) {
  # Ensure patchwork is available for combining plots
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    install.packages("patchwork")
  }
  library(patchwork)
  
  df_species <- df_average_year_species(df_all)
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  df_species$species <- factor(df_species$species, levels = species_order)
  
  ndvi_column <- if ("avg_quantile" %in% names(df_species)) "avg_quantile" else if ("avg_proportion" %in% names(df_species)) "avg_proportion" else stop("No NDVI column found.")
  ndvi_label <- if (ndvi_column == "avg_quantile") "NDVI quantiles" else "NDVI proportions"
  
  cb_palette <- c("Oak"   = "#E69F00", 
                  "Beech" = "#0072B2",
                  "Spruce"= "#009E73", 
                  "Pine"  = "#F0E442")
  
  # Regression Coefficients for NDVI ~ Soil Water Potential
  coeffs_psi <- do.call(rbind, lapply(species_order, function(sp) {
    df_sp <- df_species %>% filter(species == sp)
    mod <- lm(as.formula(paste0(ndvi_column, " ~ avg_psi")), data = df_sp)
    coef_vals <- coef(mod)
    p_vals <- summary(mod)$coefficients[,4]
    data.frame(species = sp, term = names(coef_vals), estimate = coef_vals, p.value = p_vals, stringsAsFactors = FALSE)
  }))
  coeffs_psi$term <- dplyr::recode(coeffs_psi$term, "(Intercept)" = "a", "avg_psi" = "b")
  coeffs_psi$species <- factor(coeffs_psi$species, levels = species_order)
  coeffs_psi$label <- ifelse(coeffs_psi$p.value < 0.05, "*", sprintf("%.2f", coeffs_psi$p.value))
  
  axis_title_size <- 16
  
  p_coeff_psi <- ggplot(coeffs_psi, aes(x = species, y = estimate, fill = species)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = label, y = estimate/2), position = position_dodge(width = 0.9), color = "black", size = 4) +
    scale_fill_manual(values = cb_palette, name ="") +
    facet_wrap(~ term, scales = "free_y") +
    labs(title = "NDVI quantiles (rank) ~ soil water potential (kPa)", y = "coefficient value", x = "", 
         tag = "(a)") +
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
      strip.text = element_text(face = "bold", size = 12),
      plot.tag.position = c(0.5, 0),                         # bottom center
      plot.tag = element_text(size = axis_title_size, 
                              face = "bold", hjust = 0.5)
    )
  
  # Regression Coefficients for NDVI ~ Transpiration Deficit
  coeffs_tdiff <- do.call(rbind, lapply(species_order, function(sp) {
    df_sp <- df_species %>% filter(species == sp)
    mod <- lm(as.formula(paste0(ndvi_column, " ~ avg_tdiff")), data = df_sp)
    coef_vals <- coef(mod)
    p_vals <- summary(mod)$coefficients[,4]
    data.frame(species = sp, term = names(coef_vals), estimate = coef_vals, p.value = p_vals, stringsAsFactors = FALSE)
  }))
  coeffs_tdiff$term <- dplyr::recode(coeffs_tdiff$term, "(Intercept)" = "a", "avg_tdiff" = "b")
  coeffs_tdiff$species <- factor(coeffs_tdiff$species, levels = species_order)
  coeffs_tdiff$label <- ifelse(coeffs_tdiff$p.value < 0.05, "*", sprintf("%.2f", coeffs_tdiff$p.value))
  
  p_coeff_tdiff <- ggplot(coeffs_tdiff, aes(x = species, y = estimate, fill = species)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = label, y = estimate/2), position = position_dodge(width = 0.9), color = "black", size = 4) +
    scale_fill_manual(values = cb_palette, name ="") +
    facet_wrap(~ term, scales = "free_y") +
    labs(title = "NDVI quantiles (rank) ~ transpiration deficit (mm)", y = "coefficient value", x = "", 
         tag = "(b)") +
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
      strip.text = element_text(face = "bold", size = 12),
      plot.tag.position = c(0.5, 0),                         # bottom center
      plot.tag = element_text(size = axis_title_size, 
                              face = "bold", hjust = 0.5)
    )
  
  # Combine with patchwork, collect legends
  combined_plot <- (p_coeff_psi + p_coeff_tdiff) +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "top")
  
  print(combined_plot)
  ggsave(filename = output_path, plot = combined_plot, width = 14, height = 8, dpi = 300)
}

plot_mean_linear_coeffs(data, "results_rootzone/Figures/supplementary/coeffs_yearly_mean.png")

### S10 AIC ###
NDVI_PSIbin <- function(df, bin_width = 50) {
  
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
  
  library(dplyr); library(tibble)
  
  df <- as_tibble(df)
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
  # Here the value is transpiration_deficit
  value_column <- "transpiration_deficit"
  
  species_totals <- df %>% group_by(species) %>% summarise(total_pixels = n(), .groups = "drop")
  
  psi_min <- floor(min(df$soil_water_potential, na.rm = TRUE))
  psi_max <- ceiling(max(df$soil_water_potential, na.rm = TRUE))
  bin_breaks <- seq(psi_min, psi_max, by = bin_width)
  
  df <- df %>%
    mutate(PSI_bin = cut(soil_water_potential, breaks = bin_breaks, include.lowest = TRUE, right = FALSE))
  
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
  
  # --- helpers (minimal, self-contained) ---
  get_lm_p <- function(fit) {
    fs <- summary(fit)$fstatistic
    if (is.null(fs)) return(NA_real_)
    pf(fs[1], fs[2], fs[3], lower.tail = FALSE)
  }
  get_nls_overall_p <- function(y, rss_model, k_model) {
    # Overall model test vs intercept-only using RSS reduction (approx. F-test)
    # y: response vector; rss_model: residual sum of squares of fitted model
    # k_model: number of parameters in the fitted model
    n <- length(y)
    rss_null <- sum((y - mean(y, na.rm = TRUE))^2)
    df1 <- k_model - 1
    df2 <- n - k_model
    if (any(c(df1, df2) <= 0) || is.na(rss_model) || is.na(rss_null) || rss_model >= rss_null) return(NA_real_)
    fstat <- ((rss_null - rss_model) / df1) / (rss_model / df2)
    pf(fstat, df1, df2, lower.tail = FALSE)
  }
  
  # Define common species order
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  
  #### Panel A: NDVI ~ PSIbin ####
  data_a <- NDVI_PSIbin(data)
  data_a <- na.omit(data_a)
  value_col_a <- "avg_value"
  data_a$species <- factor(data_a$species, levels = species_order)
  data_a <- data_a %>% mutate(x = -bin_median)
  
  start_list_a <- list(a = 5, b = 3, c = 0.001)
  control_params_a <- nls.control(maxiter = 1200, minFactor = 1e-9)
  
  aic_results_a <- list()
  for (sp in levels(data_a$species)) {
    sp_data <- data_a %>% filter(species == sp)
    
    # Linear model
    lm_linear <- lm(avg_value ~ x, data = sp_data)
    aic_linear <- AIC(lm_linear)
    r2_linear <- summary(lm_linear)$r.squared
    p_linear  <- get_lm_p(lm_linear)
    
    # Exponential model (nls)
    aic_exp <- NA
    r2_exp <- NA
    p_exp  <- NA
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
      p_exp  <- get_nls_overall_p(sp_data[[value_col_a]], ss_res, length(coef(nls_exp)))
    }
    
    sp_results <- data.frame(
      species = sp,
      Model   = c("linear", "exponential"),
      AIC     = c(aic_linear, aic_exp),
      R2      = c(r2_linear, r2_exp),
      P       = c(p_linear,  p_exp),
      stringsAsFactors = FALSE
    )
    sp_results$y_label_pos <- sp_results$AIC / 2
    # Label with star if p < 0.05
    sp_results$Label <- ifelse(!is.na(sp_results$P) & sp_results$P < 0.05,
                               paste0(round(sp_results$R2, 2), "*"),
                               as.character(round(sp_results$R2, 2)))
    aic_results_a[[sp]] <- sp_results
  }
  aic_df_a <- do.call(rbind, aic_results_a)
  aic_df_a$species <- factor(aic_df_a$species, levels = species_order)
  
  # Shared color palette for Panels A and B
  model_palette_shared <- c("linear" = "orange",
                            "exponential" = "#0072B2")
  
  p_a <- ggplot(aic_df_a, aes(x = species, y = AIC, fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(aes(label = Label, y = y_label_pos),
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
    p_linear  <- get_lm_p(lm_linear)
    
    aic_exp <- NA
    r2_exp <- NA
    p_exp  <- NA
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
      p_exp  <- get_nls_overall_p(sp_data[[value_col_b]], ss_res, length(coef(nls_exp)))
    }
    
    sp_results <- data.frame(
      species = sp,
      Model   = c("linear", "exponential"),
      AIC     = c(aic_linear, aic_exp),
      R2      = c(r2_linear, r2_exp),
      P       = c(p_linear,  p_exp),
      stringsAsFactors = FALSE
    )
    sp_results$y_label_pos <- sp_results$AIC / 2
    sp_results$Label <- ifelse(!is.na(sp_results$P) & sp_results$P < 0.05,
                               paste0(round(sp_results$R2, 2), "*"),
                               as.character(round(sp_results$R2, 2)))
    aic_results_b[[sp]] <- sp_results
  }
  aic_df_b <- do.call(rbind, aic_results_b)
  aic_df_b$species <- factor(aic_df_b$species, levels = species_order)
  
  p_b <- ggplot(aic_df_b, aes(x = species, y = AIC, fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(aes(label = Label, y = y_label_pos),
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
    p_linear  <- get_lm_p(lm_linear)
    
    lm_poly2 <- lm(avg_transpiration_deficit ~ x + I(x^2), data = sp_data)
    aic_poly2 <- AIC(lm_poly2)
    r2_poly2 <- summary(lm_poly2)$r.squared
    p_poly2  <- get_lm_p(lm_poly2)
    
    lm_poly3 <- lm(avg_transpiration_deficit ~ x + I(x^2) + I(x^3), data = sp_data)
    aic_poly3 <- AIC(lm_poly3)
    r2_poly3 <- summary(lm_poly3)$r.squared
    p_poly3  <- get_lm_p(lm_poly3)
    
    sp_results <- data.frame(
      species = sp,
      Model   = c("linear", "poly2", "poly3"),
      AIC     = c(aic_linear, aic_poly2, aic_poly3),
      R2      = c(r2_linear, r2_poly2, r2_poly3),
      P       = c(p_linear,  p_poly2,  p_poly3),
      stringsAsFactors = FALSE
    )
    sp_results$y_label_pos <- sp_results$AIC / 2
    sp_results$Label <- ifelse(!is.na(sp_results$P) & sp_results$P < 0.05,
                               paste0(round(sp_results$R2, 2), "*"),
                               as.character(round(sp_results$R2, 2)))
    aic_results_c[[sp]] <- sp_results
  }
  aic_df_c <- do.call(rbind, aic_results_c)
  aic_df_c$species <- factor(aic_df_c$species, levels = species_order)
  
  model_palette_c <- c("linear" = "orange",
                       "poly2"  = "#0072B2",
                       "poly3"  = "#009E73")
  
  p_c <- ggplot(aic_df_c, aes(x = species, y = AIC, fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(aes(label = Label, y = y_label_pos),
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
  
  # ---- print summary tables as requested ----
  cat("\n=== Panel A: NDVI ~ PSIbin ===\n")
  print(aic_df_a %>% select(species, Model, AIC, R2, P))
  cat("\n=== Panel B: NDVI ~ TDiffbin ===\n")
  print(aic_df_b %>% select(species, Model, AIC, R2, P))
  cat("\n=== Panel C: TDiff ~ PSIbin ===\n")
  print(aic_df_c %>% select(species, Model, AIC, R2, P))
  cat("\nNote: Asterisk (*) after R² in bars indicates p < 0.05 for the overall model test.\n")
  
  print(combined_plot)
}

plot_combined_AIC_R2(data, "results_rootzone/Figures/supplementary/AIC.png")

### S11 TDiff PSIbin ###
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

plot_TDiff_PSIbin_poly_2_slope(data, 
                               "results_rootzone/Figures/supplementary/TDiff_PSI_poly2_coeff.png",
                               "results_rootzone/Figures/supplementary/TDiff_PSI_poly2_slope.png")

plot_TDiff_PSIbin_poly_3_slope(data, 
                               "results_rootzone/Figures/supplementary/TDiff_PSI_poly3_coeff.png",
                               "results_rootzone/Figures/supplementary/TDiff_PSI_poly3_slope.png")


plot_TDiff_PSIbin_slope <- function(data,
                                    coef_output,
                                    figure_output,
                                    aic_barplot_fig) {
  # -------------------------
  # Libraries & preprocessing
  # -------------------------
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(tibble)
  library(purrr)
  library(scales)  # percent_format
  
  # Prepare data
  df <- TDiff_PSIbin(data)
  df <- na.omit(df)
  
  # Species & palette (match your first function)
  cb_palette <- c("Oak"   = "#E69F00",
                  "Beech" = "#0072B2",
                  "Spruce"= "#009E73",
                  "Pine"  = "#F0E442")
  species_levels <- c("Oak", "Beech", "Spruce", "Pine")
  df$species <- factor(df$species, levels = species_levels)
  
  # shorthand
  xcol <- "bin_median"
  ycol <- "avg_transpiration_deficit"
  
  # -------------------------
  # Fit models & choose by AIC
  # -------------------------
  degrees <- c(1, 2, 3)
  fit_one_species <- function(sp) {
    sdf <- df %>% filter(species == sp) %>% select(all_of(c(xcol, ycol)))
    if (nrow(sdf) < 5 || any(!is.finite(sdf[[xcol]])) || any(!is.finite(sdf[[ycol]]))) {
      return(list(species = sp,
                  aic_tbl = tibble(model = c("linear","poly2","poly3"), AIC = NA_real_),
                  best_deg = NA_integer_, best_fit = NULL))
    }
    fits <- lapply(degrees, function(d) {
      fm <- reformulate(termlabels = paste0("poly(", xcol, ", ", d, ")"), response = ycol)
      tryCatch(lm(fm, data = sdf), error = function(e) NULL)
    })
    aics <- sapply(fits, function(m) if (is.null(m)) NA_real_ else AIC(m))
    aic_tbl <- tibble(model = c("linear","poly2","poly3"), AIC = aics)
    best_idx <- suppressWarnings(which.min(aics))
    best_deg <- if (length(best_idx) == 0 || all(is.na(aics))) NA_integer_ else degrees[best_idx]
    best_fit <- if (is.na(best_deg)) NULL else fits[[best_idx]]
    list(species = sp, aic_tbl = aic_tbl, best_deg = best_deg, best_fit = best_fit)
  }
  fits_all <- lapply(species_levels, fit_one_species)
  
  # AIC dataframe (long) for barplot
  aic_df <- bind_rows(lapply(fits_all, function(x) {
    x$aic_tbl %>% mutate(species = x$species)
  }))
  aic_df$species <- factor(aic_df$species, levels = species_levels)
  aic_df$model <- factor(aic_df$model, levels = c("linear","poly2","poly3"))
  
  message("AIC values per species:")
  print(aic_df)
  
  # -------------------------
  # AIC barplot (saved)
  # -------------------------
  p_aic <- ggplot(aic_df, aes(x = species, y = AIC, fill = model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
    scale_fill_manual(values = c("linear" = "#0072B2",
                                 "poly2"  = "#E69F00",
                                 "poly3"  = "#009E73"),
                      name = "") +
    labs(x = "", y = "AIC") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
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
  print(p_aic)
  ggsave(filename = aic_barplot_fig, plot = p_aic, device = "png", width = 8, height = 6, dpi = 300)
  
  # ---------------------------------
  # Predictions from the best models
  # ---------------------------------
  pred_list <- lapply(fits_all, function(ff) {
    sp <- ff$species
    if (is.null(ff$best_fit) || is.na(ff$best_deg)) return(NULL)
    sp_df <- df %>% filter(species == sp)
    xr <- range(sp_df[[xcol]], na.rm = TRUE)
    if (!is.finite(xr[1]) || !is.finite(xr[2]) || xr[1] == xr[2]) return(NULL)
    xseq <- seq(xr[1], xr[2], length.out = 200)
    nd <- data.frame(bin_median = xseq)
    nd$species <- factor(sp, levels = species_levels)
    nd$predicted <- predict(ff$best_fit, newdata = nd)
    nd$degree <- ff$best_deg
    nd
  })
  pred_df <- bind_rows(pred_list)
  
  # Lines for mean / median based on observed y
  line_mean   <- mean(df[[ycol]], na.rm = TRUE)
  line_median <- median(df[[ycol]], na.rm = TRUE)
  
  print(line_median)
  
  # -------------
  # Panel (a): scatter + fitted lines + median reference
  # -------------
  plot_mixed <- ggplot() +
    geom_point(data = df,
               aes(x = .data[[xcol]],
                   y = .data[[ycol]],
                   color = species,
                   shape = species,
                   size = percentage),
               alpha = 0.7) +
    geom_line(data = pred_df,
              aes(x = bin_median, y = predicted, color = species),
              linewidth = 1) +
    geom_hline(yintercept = line_median, linetype = "dashed", color = "black", linewidth = 1) +
    annotate("text",
             x = -2000,
             y = line_median,
             label = "median",
             vjust = -0.3, fontface = "italic", size = 6) +
    scale_color_manual(values = cb_palette, name = "") +
    scale_shape_manual(values = c("Oak" = 16, "Beech" = 17, "Spruce" = 15, "Pine" = 18),
                       guide = "none") +
    scale_size_continuous(name = "Pixels per bin (%)",
                          range = c(1, 8),
                          labels = percent_format(accuracy = 1)) +
    guides(
      color = guide_legend(order = 1),
      size  = guide_legend(order = 2)
    ) +
    labs(x = "soil water potential (kPa)",
         y = "transpiration deficit") +
    ggtitle("(a)") +
    theme_minimal() +
    theme(
      plot.title.position = "plot",
      plot.title          = element_text(face = "bold", size = 18, hjust = 0.1, vjust = 1),
      axis.text.x         = element_text(angle = 0, hjust = 0.5),
      axis.title          = element_text(face = "bold", size = 16),
      axis.text           = element_text(color = "black", size = 14),
      legend.position     = "bottom",
      legend.text         = element_text(size = 14),
      legend.title        = element_text(size = 14, face = "bold"),
      panel.background    = element_rect(fill = "white"),
      plot.background     = element_rect(fill = "white", color = "white"),
      panel.grid          = element_blank(),
      panel.border        = element_blank()
    )
  
  # -------------
  # Panel (b): soil water potential at median transpiration deficit
  #   -> **root finding only** with automatic bracketing (no interpolation fallback)
  # -------------
  solve_x_at_y <- function(fit, y_target, x_range, n_grid = 400) {
    # function for root finding
    f <- function(x) predict(fit, newdata = data.frame(bin_median = x)) - y_target
    
    # quick check at the endpoints
    fx_lo <- f(x_range[1]); fx_hi <- f(x_range[2])
    
    # build a grid to locate a sign change bracket
    xs <- seq(x_range[1], x_range[2], length.out = n_grid)
    fs <- suppressWarnings(f(xs))
    
    # find first sign change interval
    sc_idx <- which(diff(sign(fs)) != 0 & is.finite(fs[-length(fs)]) & is.finite(fs[-1]))
    if (length(sc_idx) == 0) {
      # no sign change -> no root in range
      return(NA_real_)
    }
    # choose the bracket where |f| near zero is smallest
    pick <- sc_idx[which.min(pmin(abs(fs[sc_idx]), abs(fs[sc_idx+1])))]
    lo <- xs[pick]; hi <- xs[pick + 1]
    
    # refine with uniroot
    tryCatch(uniroot(function(z) f(z), interval = c(lo, hi))$root,
             error = function(e) NA_real_)
  }
  
  # compute per species (root finding only)
  x_at_median <- lapply(fits_all, function(ff) {
    sp <- ff$species
    fit <- ff$best_fit
    if (is.null(fit)) return(tibble(species = sp, x_at_median = NA_real_))
    
    sp_df <- df %>% filter(species == sp)
    xr <- range(sp_df[[xcol]], na.rm = TRUE)
    xm <- solve_x_at_y(fit, line_median, xr, n_grid = 600)
    tibble(species = sp, x_at_median = xm)
  }) %>% bind_rows() %>% mutate(species = factor(species, levels = species_levels))
  
  p_bar_x <- ggplot(x_at_median, aes(x = species, y = x_at_median, fill = species)) +
    geom_bar(stat = "identity", width = 0.7, na.rm = TRUE) +
    scale_fill_manual(values = cb_palette, guide = "none") +
    labs(x = "", y = "soil water potential (kPa)") +
    ggtitle("(b)") +
    theme_minimal() +
    theme(
      plot.title.position = "plot",
      plot.title          = element_text(face = "bold", size = 18, hjust = 0.25, vjust = 1),
      axis.text.x         = element_text(angle = 0, hjust = 0.5),
      axis.title          = element_text(face = "bold", size = 16),
      axis.text           = element_text(color = "black", size = 14),
      legend.position     = "none",
      panel.background    = element_rect(fill = "white"),
      plot.background     = element_rect(fill = "white", color = "white"),
      panel.grid          = element_blank(),
      panel.border        = element_blank()
    )
  
  # -------------
  # Panel (c): local slope near 30% of max(pred)
  # -------------
  local_slope_data <- pred_df %>%
    group_by(species) %>%
    group_modify(~{
      d <- .x
      maxp <- max(d$predicted, na.rm = TRUE)
      thresh <- 0.3 * maxp
      idx <- which.min(abs(d$predicted - thresh))
      w0 <- max(1, idx - 5)
      w1 <- min(nrow(d), idx + 5)
      dw <- d[w0:w1, , drop = FALSE]
      lm_loc <- lm(predicted ~ bin_median, data = dw)
      sm <- summary(lm_loc)
      data.frame(
        slope = coef(lm_loc)[["bin_median"]],
        slope_se = sm$coefficients["bin_median","Std. Error"],
        p_value = sm$coefficients["bin_median","Pr(>|t|)"],
        r2 = sm$r.squared
      )
    }) %>% ungroup() %>%
    mutate(slope_abs = abs(slope),
           label_text = ifelse(p_value < 0.05,
                               sprintf("%.2f*", r2),
                               sprintf("%.2f", r2)))
  
  p_bar <- ggplot(local_slope_data, aes(x = species, y = slope_abs, fill = species)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, slope_abs - slope_se), ymax = slope_abs + slope_se), width = 0.2) +
    geom_text(aes(y = slope_abs/2, label = label_text), size = 5, color = "black") +
    scale_fill_manual(values = cb_palette, guide = "none") +
    labs(x = "", y = "absolute slope") +
    ggtitle("(c)") +
    theme_minimal() +
    theme(
      plot.title.position = "plot",
      plot.title          = element_text(face = "bold", size = 18, hjust = 0.25, vjust = 1),
      axis.text.x         = element_text(angle = 0, hjust = 0.5),
      axis.title          = element_text(face = "bold", size = 16),
      axis.text           = element_text(color = "black", size = 14),
      legend.position     = "none",
      panel.background    = element_rect(fill = "white"),
      plot.background     = element_rect(fill = "white", color = "white"),
      panel.grid          = element_blank(),
      panel.border        = element_blank()
    )
  
  # -------------
  # Combine panels
  # -------------
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
  
  # -------------------------
  # Coefficients (best model per species)
  # -------------------------
  coef_rows <- lapply(fits_all, function(ff) {
    sp <- ff$species
    fit <- ff$best_fit
    deg <- ff$best_deg
    if (is.null(fit) || is.na(deg)) return(NULL)
    sm <- summary(fit)$coefficients
    out <- as.data.frame(sm)
    out$term_raw <- rownames(out)
    out$species <- sp
    out$degree <- deg
    out
  }) %>% bind_rows()
  
  if (!is.null(coef_rows) && nrow(coef_rows) > 0) {
    coef_data <- coef_rows %>%
      mutate(
        term = dplyr::case_when(
          term_raw == "(Intercept)" ~ "a",
          grepl("poly\\(bin_median,\\s*1\\)1", term_raw) ~ "b",
          grepl("poly\\(bin_median,\\s*2\\)1", term_raw) ~ "b",
          grepl("poly\\(bin_median,\\s*2\\)2", term_raw) ~ "c",
          grepl("poly\\(bin_median,\\s*3\\)1", term_raw) ~ "b",
          grepl("poly\\(bin_median,\\s*3\\)2", term_raw) ~ "c",
          grepl("poly\\(bin_median,\\s*3\\)3", term_raw) ~ "d",
          TRUE ~ term_raw
        ),
        species = factor(species, levels = species_levels)
      ) %>%
      dplyr::rename(value = Estimate, p_value = `Pr(>|t|)`) %>%
      mutate(label_text = ifelse(p_value < 0.05, "*", sprintf("%.2f", p_value)),
             term = factor(term, levels = c("a","b","c","d")))
    
    # Keep terms up to selected degree for each species
    coef_data <- coef_data %>%
      group_by(species) %>%
      filter(
        (term %in% c("a","b") & degree >= 1) |
          (term %in% c("c")    & degree >= 2) |
          (term %in% c("d")    & degree >= 3)
      ) %>% ungroup()
    
    plot_coeff <- ggplot(coef_data, aes(x = species, y = value, fill = species)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
      geom_text(aes(label = label_text, y = value/2),
                color = "black", size = 5,
                position = position_dodge(width = 0.8)) +
      facet_wrap(~ term, scales = "free_y") +
      scale_fill_manual(values = cb_palette, name = "") +
      labs(title = "transpiration deficit (mm) ~ soil water potential (kPa)",
           x = "", y = "coefficient value") +
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
    
    print(plot_coeff)
    ggsave(coef_output, plot = plot_coeff, width = 10, height = 8, dpi = 300)
    
    best_summary <- tibble(
      species = sapply(fits_all, `[[`, "species"),
      best_degree = sapply(fits_all, `[[`, "best_deg")
    )
    message("Best polynomial degree by AIC:")
    print(best_summary)
  } else {
    warning("No coefficients available to plot (fits failed or insufficient data).")
    plot_coeff <- NULL
  }
  
  invisible(list(combined_plot = final_plot,
                 coeff_plot = plot_coeff,
                 aic_plot = p_aic))
}


plot_TDiff_PSIbin_slope(
  data = data,
  coef_output = "results_rootzone/Figures/supplementary/TDiff_PSIbin_coeff.png",
  figure_output = "results_rootzone/Figures/supplementary/TDiff_PSIbin_slope.png",
  aic_barplot_fig = "results_rootzone/Figures/supplementary/TDiff_PSIbin_AIC.png")

### S7 soil texture composition ###
source("R/plot_soil_map_MODIS.R")

plot_soil_fraction_composite <- function(
    species_order = c("Oak", "Beech", "Spruce", "Pine"),
    soil_map_dir = "soil_map",
    result_data_dir = "results_rootzone/Data",
    result_fig_dir = "results_rootzone/supplementary/Figures",
    lut_path = file.path(soil_map_dir, "feinbod_lookup_with_english.csv"),
    texture_path = file.path(soil_map_dir, "soil_texture_classes.csv"),
    output_plot = file.path(result_fig_dir, "all_species_soil_fractions_composite.png"),
    output_csv  = file.path(result_data_dir, "soil_fraction_values.csv"),
    plot_width  = 12,
    plot_height = 12,
    plot_dpi    = 300
) {
  # Load required libraries
  library(terra)
  library(dplyr)
  library(ggplot2)
  library(sf)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(patchwork)
  
  # Read Germany boundary
  boundary_germany <- ne_countries(
    scale = "medium",
    country = "Germany",
    returnclass = "sf"
  )
  
  # Read lookup tables and build LUT
  soil_descrip <- read.csv(lut_path, stringsAsFactors = FALSE)
  soil_texture <- read.csv(texture_path, stringsAsFactors = FALSE)
  lut <- soil_descrip %>%
    select(feinbod_code, feinbod) %>%
    inner_join(soil_texture, by = c("feinbod" = "code")) %>%
    mutate(
      clay_mid = (clay_lower + clay_upper) / 2,
      silt_mid = (silt_lower + silt_upper) / 2,
      sand_mid = (sand_lower + sand_upper) / 2
    ) %>%
    select(feinbod_code, clay_mid, silt_mid, sand_mid)
  
  # Function to generate fraction plots and data for one raster
  make_fraction_plots <- function(raster_path) {
    soil_code <- rast(raster_path)
    clay_r <- subst(soil_code, from = lut$feinbod_code, to = lut$clay_mid); names(clay_r) <- "clay"
    silt_r <- subst(soil_code, from = lut$feinbod_code, to = lut$silt_mid); names(silt_r) <- "silt"
    sand_r <- subst(soil_code, from = lut$feinbod_code, to = lut$sand_mid); names(sand_r) <- "sand"
    
    clay_df <- as.data.frame(clay_r, xy = TRUE) %>% rename(value = clay) %>% filter(!is.na(value))
    silt_df <- as.data.frame(silt_r, xy = TRUE) %>% rename(value = silt) %>% filter(!is.na(value))
    sand_df <- as.data.frame(sand_r, xy = TRUE) %>% rename(value = sand) %>% filter(!is.na(value))
    
    base_theme <- theme_minimal() +
      theme(
        axis.text.x      = element_text(angle = 0, hjust = 0.5),
        plot.background  = element_rect(fill = "white", color = "white"),
        panel.background = element_rect(fill = "white"),
        legend.background= element_rect(fill = "white", color = "white"),
        plot.title       = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
        axis.title       = element_blank(),
        axis.text        = element_text(color = "black", size = 14),
        panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position  = "bottom",
        legend.key.width = unit(1.3, "cm"),
        legend.key.height= unit(0.5, "cm"),
        legend.text      = element_text(size = 14),
        legend.title     = element_text(size = 14),
        strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
        strip.text       = element_text(face = "bold", size = 12)
      )
    
    list(
      clay = ggplot(clay_df, aes(x = x, y = y, color = value)) +
        geom_point(size = 0.5) +
        geom_sf(data = boundary_germany, fill = NA, color = "black", inherit.aes = FALSE) +
        coord_sf(expand = FALSE) +
        scale_color_gradientn(colours = c("white", "lightblue", "dodgerblue", "#0072B2"), name = "clay (%)") +
        base_theme,
      silt = ggplot(silt_df, aes(x = x, y = y, color = value)) +
        geom_point(size = 0.5) +
        geom_sf(data = boundary_germany, fill = NA, color = "black", inherit.aes = FALSE) +
        coord_sf(expand = FALSE) +
        scale_color_gradientn(colours = c("white", "#b9f6ca", "#00bfae", "#00675b"), name = "silt (%)") +
        base_theme,
      sand = ggplot(sand_df, aes(x = x, y = y, color = value)) +
        geom_point(size = 0.5) +
        geom_sf(data = boundary_germany, fill = NA, color = "black", inherit.aes = FALSE) +
        coord_sf(expand = FALSE) +
        scale_color_gradientn(colours = c("white", "#ffe0b2", "orange", "#ff6600"), name = "sand (%)") +
        base_theme,
      data = bind_rows(
        clay_df %>% mutate(fraction = "clay"),
        silt_df %>% mutate(fraction = "silt"),
        sand_df %>% mutate(fraction = "sand")
      )
    )
  }
  
  # Prepare storage for plots and data
  plot_list <- list()
  soil_fraction_df <- tibble()
  
  # Generate plots per species
  for (sp in species_order) {
    raster_file <- file.path(soil_map_dir, paste0(sp, "_soilCode_MODIS.tif"))
    res <- make_fraction_plots(raster_file)
    plot_list[[sp]] <- res[c("clay", "silt", "sand")]
    soil_fraction_df <- bind_rows(soil_fraction_df, res$data %>% mutate(species = sp))
  }
  
  # Define desired row (fraction) order
  fraction_order <- c("silt", "clay", "sand")
  
  # Build composite list: row-by-row, adding species title only on first row
  all_plots_list <- list()
  for (fr in fraction_order) {
    for (i in seq_along(species_order)) {
      sp <- species_order[i]
      p <- plot_list[[sp]][[fr]] +
        theme(plot.margin = margin(5, 5, 2, 5))
      # Add title only on first row
      if (fr == fraction_order[1]) {
        p <- p + labs(title = sp)
      }
      all_plots_list <- c(all_plots_list, list(p))
    }
  }
  
  # Assemble and save composite
  composite <- wrap_plots(all_plots_list,
                          ncol = length(species_order),
                          nrow = length(fraction_order),
                          guides = "collect") &
    theme(
      legend.position = "bottom",
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin     = margin(5, 5, 5, 5),
      panel.spacing   = unit(0, "cm")
    )
  
  ggsave(output_plot, composite,
         width = plot_width,
         height = plot_height,
         dpi = plot_dpi,
         bg = "white")
  
  # Save raw data
  write.csv(soil_fraction_df, output_csv, row.names = FALSE)
  
  message("Composite plot saved to: ", output_plot)
  message("Soil fraction data saved to: ", output_csv)
}

plot_soil_fraction_composite()

### S8 soil texture composition bar ###
plot_bar_soil_composition <- function(
    data_csv     = "results_rootzone/Data/soil_fraction_values.csv",
    result_fig_dir = "results_rootzone/Figures/supplementary",
    output_file  = file.path(result_fig_dir, "soil_composition_bar.png"),
    widths       = 12,
    height       = 4,
    dpi          = 300
) {
  library(readr); library(dplyr); library(ggplot2)
  
  # Read raw fraction data
  df <- read_csv(data_csv)
  
  # Compute mean & proportions
  mean_df <- df %>%
    group_by(species, fraction) %>%
    summarise(mean_value = mean(value, na.rm=TRUE), .groups="drop") %>%
    group_by(species) %>%
    mutate(proportion = 100 * mean_value / sum(mean_value)) %>%
    ungroup() %>%
    mutate(
      species  = factor(species, levels=c("Oak","Beech","Spruce","Pine")),
      fraction = factor(fraction, levels=c("clay","silt","sand"))
    )
  
  # Bar theme (reuse bar_theme defined earlier if in global env)
  bar_theme <- theme_minimal() +
    theme(axis.text.x=element_text(size=14, color="black"),
          axis.text.y=element_text(size=14, color="black"),
          axis.title.x=element_blank(),
          axis.title.y=element_text(size=14, color="black"),
          legend.position="none",
          plot.background=element_rect(fill="white",color="white"),
          panel.border=element_rect(color="black",fill=NA),
          panel.grid=element_blank(),
          strip.text=element_text(size=16,face="bold"))
  
  # Plot three-panel bar chart
  p <- ggplot(mean_df, aes(x=species, y=proportion, fill=species)) +
    geom_col(width=0.7) +
    facet_wrap(~ fraction, nrow=1, scales="free_y") +
    scale_fill_manual(values=c(
      "Oak"    = "#E69F00",
      "Beech"  = "#0072B2",
      "Spruce" = "#009E73",
      "Pine"   = "#F0E442"
    )) +
    ylim(0,80) +
    labs(y="percentage (%)") +
    bar_theme
  
  print(p)
  ggsave(output_file, p, width=widths, height=height, dpi=dpi)
}
plot_bar_soil_composition()

### S9 LAI MODIS ###
source("R/plot_LAI.R")
plot_lai_map <- function(df, var = "LAI", title,
                         boundary = boundary_germany,
                         palette = c("white", "#e8fce8", "#bafcb6", "lightgreen", "#59d269", "#228b22", "darkgreen"),
                         limits = c(0, 5)) {
  ggplot(df) +
    geom_point(aes(x = x, y = y, color = .data[[var]]), size = 0.5) +
    geom_sf(data = boundary, fill = NA, color = "black", inherit.aes = FALSE, linewidth = 0.8) +
    scale_color_gradientn(
      colours = palette,
      values = scales::rescale(seq(limits[1], limits[2], length.out = length(palette)), from = limits),
      limits = limits,
      oob = scales::squish,
      name = "leaf area index (LAI)"
    ) +
    facet_wrap(~species, nrow = 1) +
    coord_sf(crs = 4326, expand = FALSE) +
    labs(x = "longitude", y = "latitude") +
    guides(color = guide_colorbar(title.position = "left", title.hjust = 0.5)) +
    custom_theme +
    theme(axis.text.x = element_text(angle = 0, hjust = 1), 
          legend.title = element_text(size = 16, face = "bold"))
}

plot_mean_lai_bar <- function(df,
                              var = "LAI",
                              species_levels = c("Oak", "Beech", "Spruce", "Pine"),
                              title = "",
                              palette = cb_palette) {
  summary_df <- df %>%
    group_by(species) %>%
    summarize(mean_value = mean(.data[[var]], na.rm = TRUE)) %>%
    mutate(species = factor(species, levels = species_levels))
  
  ggplot(summary_df, aes(x = species, y = mean_value, fill = species)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = palette, name = "") +
    labs(x = "", y = "mean leaf area index (LAI)", title = title) +
    custom_theme +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      legend.title = element_text(size = 16, face = "bold")
    )
}

p_map <- plot_lai_map(lai_df)
ggsave("results_rootzone/Figures/mean_LAI_map.png", p_map, width = 14, height = 6)

p_bar <- plot_mean_lai_bar(lai_df)
ggsave("results_rootzone/Figures/mean_LAI_bar.png", p_bar, width = 8, height = 6)

### Different match ###
source("R/plot_rootzone_diff_NDVI_PSI.R")
source("R/plot_diff_NDVI_TDiff.R")

### Different depth analysis ###
