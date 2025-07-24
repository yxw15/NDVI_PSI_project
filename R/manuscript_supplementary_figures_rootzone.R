setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")
load("results/Data/AllSpecies_AllMonths_rootzone.RData")

data_all <- combined %>% filter(Quantiles > 0) %>% filter(month %in% c("May", "June", "July", "August"))
data_all <- na.omit(data_all)
data <- data_all %>% filter(month %in% c("July", "August"))

library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)

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
      title = "NDVI quantile (rank) ~ soil water potential (kPa)",
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
                                        "results_rootzone/Figures/NDVI_PSI_linear_slope_monthly_orig.png",
                                        "results_rootzone/Figures/NDVI_PSI_linear_coeff_monthly_orig.png")

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

plot_NDVI_TDiff_month_species_linear_orig(data_all, 
                                          "results_rootzone/Figures/NDVI_TDiff_linear_slope_monthly_orig.png",
                                          "results_rootzone/Figures/NDVI_TDiff_linear_coeff_monthly_orig.png")

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
                             output_path = "results_rootzone/Figures/NDVI_PSI_density.png")

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

plot_density_NDVI_TDiff_linear(data,
                               tdiff_bin_width = 3,
                               output_path = "results_rootzone/Figures/NDVI_Tdiff_density.png")

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
  output_file  = "results_rootzone/Figures/mean_soil_water_potential.png"
)

plot_distribution_PSI_TDiff(
  data         = data,
  var          = "transpiration_deficit",
  legend_title = "transpiration deficit (mm)",
  output_file  = "results_rootzone/Figures/mean_transpiration_deficit.png"
)

### S5 mean bar PSI TDiff ###
plot_bar_mean_rootzone <- function(df, output_file) {
  
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

plot_bar_mean_rootzone(data, "results_rootzone/Figures/bar_mean_time_series_rootzone.png")

### S6 year mean linear regression coefficients ###
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

plot_yearly_mean_linear_coeffs <- function(df_all, output_path) {
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
  
  p_coeff_psi <- ggplot(coeffs_psi, aes(x = species, y = estimate, fill = species)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = label, y = estimate/2), position = position_dodge(width = 0.9), color = "black", size = 4) +
    scale_fill_manual(values = cb_palette, name ="") +
    facet_wrap(~ term, scales = "free_y") +
    labs(title = "NDVI (rank) ~ soil water potential (kPa)", y = "coefficient value", caption = "(a)") +
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
    labs(title = "NDVI (rank) ~ transpiration deficit (mm)", caption = "(b)") +
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
  
  # Combine with patchwork, collect legends
  combined_plot <- (p_coeff_psi + p_coeff_tdiff) +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "top")
  
  print(combined_plot)
  ggsave(filename = output_path, plot = combined_plot, width = 14, height = 8, dpi = 300)
}

plot_yearly_mean_linear_coeffs(data, "results_rootzone/Figures/coeffs_yearly_mean.png")

### S7 AIC ###
NDVI_PSIbin <- function(df, bin_width = 50) {
  value_column <- if ("Quantiles" %in% names(df)) "Quantiles" else "Proportions"
  
  species_totals <- df %>% group_by(species) %>% summarise(total_pixels = n(), .groups = "drop")
  
  psi_min <- floor(min(df$soil_water_potential, na.rm = TRUE))
  psi_max <- ceiling(max(df$soil_water_potential, na.rm = TRUE))
  bin_breaks <- seq(psi_min, psi_max, by = bin_width)
  
  df <- df %>%
    mutate(PSI_bin = cut(soil_water_potential, breaks = bin_breaks, include.lowest = TRUE, right = FALSE))
  
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
  value_col <- if ("Quantiles" %in% names(df)) "Quantiles" else "Proportions"
  
  species_totals <- df %>% group_by(species) %>% summarise(total_pixels = n(), .groups = "drop")
  
  tdiff_min <- floor(min(df$transpiration_deficit, na.rm = TRUE))
  tdiff_max <- ceiling(max(df$transpiration_deficit, na.rm = TRUE))
  bin_breaks <- seq(tdiff_min, tdiff_max, by = bin_width)
  
  df <- df %>%
    mutate(TDiff_bin = cut(transpiration_deficit, breaks = bin_breaks, include.lowest = TRUE, right = FALSE))
  
  get_bin_median <- function(bin_label) {
    nums <- as.numeric(strsplit(gsub("\\[|\\]|\\(|\\)", "", bin_label), ",")[[1]])
    mean(nums)
  }
  
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
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  
  # Panel A: NDVI ~ PSIbin
  data_a <- NDVI_PSIbin(data)
  data_a <- na.omit(data_a)
  data_a$species <- factor(data_a$species, levels = species_order)
  data_a <- data_a %>% mutate(x = -bin_median)
  
  start_list_a <- list(a = 5, b = 3, c = 0.001)
  control_params_a <- nls.control(maxiter = 200, minFactor = 1e-4)
  
  aic_results_a <- list()
  for (sp in levels(data_a$species)) {
    sp_data <- data_a %>% filter(species == sp)
    lm_linear <- lm(avg_value ~ x, data = sp_data)
    aic_linear <- AIC(lm_linear)
    r2_linear <- summary(lm_linear)$r.squared
    
    aic_exp <- NA; r2_exp <- NA
    nls_exp <- tryCatch({
      nls(avg_value ~ a + b * exp(-c * x),
          data = sp_data, start = start_list_a, control = control_params_a)
    }, error = function(e) NULL)
    if (!is.null(nls_exp)) {
      aic_exp <- AIC(nls_exp)
      res <- resid(nls_exp)
      ss_res <- sum(res^2)
      ss_tot <- sum((sp_data$avg_value - mean(sp_data$avg_value))^2)
      r2_exp <- 1 - ss_res / ss_tot
    }
    
    sp_results <- data.frame(species = sp, Model = c("Linear", "Exponential"),
                             AIC = c(aic_linear, aic_exp), R2 = c(r2_linear, r2_exp))
    sp_results$y_label_pos <- sp_results$AIC / 2
    aic_results_a[[sp]] <- sp_results
  }
  aic_df_a <- do.call(rbind, aic_results_a)
  aic_df_a$species <- factor(aic_df_a$species, levels = species_order)
  
  model_palette_shared <- c("Linear" = "orange", "Exponential" = "dodgerblue")
  
  p_a <- ggplot(aic_df_a, aes(x = species, y = AIC, fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(aes(label = round(R2, 2), y = y_label_pos),
              position = position_dodge(width = 0.9), vjust = 0.5, size = 8, color = "black") +
    labs(x = "", y = "AIC", title = "NDVI ~ PSIbin", caption = "(a)") +
    scale_fill_manual(values = model_palette_shared) +
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
      legend.position = "none",
      legend.text = element_text(size = 14),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    )
  
  # Panel B: NDVI ~ TDiffbin
  data_b <- NDVI_TDiffbin(data)
  data_b <- na.omit(data_b)
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
    
    aic_exp <- NA; r2_exp <- NA
    nls_exp <- tryCatch({
      nls(avg_value ~ a + b * exp(-c * x),
          data = sp_data, start = start_list_b, control = control_params_b)
    }, error = function(e) NULL)
    if (!is.null(nls_exp)) {
      aic_exp <- AIC(nls_exp)
      res <- resid(nls_exp)
      ss_res <- sum(res^2)
      ss_tot <- sum((sp_data$avg_value - mean(sp_data$avg_value))^2)
      r2_exp <- 1 - ss_res / ss_tot
    }
    
    sp_results <- data.frame(species = sp, Model = c("Linear", "Exponential"),
                             AIC = c(aic_linear, aic_exp), R2 = c(r2_linear, r2_exp))
    sp_results$y_label_pos <- sp_results$AIC / 2
    aic_results_b[[sp]] <- sp_results
  }
  aic_df_b <- do.call(rbind, aic_results_b)
  aic_df_b$species <- factor(aic_df_b$species, levels = species_order)
  
  p_b <- ggplot(aic_df_b, aes(x = species, y = AIC, fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(aes(label = round(R2, 2), y = y_label_pos),
              position = position_dodge(width = 0.9), vjust = 0.5, size = 8, color = "black") +
    labs(x = "", y = "AIC", title = "NDVI ~ TDiffbin", caption = "(b)") +
    scale_fill_manual(values = model_palette_shared) +
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
      legend.position = "none",
      legend.text = element_text(size = 14),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    )
  
  # Panel C: TDiff ~ PSIbin (with polynomial fits)
  data_c <- TDiff_PSIbin(data)
  data_c <- na.omit(data_c)
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
    
    sp_results <- data.frame(species = sp, Model = c("Linear", "Poly2", "Poly3"),
                             AIC = c(aic_linear, aic_poly2, aic_poly3),
                             R2 = c(r2_linear, r2_poly2, r2_poly3))
    sp_results$y_label_pos <- sp_results$AIC / 2
    aic_results_c[[sp]] <- sp_results
  }
  aic_df_c <- do.call(rbind, aic_results_c)
  aic_df_c$species <- factor(aic_df_c$species, levels = species_order)
  
  model_palette_c <- c("Linear" = "orange", "Poly2" = "dodgerblue", "Poly3" = "green4")
  
  p_c <- ggplot(aic_df_c, aes(x = species, y = AIC, fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(aes(label = round(R2, 2), y = y_label_pos),
              position = position_dodge(width = 0.9), vjust = 0.5, size = 8, color = "black") +
    labs(x = "", y = "AIC", title = "TDiff ~ PSIbin", caption = "(c)") +
    scale_fill_manual(values = model_palette_c) +
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
      legend.position = "none",
      legend.text = element_text(size = 14),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    )
  
  combined_top <- (p_a | p_b) + plot_layout(guides = "collect") +
    theme(legend.position = "top")
  combined_plot <- combined_top / p_c
  
  dir.create(dirname(save_combined_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = save_combined_fig, plot = combined_plot, width = 14, height = 12, dpi = 300)
  print(combined_plot)
}

plot_combined_AIC_R2(data, "results_rootzone/Figures/combined_AIC.png")

