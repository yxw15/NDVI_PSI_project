setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")
load("results/Data/AllSpecies_AllMonths_rootzone.RData")

library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(patchwork)

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

### S7 soil texture composition ###
source("R/plot_soil_map_MODIS.R")

plot_soil_fraction_composite <- function(
    species_order = c("Oak", "Beech", "Spruce", "Pine"),
    soil_map_dir = "soil_map",
    result_data_dir = "results_rootzone/Data",
    result_fig_dir = "results_rootzone/Figures",
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
    result_fig_dir = "results_rootzone/Figures",
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

### TDiff - PSI ###
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
                               "results_rootzone/Figures/TDiff_PSI_poly2_coeff.png",
                               "results_rootzone/Figures/TDiff_PSI_poly2_slope.png")

plot_TDiff_PSIbin_poly_3_slope(data, 
                               "results_rootzone/Figures/TDiff_PSI_poly3_coeff.png",
                               "results_rootzone/Figures/TDiff_PSI_poly3_slope.png")

### Parameterize analysis ###
