setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")
library(dplyr)
library(tidyverse)

### load data ###
# depth 150 #
load("results/Data/AllSpecies_AllMonths_depth150.RData")
combined_depth150 <- combined
combined_depth150 <- combined_depth150 %>% filter(Quantiles > 0)
combined_depth150 <- combined_depth150 %>% filter(month %in% c("May", "June", "July", "August"))
combined_depth150 <- combined_depth150 %>% filter(Quantiles > 0)
combined_depth150_July_August <- combined_depth150 %>% filter(month %in% c("July", "August"))

# depth 100 #
load("results/Data/AllSpecies_AllMonths_depth100.RData")
combined_depth100 <- combined
combined_depth100 <- combined_depth100 %>% filter(Quantiles > 0)
combined_depth100 <- combined_depth100 %>% filter(month %in% c("May", "June", "July", "August"))
combined_depth100 <- combined_depth100 %>% filter(Quantiles > 0)
combined_depth100_July_August <- combined_depth100 %>% filter(month %in% c("July", "August"))


# depth 50 #
load("results/Data/AllSpecies_AllMonths_depth50.RData")
combined_depth50 <- combined
combined_depth50 <- combined_depth50 %>% filter(Quantiles > 0)
combined_depth50 <- combined_depth50 %>% filter(month %in% c("May", "June", "July", "August"))
combined_depth50 <- combined_depth50 %>% filter(Quantiles > 0)
combined_depth50_July_August <- combined_depth50 %>% filter(month %in% c("July", "August"))

# depth all #
combined_depth_all <- rbind(combined_depth50, combined_depth100, combined_depth150)
combined_depth_all <- combined_depth_all %>% filter(Quantiles > 0)
combined_depth_all_July_August <- combined_depth_all %>% filter(month %in% c("July", "August"))

### S1 original ###
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

# depth 150 #
plot_NDVI_PSI_month_species_linear_orig(combined_depth150, 
                                        "results/Figures/depth150/NDVI_PSI_linear_slope_monthly_depth150.png",
                                        "results/Figures/depth150/NDVI_PSI_linear_coeff_monthly_depth150.png")

plot_NDVI_TDiff_month_species_linear_orig(combined_depth150, 
                                        "results/Figures/depth150/NDVI_TDiff_linear_slope_monthly_depth150.png",
                                        "results/Figures/depth150/NDVI_TDiff_linear_coeff_monthly_depth150.png")

# depth 100 #
plot_NDVI_PSI_month_species_linear_orig(combined_depth100, 
                                        "results/Figures/depth100/NDVI_PSI_linear_slope_monthly_depth100.png",
                                        "results/Figures/depth100/NDVI_PSI_linear_coeff_monthly_depth100.png")

plot_NDVI_TDiff_month_species_linear_orig(combined_depth100, 
                                          "results/Figures/depth100/NDVI_TDiff_linear_slope_monthly_depth100.png",
                                          "results/Figures/depth100/NDVI_TDiff_linear_coeff_monthly_depth100.png")

# depth 50 #
plot_NDVI_PSI_month_species_linear_orig(combined_depth50, 
                                        "results/Figures/depth50/NDVI_PSI_linear_slope_monthly_depth50.png",
                                        "results/Figures/depth50/NDVI_PSI_linear_coeff_monthly_depth50.png")

plot_NDVI_TDiff_month_species_linear_orig(combined_depth50, 
                                          "results/Figures/depth50/NDVI_TDiff_linear_slope_monthly_depth50.png",
                                          "results/Figures/depth50/NDVI_TDiff_linear_coeff_monthly_depth50.png")

# depth all #
plot_NDVI_PSI_month_species_linear_orig(combined_depth_all, 
                                        "results/Figures/depth_all/NDVI_PSI_linear_slope_monthly_depth_all.png",
                                        "results/Figures/depth_all/NDVI_PSI_linear_coeff_monthly_depth_all.png")

plot_NDVI_TDiff_month_species_linear_orig(combined_depth_all, 
                                          "results/Figures/depth_all/NDVI_TDiff_linear_slope_monthly_depth_all.png",
                                          "results/Figures/depth_all/NDVI_TDiff_linear_coeff_monthly_depth_all.png")


### density ###
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

# depth 150 #
plot_density_NDVI_PSI_linear(combined_depth150_July_August,
                             swp_bin_width = 50,
                             output_path = "results/Figures/depth150/NDVI_PSI_density_depth150.png")

plot_density_NDVI_TDiff_linear(combined_depth150_July_August,
                               tdiff_bin_width = 3,
                               output_path = "results/Figures/depth150/NDVI_Tdiff_density_depth150.png")

# depth 100 #
plot_density_NDVI_PSI_linear(combined_depth100_July_August,
                             swp_bin_width = 50,
                             output_path = "results/Figures/depth100/NDVI_PSI_density_depth100.png")

plot_density_NDVI_TDiff_linear(combined_depth100_July_August,
                               tdiff_bin_width = 3,
                               output_path = "results/Figures/depth100/NDVI_Tdiff_density_depth100.png")

# depth 50 #
plot_density_NDVI_PSI_linear(combined_depth50_July_August,
                             swp_bin_width = 50,
                             output_path = "results/Figures/depth50/NDVI_PSI_density_depth50.png")

plot_density_NDVI_TDiff_linear(combined_depth50_July_August,
                               tdiff_bin_width = 3,
                               output_path = "results/Figures/depth50/NDVI_Tdiff_density_depth50.png")

# depth all #
plot_density_NDVI_PSI_linear(combined_depth_all_July_August,
                             bin_width = 50,
                             output_path = "results/Figures/depth_all/NDVI_PSI_density_depth_all.png")

plot_density_NDVI_TDiff_linear(combined_depth_all_July_August,
                               tdiff_bin_width = 3,
                               output_path = "results/Figures/depth_all/NDVI_Tdiff_density_depth_all.png")

### distribution PSI and TDiff ### 
plot_distribution_PSI_TDiff <- function(
    data,
    var,            # column name to average & plot, e.g. "soil_water_potential"
    legend_title,   # legend title, e.g. "soil water potential (kPa)"
    output_file,    # where to save the PNG
    limits = NULL   # optional numeric vector c(min, max) for the color scale
) {
  
  library(dplyr)
  library(ggplot2)
  library(grid)          # for unit()
  library(sf)            # if boundary_germany is an sf object
  library(rnaturalearth)
  library(rnaturalearthdata)
  
  # fetch Germany border
  boundary_germany <- ne_countries(
    scale      = "medium", 
    country    = "Germany", 
    returnclass= "sf"
  )
  
  # compute per‐species mean
  dat_mean <- data %>%
    group_by(x, y, species) %>%
    summarise(
      mean_val = mean(.data[[var]], na.rm = TRUE),
      .groups  = "drop"
    )
  
  # enforce species order
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  dat_mean$species <- factor(dat_mean$species, levels = species_order)
  
  # choose palette (reverse only for SWP)
  base_pal <- c("blue","dodgerblue","cyan","yellow","orange","red")
  pal      <- if (var == "soil_water_potential") rev(base_pal) else base_pal
  
  # build plot
  p <- ggplot(dat_mean, aes(x=x, y=y, color=mean_val)) +
    geom_point(size = 0.5) +
    geom_sf(data = boundary_germany, fill = NA, color = "black", inherit.aes = FALSE) +
    scale_color_gradientn(
      colours = pal,
      name    = legend_title,
      limits  = limits  # here we use the shared limits
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
  
  # save
  ggsave(
    filename = output_file,
    plot     = p,
    width    = 14,
    height   = 6,
    dpi      = 300
  )
  
  invisible(p)
}

compute_swp_means <- function(df) {
  df %>%
    group_by(x,y,species) %>%
    summarise(mean_val = mean(soil_water_potential, na.rm=TRUE), .groups="drop") %>%
    pull(mean_val)
}

all_vals <- c(
  compute_swp_means(combined_depth150_July_August),
  compute_swp_means(combined_depth100_July_August),
  compute_swp_means(combined_depth50_July_August)
)
common_limits <- range(all_vals, na.rm = TRUE)

plot_distribution_PSI_TDiff(
  data         = combined_depth150_July_August,
  var          = "soil_water_potential",
  legend_title = "soil water potential (kPa)",
  output_file  = "results/Figures/depth150/mean_soil_water_potential_July_August_depth150.png",
  limits       = common_limits
)

plot_distribution_PSI_TDiff(
  data         = combined_depth100_July_August,
  var          = "soil_water_potential",
  legend_title = "soil water potential (kPa)",
  output_file  = "results/Figures/depth100/mean_soil_water_potential_July_August_depth100.png",
  limits       = common_limits
)

plot_distribution_PSI_TDiff(
  data         = combined_depth50_July_August,
  var          = "soil_water_potential",
  legend_title = "soil water potential (kPa)",
  output_file  = "results/Figures/depth50/mean_soil_water_potential_July_August_depth50.png",
  limits       = common_limits
)

# Transpiration deficit map (forward palette)
plot_distribution_PSI_TDiff(
  data         = combined_depth150_July_August,
  var          = "transpiration_deficit",
  legend_title = "transpiration deficit (mm)",
  output_file  = "results/Figures/depth150/mean_transpiration_deficit_July_August.png"
)

### distribution Temp and VPD ###
plot_distribution_Temp_VPD <- function(
    df,
    var,            # column name in df, e.g. "temp" or "vpd"
    legend_title,   # e.g. "temperature (°C)" or "VPD (kPa)"
    output_file,    # full path to save PNG
    limits = NULL   # optional two‐element vector to lock scale
) {
  
  library(ggplot2)
  library(sf)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(grid)  
  
  # fetch Germany boundary
  boundary_germany <- ne_countries(
    scale       = "medium",
    country     = "Germany",
    returnclass = "sf"
  )
  
  species_order = c("Oak","Beech","Spruce","Pine")
  palette = c("blue","dodgerblue","cyan","yellow","orange","red")
  
  # ensure species factor
  df$species <- factor(df$species, levels = species_order)
  
  p <- ggplot(df, aes_string(x = "x", y = "y", color = var)) +
    geom_point(size = 0.5) +
    geom_sf(data = boundary_germany, fill = NA, color = "black", inherit.aes = FALSE) +
    scale_color_gradientn(
      colours = palette,
      name    = legend_title,
      limits  = limits
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
  
  ggsave(
    filename = output_file,
    plot     = p,
    width    = 14,
    height   = 6,
    dpi      = 300
  )
  
  print(p)
}


plot_mean_map(
  df           = all_df,
  var          = "temp",
  legend_title = "temperature (°C)",
  output_file  = "results/Figures/environment/mean_temperature_July_August.png"
  # limits     = temp_limits
)

plot_mean_map(
  df           = all_df,
  var          = "vpd",
  legend_title = "VPD (kPa)",
  output_file  = "results/Figures/environment/mean_vpd_July_August.png"
  # limits     = vpd_limits
)

### mean time series ###
plot_mean_time_series_PSI_TDiff <- function(
    df50,                # data.frame for 50 cm depth (must have species, soil_water_potential)
    df100,               # data.frame for 100 cm depth
    df150,               # data.frame for 150 cm depth (needs both soil_water_potential and transpiration_deficit)
    output_file
) {
  
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  
  
  cb_palette <- c(
    "Oak"    = "#E69F00",
    "Beech"  = "#0072B2",
    "Spruce" = "#009E73",
    "Pine"   = "#F0E442"
  )
  
  species_order = c("Oak","Beech","Spruce","Pine")
  
  # helper to summarise SWP at a given depth
  summarize_swp <- function(df, depth_cm) {
    df %>%
      group_by(species) %>%
      summarise(Mean_Value = mean(soil_water_potential, na.rm = TRUE), .groups="drop") %>%
      mutate(Parameter = paste0("Mean SWP at ", depth_cm, " cm"))
  }
  
  # and for transpiration deficit (only from df150)
  summarize_tdef <- function(df, depth_cm) {
    df %>%
      group_by(species) %>%
      summarise(Mean_Value = mean(transpiration_deficit, na.rm = TRUE), .groups="drop") %>%
      mutate(Parameter = paste0("Mean Transp. Deficit at ", depth_cm, " cm"))
  }
  
  # 1. build each summary
  swp50_df   <- summarize_swp(df50,  50)
  swp100_df  <- summarize_swp(df100, 100)
  swp150_df  <- summarize_swp(df150, 150)
  tdef150_df <- summarize_tdef(df150,150)
  
  # 2. bind & factor
  summary_all <- bind_rows(swp50_df, swp100_df, swp150_df, tdef150_df) %>%
    mutate(
      species = str_to_title(species),
      species = factor(species, levels = species_order),
      Parameter = factor(
        Parameter,
        levels = c(
          "Mean SWP at 50 cm",
          "Mean SWP at 100 cm",
          "Mean SWP at 150 cm",
          "Mean Transp. Deficit at 150 cm"
        ),
        labels = c(
          expression(Phi[~soil]~(50*cm)),
          expression(Phi[~soil]~(100*cm)),
          expression(Phi[~soil]~(150*cm)),
          expression(T[diff])
        )
      )
    )
  
  # 3. plot
  p <- ggplot(summary_all, aes(x = species, y = Mean_Value, fill = species)) +
    geom_col(width = 0.7, show.legend = FALSE) +
    facet_wrap(
      ~ Parameter, nrow = 1, scales = "free_y",
      labeller = label_parsed        # ← Tell ggplot to parse those expressions
    ) +
    scale_fill_manual(values = cb_palette) +
    labs(
      x        = NULL,
      y        = "mean value",
      title    = "",
      subtitle = ""
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
      # axis.text.x = element_text(angle = 45, hjust = 1),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.text = element_text(size = 14),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 12)
    )
  
  # 4. save
  ggsave(
    filename = output_file,
    plot     = p,
    width    = 16,
    height   = 4,
    dpi      = 300
  )
  
  print(p)
}


plot_mean_time_series_PSI_TDiff(
  df50         = combined_depth50_July_August,
  df100        = combined_depth100_July_August,
  df150        = combined_depth150_July_August,
  output_file  = "results/Figures/depth_all/time_series_PSI_TDiff.png")

### AIC ###
plot_combined_AIC_R2 <- function(data, depths = c(50, 100, 150), save_path) {
  # Load required libraries
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(patchwork)
  
  # Color-blind–friendly palette
  cb_palette <- c(
    linear      = "#E69F00",
    exponential = "#0072B2",
    poly2       = "#009E73",
    poly3       = "#F0E442"
  )
  
  # Base theme with large tags nudged down & right
  base_theme <- theme_minimal() +
    theme(
      axis.text.x        = element_text(angle = 0, hjust = 0.5, size = 14),
      plot.background    = element_rect(fill = "white", color = "white"),
      panel.background   = element_rect(fill = "white"),
      legend.background  = element_rect(fill = "white", color = "white"),
      legend.text        = element_text(color = "black", size = 14),
      plot.title         = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title         = element_text(face = "bold", size = 16),
      axis.text          = element_text(color = "black", size = 14),
      panel.border       = element_blank(),
      panel.grid.major   = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.tag           = element_text(size = 16, face = "bold", color = "black"),
      plot.tag.position  = c(0.03, 0.96)
    )
  
  # Helper: compute AIC & R² for NDVI ~(PSI or TD)
  compute_aic_ndvi <- function(df, type = c("PSI", "TD"), species_order) {
    type <- match.arg(type)
    if (type == "PSI") {
      d      <- NDVI_PSIbin(df)   %>% na.omit() %>% mutate(x = -bin_median)
      starts <- list(a = 5, b = 3,   c = 0.001)
      ctrl   <- nls.control(maxiter = 200, minFactor = 1e-4)
    } else {
      d      <- NDVI_TDiffbin(df) %>% na.omit() %>% mutate(x =  bin_median)
      starts <- list(a = 5, b = 7,   c = 0.04)
      ctrl   <- nls.control(maxiter = 1200, minFactor = 1e-9)
    }
    d$species <- factor(d$species, levels = species_order)
    map_df(species_order, function(sp) {
      spd <- filter(d, species == sp)
      # linear model
      m1 <- lm(avg_value ~ x, data = spd)
      a1 <- AIC(m1); r1 <- summary(m1)$r.squared
      # exponential model
      expo <- tryCatch(
        nls(avg_value ~ a + b * exp(-c * x),
            data    = spd,
            start   = starts,
            control = ctrl),
        error = function(e) NULL
      )
      if (!is.null(expo)) {
        a2 <- AIC(expo)
        r2 <- 1 - sum(resid(expo)^2) / sum((spd$avg_value - mean(spd$avg_value))^2)
      } else {
        a2 <- NA; r2 <- NA
      }
      tibble(
        species     = factor(sp, levels = species_order),
        Model       = c("linear", "exponential"),
        AIC         = c(a1, a2),
        R2          = c(r1, r2),
        y_label_pos = c(a1/2, a2/2)
      )
    })
  }
  
  # Helper: compute AIC & R² for transpiration deficit ~ PSI (linear, poly2, poly3)
  compute_aic_td_psi <- function(df, species_order) {
    d <- TDiff_PSIbin(df) %>% na.omit() %>% mutate(x = bin_median)
    d$species <- factor(d$species, levels = species_order)
    map_df(species_order, function(sp) {
      spd <- filter(d, species == sp)
      m1 <- lm(avg_transpiration_deficit ~ x,                    data = spd)
      m2 <- lm(avg_transpiration_deficit ~ x + I(x^2),           data = spd)
      m3 <- lm(avg_transpiration_deficit ~ x + I(x^2) + I(x^3),  data = spd)
      tibble(
        species     = factor(sp, levels = species_order),
        Model       = c("linear","poly2","poly3"),
        AIC         = c(AIC(m1), AIC(m2), AIC(m3)),
        R2          = c(summary(m1)$r.squared,
                        summary(m2)$r.squared,
                        summary(m3)$r.squared),
        y_label_pos = c(AIC(m1)/2, AIC(m2)/2, AIC(m3)/2)
      )
    })
  }
  
  # Main plotting
  species_order <- c("Oak","Beech","Spruce","Pine")
  
  # (a)-(c): NDVI ~ PSI
  p_abc <- map2(depths, letters[1:3], ~{
    df_d   <- filter(data, depth == .x)
    aic_df <- compute_aic_ndvi(df_d, "PSI", species_order)
    ggplot(aic_df, aes(x = species, y = AIC, fill = Model)) +
      geom_col(position = position_dodge(0.9)) +
      geom_text(aes(label = round(R2, 2), y = y_label_pos),
                position = position_dodge(0.9), vjust = 0.5, size = 5) +
      scale_fill_manual(values = cb_palette, name = "model (a)-(c)") +
      labs(title = paste0("depth ", .x, " cm"),
           tag   = paste0("(", .y, ")"),
           x     = NULL, y = "AIC") +
      base_theme +
      theme(legend.position = "none",
            plot.tag.position = c(0.12, 0.74))
  })
  
  # (g): NDVI ~ transpiration deficit at last depth
  df150  <- filter(data, depth == depths[length(depths)])
  aic_td <- compute_aic_ndvi(df150, "TD", species_order)
  p_g    <- ggplot(aic_td, aes(x = species, y = AIC, fill = Model)) +
    geom_col(position = position_dodge(0.9)) +
    geom_text(aes(label = round(R2, 2), y = y_label_pos),
              position = position_dodge(0.9), vjust = 0.5, size = 5) +
    scale_fill_manual(values = cb_palette, name = "model (g)") +
    labs(tag   = "(g)", x     = NULL, y = "AIC") +
    base_theme +
    theme(legend.position = "none",
          plot.tag.position = c(0.04, 0.95))
  
  # (d)-(f): transpiration deficit ~ PSI
  p_def <- map2(depths, letters[4:6], ~{
    df_d   <- filter(data, depth == .x)
    aic_df <- compute_aic_td_psi(df_d, species_order)
    ggplot(aic_df, aes(x = species, y = AIC, fill = Model)) +
      geom_col(position = position_dodge(0.9)) +
      geom_text(aes(label = round(R2, 2), y = y_label_pos),
                position = position_dodge(0.9), vjust = 0.5, size = 5) +
      scale_fill_manual(values = cb_palette, name = "model (d)-(f)") +
      labs(tag   = paste0("(", .y, ")"), x = NULL, y = "AIC") +
      base_theme +
      theme(legend.position = "none",
            plot.tag.position = c(0.12, 0.95))
  })
  
  # Combine panels
  combined <-
    (p_abc[[1]] | p_abc[[2]] | p_abc[[3]]) /
    (p_def[[1]] | p_def[[2]] | p_def[[3]]) /
    p_g +
    plot_layout(guides = "collect") &
    theme(legend.position = "top")
  
  # Save & display
  dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
  ggsave(save_path, combined, width = 20, height = 10, dpi = 300)
  print(combined)
}

plot_combined_AIC_R2(data = combined_depth_all_July_August,
                     save_path = "results/Figures/depth_all/combined_AIC_R2.png")

