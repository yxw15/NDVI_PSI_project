# =============================================================================
# NDVI and TDiff Analysis Script for PSI and TDiff Binning Methods
#
# This script performs the following steps:
# 1. Sets the working directory and sources custom functions.
# 2. Loads processed NDVI data (Quantiles and Proportions) for PSI and TDiff.
# 3. Generates plots for NDVI Quantiles and Proportions using:
#    - PSI binning method
#    - TDiff binning method
# 4. For each case, the script fits several models (Linear, Poly2, Poly3,
#    Exponential) to the data per species, computes the AIC values for model 
#    comparison, and saves the corresponding grouped bar plots.
#
# Required external R functions are loaded from:
# - "R/functions_PSI.R"
# - "R/functions_TDiff.R"
#
# Output plots are saved to the "results/key_displays/" directory.
# =============================================================================

# ---------------------------
# Set working directory and load functions
# ---------------------------
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

# Source custom functions for PSI and TDiff data processing
source("R/functions_PSI.R")
source("R/functions_TDiff.R")

# =============================================================================
# SECTION 1: NDVI Quantiles - PSI Analysis
# =============================================================================

# Load pre-processed NDVI Quantiles data (contains both PSI and TDiff data)
load("results/Data/All_Species_Quantiles_PSI_TDiff.RData")

# Generate NDVI Quantiles plots for PSI using an exponential model
plot_NDVI_Q_PSIbin_exp_slope(all_results_df, 
                             "results/key_displays/NDVI_Q_PSIbin_exp_coeff.png",
                             "results/key_displays/NDVI_Q_PSIbin_exp_slope.png")

# -----------------------------------------------------------------------------
# Function: plot_NDVI_Q_PSIbin_AIC
#
# This function:
# - Processes data with the custom NDVI_PSIbin() function.
# - Fits four models (Linear, Poly2, Poly3, Exponential) per species.
# - Computes the AIC values for each model.
# - Plots a grouped bar chart comparing AIC values and saves the plot.
# -----------------------------------------------------------------------------
plot_NDVI_Q_PSIbin_AIC <- function(data, save_aic_fig) {
  
  # Process data using the custom PSI binning function and remove missing values
  data <- NDVI_PSIbin(data)
  data <- na.omit(data)
  
  # Load necessary libraries
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(purrr)
  
  # Define response variable and species order
  value_col <- "avg_value"
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  data$species <- factor(data$species, levels = species_order)
  
  # Define a color palette for species (if needed later)
  cb_palette <- c("Oak"   = "#E69F00",   # Orange
                  "Beech" = "#0072B2",   # Deep blue
                  "Spruce"= "#009E73",   # Bluish-green
                  "Pine"  = "#F0E442")   # Yellow
  
  # Create a positive soil water potential variable (x = -bin_median)
  data <- data %>% mutate(x = -bin_median)
  
  # Clean data: remove rows with missing or non-finite values for avg_value or x
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  
  # Define starting values and control parameters for the exponential (nls) model
  start_list <- list(a = 5, b = 3, c = 0.001)
  control_params <- nls.control(maxiter = 200, minFactor = 1e-4)
  
  # Initialize list to store AIC results per species and model
  aic_results <- list()
  
  # Loop through each species to fit models and compute AIC values
  for (sp in levels(data_clean$species)) {
    sp_data <- data_clean %>% filter(species == sp)
    
    # Linear model: avg_value ~ x
    lm_linear <- lm(avg_value ~ x, data = sp_data)
    aic_linear <- AIC(lm_linear)
    
    # Quadratic model (Poly2): avg_value ~ x + I(x^2)
    lm_poly2 <- lm(avg_value ~ x + I(x^2), data = sp_data)
    aic_poly2 <- AIC(lm_poly2)
    
    # Cubic model (Poly3): avg_value ~ x + I(x^2) + I(x^3)
    lm_poly3 <- lm(avg_value ~ x + I(x^2) + I(x^3), data = sp_data)
    aic_poly3 <- AIC(lm_poly3)
    
    # Exponential model: avg_value ~ a + b * exp(-c * x)
    aic_exp <- NA
    nls_exp <- tryCatch({
      nls(avg_value ~ a + b * exp(-c * x), 
          data = sp_data, 
          start = start_list, 
          control = control_params)
    }, error = function(e) { NULL })
    if (!is.null(nls_exp)) {
      aic_exp <- AIC(nls_exp)
    }
    
    # Save AIC values for current species
    sp_results <- data.frame(species = sp,
                             Model = c("Linear", "Poly2", "Poly3", "Exponential"),
                             AIC = c(aic_linear, aic_poly2, aic_poly3, aic_exp))
    aic_results[[sp]] <- sp_results
  }
  
  # Combine AIC results for all species into one data frame
  aic_df <- do.call(rbind, aic_results)
  aic_df$species <- factor(aic_df$species, levels = species_order)
  
  # Print AIC table to console
  print(aic_df)
  
  # Define a palette for the models
  model_palette <- c("Linear"      = "#E69F00",
                     "Poly2"       = "#0072B2",
                     "Poly3"       = "#009E73",
                     "Exponential" = "#F0E442")
  
  # Create grouped bar plot of AIC values for each model per species
  p_aic <- ggplot(aic_df, aes(x = species, y = AIC, fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(aes(label = round(AIC, 2)), 
              position = position_dodge(width = 0.9), 
              vjust = -0.5, size = 3) +
    labs(x = "Species", y = "AIC", 
         title = "Model Comparison via AIC for NDVI Quantiles - PSI") +
    scale_fill_manual(values = model_palette) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid = element_blank(),
      legend.position = "top"
    )
  
  print(p_aic)
  
  # Save the AIC comparison plot to file
  dir.create(dirname(save_aic_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = save_aic_fig, plot = p_aic, width = 10, height = 8, dpi = 300)
}

# Call the function to generate and save the AIC plot for NDVI Quantiles - PSI data
plot_NDVI_Q_PSIbin_AIC(all_results_df, "results/key_displays/NDVI_Q_PSIbin_AIC.png")

# =============================================================================
# SECTION 2: NDVI Quantiles - TDiff Analysis
# =============================================================================

# Generate NDVI Quantiles plots for TDiff using a quadratic (Poly2) model
plot_Quantiles_TDiff_poly_2_slope_coeff(all_results_df, 
                                        "results/key_displays/NDVI_Q_TDiffbin_poly2_coeff.png",
                                        "results/key_displays/NDVI_Q_TDiffbin_poly2_slope.png")

# -----------------------------------------------------------------------------
# Function: plot_NDVI_Q_TDiffbin_AIC
#
# This function processes TDiff data using NDVI_TDiffbin(), fits four models 
# per species, calculates their AIC values, and plots the grouped bar chart.
# -----------------------------------------------------------------------------
plot_NDVI_Q_TDiffbin_AIC <- function(data, save_aic_fig) {
  
  # Process data with the custom TDiff binning function and remove missing values
  data <- NDVI_TDiffbin(data)
  data <- na.omit(data)
  
  # Load required libraries
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(purrr)
  
  # Define response variable and species order
  value_col <- "avg_value"
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  data$species <- factor(data$species, levels = species_order)
  
  # For TDiff analysis, use bin_median directly as x
  data <- data %>% mutate(x = bin_median)
  
  # Clean data: remove rows with missing or non-finite values
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  
  # Define starting values and control parameters for the exponential model
  start_list <- list(a = 5, b = 3, c = 0.001)
  control_params <- nls.control(maxiter = 200, minFactor = 1e-4)
  
  # Initialize list to store AIC results
  aic_results <- list()
  
  # Loop over each species to fit models and calculate AIC
  for (sp in levels(data_clean$species)) {
    sp_data <- data_clean %>% filter(species == sp)
    
    # Linear model
    lm_linear <- lm(avg_value ~ x, data = sp_data)
    aic_linear <- AIC(lm_linear)
    
    # Quadratic model (Poly2)
    lm_poly2 <- lm(avg_value ~ x + I(x^2), data = sp_data)
    aic_poly2 <- AIC(lm_poly2)
    
    # Cubic model (Poly3)
    lm_poly3 <- lm(avg_value ~ x + I(x^2) + I(x^3), data = sp_data)
    aic_poly3 <- AIC(lm_poly3)
    
    # Exponential model
    aic_exp <- NA
    nls_exp <- tryCatch({
      nls(avg_value ~ a + b * exp(-c * x), 
          data = sp_data, 
          start = start_list, 
          control = control_params)
    }, error = function(e) { NULL })
    if (!is.null(nls_exp)) {
      aic_exp <- AIC(nls_exp)
    }
    
    # Save the AIC values for the current species
    sp_results <- data.frame(species = sp,
                             Model = c("Linear", "Poly2", "Poly3", "Exponential"),
                             AIC = c(aic_linear, aic_poly2, aic_poly3, aic_exp))
    aic_results[[sp]] <- sp_results
  }
  
  # Combine all species AIC results
  aic_df <- do.call(rbind, aic_results)
  aic_df$species <- factor(aic_df$species, levels = species_order)
  
  # Define color palette for models
  model_palette <- c("Linear"      = "#E69F00",
                     "Poly2"       = "#0072B2",
                     "Poly3"       = "#009E73",
                     "Exponential" = "#F0E442")
  
  # Create the grouped bar plot for AIC values
  p_aic <- ggplot(aic_df, aes(x = species, y = AIC, fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
    labs(x = "Species", y = "AIC", 
         title = "Model Comparison via AIC for NDVI Quantiles - TDiff") +
    scale_fill_manual(values = model_palette) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid = element_blank(),
      legend.position = "top"
    )
  
  print(p_aic)
  
  # Save the plot as a PNG file
  dir.create(dirname(save_aic_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = save_aic_fig, plot = p_aic, width = 10, height = 8, dpi = 300)
}

# Call the function for NDVI Quantiles - TDiff AIC plot
plot_NDVI_Q_TDiffbin_AIC(all_results_df, "results/key_displays/NDVI_Q_TDiffbin_AIC.png")

# =============================================================================
# SECTION 3: NDVI Proportions Analysis - PSI
# =============================================================================

# Load pre-processed NDVI Proportions data
load("results/Data/All_Species_Proportions_PSI_TDiff.RData")

# Generate NDVI Proportions plots for PSI using an exponential model
plot_NDVI_PDM_PSIbin_exp_slope(all_results_df, 
                               "results/key_displays/NDVI_PDM_PSIbin_exp_coeff.png",
                               "results/key_displays/NDVI_PDM_PSIbin_exp_slope.png")

# -----------------------------------------------------------------------------
# Function: plot_NDVI_PDM_PSIbin_AIC
#
# This function:
# - Processes NDVI Proportions data using NDVI_PSIbin().
# - Fits four models per species and calculates AIC values.
# - Creates and saves a grouped bar plot for model comparison.
# -----------------------------------------------------------------------------
plot_NDVI_PDM_PSIbin_AIC <- function(data, save_aic_fig) {
  
  # Process data using the custom PSI binning function and remove missing rows
  data <- NDVI_PSIbin(data)
  data <- na.omit(data)
  
  # Load required libraries
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(purrr)
  
  # Define response variable and species order
  value_col <- "avg_value"
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  data$species <- factor(data$species, levels = species_order)
  
  # Define species color palette
  cb_palette <- c("Oak"   = "#E69F00",
                  "Beech" = "#0072B2",
                  "Spruce"= "#009E73",
                  "Pine"  = "#F0E442")
  
  # Create a positive soil water potential variable (x = -bin_median)
  data <- data %>% mutate(x = -bin_median)
  
  # Clean data: remove rows with missing or non-finite values
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  
  # Define starting values and control parameters for the exponential model
  start_list <- list(a = 5, b = 3, c = 0.001)
  control_params <- nls.control(maxiter = 200, minFactor = 1e-4)
  
  # Initialize list to store AIC results
  aic_results <- list()
  
  # Loop through each species to fit the models and calculate AIC values
  for (sp in levels(data_clean$species)) {
    sp_data <- data_clean %>% filter(species == sp)
    
    # Linear model
    lm_linear <- lm(avg_value ~ x, data = sp_data)
    aic_linear <- AIC(lm_linear)
    
    # Quadratic model (Poly2)
    lm_poly2 <- lm(avg_value ~ x + I(x^2), data = sp_data)
    aic_poly2 <- AIC(lm_poly2)
    
    # Cubic model (Poly3)
    lm_poly3 <- lm(avg_value ~ x + I(x^2) + I(x^3), data = sp_data)
    aic_poly3 <- AIC(lm_poly3)
    
    # Exponential model
    aic_exp <- NA
    nls_exp <- tryCatch({
      nls(avg_value ~ a + b * exp(-c * x), 
          data = sp_data, 
          start = start_list, 
          control = control_params)
    }, error = function(e) { NULL })
    if (!is.null(nls_exp)) {
      aic_exp <- AIC(nls_exp)
    }
    
    # Store the AIC values for current species
    sp_results <- data.frame(species = sp,
                             Model = c("Linear", "Poly2", "Poly3", "Exponential"),
                             AIC = c(aic_linear, aic_poly2, aic_poly3, aic_exp))
    aic_results[[sp]] <- sp_results
  }
  
  # Combine all species AIC results into one data frame
  aic_df <- do.call(rbind, aic_results)
  aic_df$species <- factor(aic_df$species, levels = species_order)
  
  # Define model color palette
  model_palette <- c("Linear"   = "#E69F00",
                     "Poly2" = "#0072B2",
                     "Poly3"= "#009E73",
                     "Exponential"  = "#F0E442")
  
  # Create grouped bar plot for the AIC comparison
  p_aic <- ggplot(aic_df, aes(x = species, y = AIC, fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
    labs(x = "Species", y = "AIC", 
         title = "Model Comparison via AIC for NDVI Proportions - PSI") +
    scale_fill_manual(values = model_palette) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid = element_blank(),
      legend.position = "top"
    )
  
  print(p_aic)
  
  # Save the plot to file
  dir.create(dirname(save_aic_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = save_aic_fig, plot = p_aic, width = 10, height = 8, dpi = 300)
}

# Call the function to generate and save the AIC plot for NDVI Proportions - PSI data
plot_NDVI_PDM_PSIbin_AIC(all_results_df, "results/key_displays/NDVI_PDM_PSIbin_AIC.png")

# =============================================================================
# SECTION 4: NDVI Proportions Analysis - TDiff
# =============================================================================

# Generate NDVI Proportions plots for TDiff using a quadratic (Poly2) model
plot_Proportions_TDiff_poly_2_slope_coeff(all_results_df, 
                                          "results/key_displays/NDVI_PDM_TDiffbin_poly2_coeff.png",
                                          "results/key_displays/NDVI_PDM_TDiffbin_poly2_slope.png")

# -----------------------------------------------------------------------------
# Function: plot_NDVI_PDM_TDiffbin_AIC
#
# This function processes NDVI Proportions data for TDiff, fits the four models,
# calculates their AIC values per species, and creates a grouped bar plot.
# -----------------------------------------------------------------------------
plot_NDVI_PDM_TDiffbin_AIC <- function(data, save_aic_fig) {
  
  # Process data using the custom TDiff binning function and remove missing values
  data <- NDVI_TDiffbin(data)
  data <- na.omit(data)
  
  # Load required libraries
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(purrr)
  
  # Define response variable and species order
  value_col <- "avg_value"
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  data$species <- factor(data$species, levels = species_order)
  
  # Use bin_median directly as x for TDiff analysis
  data <- data %>% mutate(x = bin_median)
  
  # Clean data: remove rows with missing or non-finite values
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  
  # Define starting values and control parameters for the exponential model
  start_list <- list(a = 5, b = 3, c = 0.001)
  control_params <- nls.control(maxiter = 200, minFactor = 1e-4)
  
  # Initialize list to store AIC results
  aic_results <- list()
  
  # Loop through each species to fit models and compute AIC values
  for (sp in levels(data_clean$species)) {
    sp_data <- data_clean %>% filter(species == sp)
    
    # Linear model
    lm_linear <- lm(avg_value ~ x, data = sp_data)
    aic_linear <- AIC(lm_linear)
    
    # Quadratic model (Poly2)
    lm_poly2 <- lm(avg_value ~ x + I(x^2), data = sp_data)
    aic_poly2 <- AIC(lm_poly2)
    
    # Cubic model (Poly3)
    lm_poly3 <- lm(avg_value ~ x + I(x^2) + I(x^3), data = sp_data)
    aic_poly3 <- AIC(lm_poly3)
    
    # Exponential model
    aic_exp <- NA
    nls_exp <- tryCatch({
      nls(avg_value ~ a + b * exp(-c * x), 
          data = sp_data, 
          start = start_list, 
          control = control_params)
    }, error = function(e) { NULL })
    if (!is.null(nls_exp)) {
      aic_exp <- AIC(nls_exp)
    }
    
    # Save AIC values for the species
    sp_results <- data.frame(species = sp,
                             Model = c("Linear", "Poly2", "Poly3", "Exponential"),
                             AIC = c(aic_linear, aic_poly2, aic_poly3, aic_exp))
    aic_results[[sp]] <- sp_results
  }
  
  # Combine AIC results for all species
  aic_df <- do.call(rbind, aic_results)
  aic_df$species <- factor(aic_df$species, levels = species_order)
  
  # Define the color palette for models
  model_palette <- c("Linear"      = "#E69F00",
                     "Poly2"       = "#0072B2",
                     "Poly3"       = "#009E73",
                     "Exponential" = "#F0E442")
  
  # Create grouped bar plot of AIC values
  p_aic <- ggplot(aic_df, aes(x = species, y = AIC, fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
    labs(x = "Species", y = "AIC", 
         title = "Model Comparison via AIC for NDVI Proportions - TDiff") +
    scale_fill_manual(values = model_palette) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid = element_blank(),
      legend.position = "top"
    )
  
  print(p_aic)
  
  # Save the plot as a PNG file
  dir.create(dirname(save_aic_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = save_aic_fig, plot = p_aic, width = 10, height = 8, dpi = 300)
}

# Call the function for NDVI Proportions - TDiff AIC plot
plot_NDVI_PDM_TDiffbin_AIC(all_results_df, "results/key_displays/NDVI_PDM_TDiff_AIC.png")

# =============================================================================
# SECTION 5: TDiff - PSI Analysis
# =============================================================================

# Generate TDiff PSI plots using a cubic (Poly3) model for slope coefficients
plot_TDiff_PSIbin_poly_3_slope(all_results_df, 
                               "results/key_displays/NDVI_TDiff_PSIbin_poly3_coeff.png",
                               "results/key_displays/NDVI_TDiff_PSIbin_poly3_slope.png")

# -----------------------------------------------------------------------------
# Function: plot_TDiff_PSIbin_AIC
#
# This function processes TDiff PSI data using TDiff_PSIbin(), fits several 
# models (Linear, Poly2, Poly3, Exponential) for each species, calculates AIC, 
# and generates a grouped bar plot for model comparison.
# -----------------------------------------------------------------------------
plot_TDiff_PSIbin_AIC <- function(data, save_aic_fig) {
  
  # Process data using the custom TDiff PSI binning function and remove missing values
  data <- TDiff_PSIbin(data)
  data <- na.omit(data)
  
  # Load required libraries
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(purrr)
  
  # Define response variable (average transpiration deficit) and species order
  value_col <- "avg_transpiration_deficit"
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  data$species <- factor(data$species, levels = species_order)
  
  # For TDiff, use bin_median directly as x
  data <- data %>% mutate(x = bin_median)
  
  # Clean data: remove rows with missing or non-finite values
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  
  # Define starting values and control parameters for the exponential model
  start_list <- list(a = 5, b = 3, c = 0.001)
  control_params <- nls.control(maxiter = 200, minFactor = 1e-4)
  
  # Initialize list to store AIC results per species
  aic_results <- list()
  
  # Loop over each species to fit models and calculate AIC
  for (sp in levels(data_clean$species)) {
    sp_data <- data_clean %>% filter(species == sp)
    
    # Linear model
    lm_linear <- lm(avg_transpiration_deficit ~ x, data = sp_data)
    aic_linear <- AIC(lm_linear)
    
    # Quadratic model (Poly2)
    lm_poly2 <- lm(avg_transpiration_deficit ~ x + I(x^2), data = sp_data)
    aic_poly2 <- AIC(lm_poly2)
    
    # Cubic model (Poly3)
    lm_poly3 <- lm(avg_transpiration_deficit ~ x + I(x^2) + I(x^3), data = sp_data)
    aic_poly3 <- AIC(lm_poly3)
    
    # Exponential model
    aic_exp <- NA
    nls_exp <- tryCatch({
      nls(avg_transpiration_deficit ~ a + b * exp(-c * x),
          data = sp_data,
          start = start_list,
          control = control_params)
    }, error = function(e) { NULL })
    if (!is.null(nls_exp)) {
      aic_exp <- AIC(nls_exp)
    }
    
    # Store AIC values for current species
    sp_results <- data.frame(species = sp,
                             Model = c("Linear", "Poly2", "Poly3", "Exponential"),
                             AIC = c(aic_linear, aic_poly2, aic_poly3, aic_exp))
    aic_results[[sp]] <- sp_results
  }
  
  # Combine AIC results for all species
  aic_df <- do.call(rbind, aic_results)
  aic_df$species <- factor(aic_df$species, levels = species_order)
  
  # Define color palette for models
  model_palette <- c("Linear"      = "#E69F00",
                     "Poly2"       = "#0072B2",
                     "Poly3"       = "#009E73",
                     "Exponential" = "#F0E442")
  
  # Create grouped bar plot for the AIC comparison
  p_aic <- ggplot(aic_df, aes(x = species, y = AIC, fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
    labs(x = "Species", y = "AIC",
         title = "Model Comparison via AIC for TDiff PSI Data") +
    scale_fill_manual(values = model_palette) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid = element_blank(),
      legend.position = "top"
    )
  
  print(p_aic)
  
  # Save the AIC plot to file
  dir.create(dirname(save_aic_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = save_aic_fig, plot = p_aic, width = 10, height = 8, dpi = 300)
}

# Call the function for TDiff PSI AIC plot
plot_TDiff_PSIbin_AIC(all_results_df, "results/key_displays/TDiff_PSIbin_AIC.png")
