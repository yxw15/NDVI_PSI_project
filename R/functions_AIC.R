setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

source("R/functions_PSI.R")
source("R/functions_TDiff.R")

plot_NDVI_Q_PSIbin_AIC <- function(data, save_aic_fig) {
  
  # Process data with your custom NDVI_PSIbin function and remove missing values
  data <- NDVI_PSIbin(data)
  data <- na.omit(data)
  
  # Load required libraries
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(purrr)
  
  # Identify the value column and order species; define the color palette
  value_col <- "avg_value"
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  data$species <- factor(data$species, levels = species_order)
  
  # Define a palette for species (for consistency, if needed later)
  cb_palette <- c("Oak"   = "#E69F00",   # Orange
                  "Beech" = "#0072B2",   # Deep blue
                  "Spruce"= "#009E73",   # Bluish-green
                  "Pine"  = "#F0E442")   # Yellow
  
  # Create a positive soil water potential variable (x = -bin_median)
  data <- data %>% mutate(x = -bin_median)
  
  # Clean data: remove rows with missing or non-finite values
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  
  # Define starting values and control parameters for the exponential model (nls)
  start_list <- list(a = 5, b = 3, c = 0.001)
  control_params <- nls.control(maxiter = 200, minFactor = 1e-4)
  
  # Initialize an empty list to store AIC results for each species and model
  aic_results <- list()
  
  # Loop through each species to fit the models and compute AIC
  for (sp in levels(data_clean$species)) {
    sp_data <- data_clean %>% filter(species == sp)
    
    # Fit linear model: avg_value ~ x
    lm_linear <- lm(avg_value ~ x, data = sp_data)
    aic_linear <- AIC(lm_linear)
    
    # Fit quadratic (poly2) model: avg_value ~ x + I(x^2)
    lm_poly2 <- lm(avg_value ~ x + I(x^2), data = sp_data)
    aic_poly2 <- AIC(lm_poly2)
    
    # Fit cubic (poly3) model: avg_value ~ x + I(x^2) + I(x^3)
    lm_poly3 <- lm(avg_value ~ x + I(x^2) + I(x^3), data = sp_data)
    aic_poly3 <- AIC(lm_poly3)
    
    # Fit exponential model: avg_value ~ a + b * exp(-c * x)
    aic_exp <- NA
    nls_exp <- tryCatch({
      nls(avg_value ~ a + b * exp(-c * x), 
          data = sp_data, 
          start = start_list, 
          control = control_params)
    }, error = function(e) {
      NULL
    })
    if (!is.null(nls_exp)) {
      aic_exp <- AIC(nls_exp)
    }
    
    # Combine the results for this species
    sp_results <- data.frame(species = sp,
                             Model = c("Linear", "Poly2", "Poly3", "Exponential"),
                             AIC = c(aic_linear, aic_poly2, aic_poly3, aic_exp))
    aic_results[[sp]] <- sp_results
  }
  
  # Combine AIC results for all species into one data frame
  aic_df <- do.call(rbind, aic_results)
  aic_df$species <- factor(aic_df$species, levels = species_order)
  
  # Print the AIC values for each case to the console
  print(aic_df)
  
  # Define a distinct color palette for models
  model_palette <- c("Linear"      = "#E69F00",   # Orange
                     "Poly2"       = "#0072B2",   # Deep blue
                     "Poly3"       = "#009E73",   # Bluish-green
                     "Exponential" = "#F0E442")   # Yellow
  
  # Create a grouped bar plot of AIC values per species for each model
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
  
  print(p_aic)
  
  # Save the AIC comparison plot to file
  dir.create(dirname(save_aic_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = save_aic_fig, plot = p_aic, width = 10, height = 8, dpi = 300)
}

plot_NDVI_Q_TDiffbin_AIC <- function(data, save_aic_fig) {
  
  # Process data with your custom NDVI_TDiffbin function and remove missing rows
  data <- NDVI_TDiffbin(data)
  data <- na.omit(data)
  
  # Load required libraries
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(purrr)
  
  # Identify the value column and order species; define the species palette if needed
  value_col <- "avg_value"
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  data$species <- factor(data$species, levels = species_order)
  
  # For TDiff, we assume x is already positive, so simply use bin_median
  data <- data %>% mutate(x = bin_median)
  
  # Clean data: remove rows with missing or non-finite values in avg_value and x
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  
  # Define starting values and control parameters for the exponential model (nls)
  start_list <- list(a = 5, b = 3, c = 0.001)
  control_params <- nls.control(maxiter = 200, minFactor = 1e-4)
  
  # Initialize an empty list to store AIC results for each species and model
  aic_results <- list()
  
  # Loop through each species to fit the models and compute AIC
  for (sp in levels(data_clean$species)) {
    sp_data <- data_clean %>% filter(species == sp)
    
    # Fit linear model: avg_value ~ x
    lm_linear <- lm(avg_value ~ x, data = sp_data)
    aic_linear <- AIC(lm_linear)
    
    # Fit quadratic (poly2) model: avg_value ~ x + I(x^2)
    lm_poly2 <- lm(avg_value ~ x + I(x^2), data = sp_data)
    aic_poly2 <- AIC(lm_poly2)
    
    # Fit cubic (poly3) model: avg_value ~ x + I(x^2) + I(x^3)
    lm_poly3 <- lm(avg_value ~ x + I(x^2) + I(x^3), data = sp_data)
    aic_poly3 <- AIC(lm_poly3)
    
    # Fit exponential model: avg_value ~ a + b * exp(-c * x)
    aic_exp <- NA
    nls_exp <- tryCatch({
      nls(avg_value ~ a + b * exp(-c * x), 
          data = sp_data, 
          start = start_list, 
          control = control_params)
    }, error = function(e) {
      NULL
    })
    if (!is.null(nls_exp)) {
      aic_exp <- AIC(nls_exp)
    }
    
    # Combine the results for this species
    sp_results <- data.frame(species = sp,
                             Model = c("Linear", "Poly2", "Poly3", "Exponential"),
                             AIC = c(aic_linear, aic_poly2, aic_poly3, aic_exp))
    aic_results[[sp]] <- sp_results
  }
  
  # Combine AIC results for all species into one data frame
  aic_df <- do.call(rbind, aic_results)
  aic_df$species <- factor(aic_df$species, levels = species_order)
  
  # Define the model palette as specified
  model_palette <- c("Linear"      = "#E69F00",   # Orange
                     "Poly2"       = "#0072B2",   # Deep blue
                     "Poly3"       = "#009E73",   # Bluish-green
                     "Exponential" = "#F0E442")   # Yellow
  
  # Create a grouped bar plot of AIC values per species for each model
  p_aic <- ggplot(aic_df, aes(x = species, y = AIC, fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
    labs(x = "Species", y = "AIC", 
         title = "Model Comparison via AIC for NDVI Quantiles - TDiff") +
    scale_fill_manual(values = model_palette) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    )
  
  print(p_aic)
  
  # Save the AIC comparison plot to file
  dir.create(dirname(save_aic_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = save_aic_fig, plot = p_aic, width = 10, height = 8, dpi = 300)
}

plot_NDVI_PDM_PSIbin_AIC <- function(data, save_aic_fig) {
  
  # Process data with your custom NDVI_PSIbin function and remove missing rows
  data <- NDVI_PSIbin(data)
  data <- na.omit(data)
  
  # Load required libraries
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(purrr)
  
  # Identify the value column and order species; define the species color palette
  value_col <- "avg_value"
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  data$species <- factor(data$species, levels = species_order)
  
  cb_palette <- c("Oak"   = "#E69F00",   # Orange
                  "Beech" = "#0072B2",   # Deep blue
                  "Spruce"= "#009E73",   # Bluish-green
                  "Pine"  = "#F0E442")   # Yellow
  
  # Create a positive soil water potential variable (x = -bin_median)
  data <- data %>% mutate(x = -bin_median)
  
  # Clean data: remove rows with missing or non-finite values
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  
  # Define starting values and control parameters for the exponential model (nls)
  start_list <- list(a = 5, b = 3, c = 0.001)
  control_params <- nls.control(maxiter = 200, minFactor = 1e-4)
  
  # Initialize an empty list to store AIC results for each species and model
  aic_results <- list()
  
  # Loop through each species to fit the models and compute AIC
  for (sp in levels(data_clean$species)) {
    sp_data <- data_clean %>% filter(species == sp)
    
    # Fit linear model: avg_value ~ x
    lm_linear <- lm(avg_value ~ x, data = sp_data)
    aic_linear <- AIC(lm_linear)
    
    # Fit quadratic (poly2) model: avg_value ~ x + I(x^2)
    lm_poly2 <- lm(avg_value ~ x + I(x^2), data = sp_data)
    aic_poly2 <- AIC(lm_poly2)
    
    # Fit cubic (poly3) model: avg_value ~ x + I(x^2) + I(x^3)
    lm_poly3 <- lm(avg_value ~ x + I(x^2) + I(x^3), data = sp_data)
    aic_poly3 <- AIC(lm_poly3)
    
    # Fit exponential model: avg_value ~ a + b * exp(-c * x)
    aic_exp <- NA
    nls_exp <- tryCatch({
      nls(avg_value ~ a + b * exp(-c * x), 
          data = sp_data, 
          start = start_list, 
          control = control_params)
    }, error = function(e) {
      NULL
    })
    if (!is.null(nls_exp)) {
      aic_exp <- AIC(nls_exp)
    }
    
    # Combine the results for this species
    sp_results <- data.frame(species = sp,
                             Model = c("Linear", "Poly2", "Poly3", "Exponential"),
                             AIC = c(aic_linear, aic_poly2, aic_poly3, aic_exp))
    aic_results[[sp]] <- sp_results
  }
  
  # Combine AIC results for all species into one data frame
  aic_df <- do.call(rbind, aic_results)
  aic_df$species <- factor(aic_df$species, levels = species_order)
  
  # Define a distinct color palette for models
  model_palette <- c("Linear"   = "#E69F00",   # Orange
                     "Poly2" = "#0072B2",   # Deep blue
                     "Poly3"= "#009E73",   # Bluish-green
                     "Exponential"  = "#F0E442")   # Yellow
  
  # Create a grouped bar plot of AIC values per species for each model
  p_aic <- ggplot(aic_df, aes(x = species, y = AIC, fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
    labs(x = "Species", y = "AIC", 
         title = "Model Comparison via AIC for NDVI Proportions - PSI") +
    scale_fill_manual(values = model_palette) +
    theme_minimal() +
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
  
  print(p_aic)
  
  # Save the AIC comparison plot to file
  dir.create(dirname(save_aic_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = save_aic_fig, plot = p_aic, width = 10, height = 8, dpi = 300)
}

plot_NDVI_PDM_TDiffbin_AIC <- function(data, save_aic_fig) {
  
  # Process data with your custom NDVI_TDiffbin function and remove missing rows
  data <- NDVI_TDiffbin(data)
  data <- na.omit(data)
  
  # Load required libraries
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(purrr)
  
  # Identify the value column and order species; use the "avg_value" column
  value_col <- "avg_value"
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  data$species <- factor(data$species, levels = species_order)
  
  # For TDiff, assume x is already positive; set x as bin_median
  data <- data %>% mutate(x = bin_median)
  
  # Clean data: remove rows with missing values or non-finite x
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  
  # Define starting values and control parameters for the exponential model (nls)
  start_list <- list(a = 5, b = 3, c = 0.001)
  control_params <- nls.control(maxiter = 200, minFactor = 1e-4)
  
  # Initialize an empty list to store AIC results for each species and model
  aic_results <- list()
  
  # Loop through each species to fit models and compute AIC
  for (sp in levels(data_clean$species)) {
    sp_data <- data_clean %>% filter(species == sp)
    
    # Fit linear model: avg_value ~ x
    lm_linear <- lm(avg_value ~ x, data = sp_data)
    aic_linear <- AIC(lm_linear)
    
    # Fit quadratic (Poly2) model: avg_value ~ x + I(x^2)
    lm_poly2 <- lm(avg_value ~ x + I(x^2), data = sp_data)
    aic_poly2 <- AIC(lm_poly2)
    
    # Fit cubic (Poly3) model: avg_value ~ x + I(x^2) + I(x^3)
    lm_poly3 <- lm(avg_value ~ x + I(x^2) + I(x^3), data = sp_data)
    aic_poly3 <- AIC(lm_poly3)
    
    # Fit exponential model: avg_value ~ a + b * exp(-c * x)
    aic_exp <- NA
    nls_exp <- tryCatch({
      nls(avg_value ~ a + b * exp(-c * x), 
          data = sp_data, 
          start = start_list, 
          control = control_params)
    }, error = function(e) {
      NULL
    })
    if (!is.null(nls_exp)) {
      aic_exp <- AIC(nls_exp)
    }
    
    # Combine the results for this species
    sp_results <- data.frame(species = sp,
                             Model = c("Linear", "Poly2", "Poly3", "Exponential"),
                             AIC = c(aic_linear, aic_poly2, aic_poly3, aic_exp))
    aic_results[[sp]] <- sp_results
  }
  
  # Combine AIC results for all species into one data frame
  aic_df <- do.call(rbind, aic_results)
  aic_df$species <- factor(aic_df$species, levels = species_order)
  
  # Define the model color palette as specified
  model_palette <- c("Linear"      = "#E69F00",   # Orange
                     "Poly2"       = "#0072B2",   # Deep blue
                     "Poly3"       = "#009E73",   # Bluish-green
                     "Exponential" = "#F0E442")   # Yellow
  
  # Create a grouped bar plot of AIC values per species for each model
  p_aic <- ggplot(aic_df, aes(x = species, y = AIC, fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
    labs(x = "Species", y = "AIC", 
         title = "Model Comparison via AIC for NDVI Proportions - TDiff") +
    scale_fill_manual(values = model_palette) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    )
  
  print(p_aic)
  
  # Save the AIC comparison plot to file
  dir.create(dirname(save_aic_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = save_aic_fig, plot = p_aic, width = 10, height = 8, dpi = 300)
}

plot_TDiff_PSIbin_AIC <- function(data, save_aic_fig) {
  
  # Process the data with your custom TDiff_PSIbin function and remove missing rows
  data <- TDiff_PSIbin(data)
  data <- na.omit(data)
  
  # Load required libraries
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(purrr)
  
  # Identify the response column (average transpiration deficit) and order species.
  value_col <- "avg_transpiration_deficit"
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  data$species <- factor(data$species, levels = species_order)
  
  # For TDiff, we use the bin_median directly as x.
  data <- data %>% mutate(x = bin_median)
  
  # Clean data: remove rows with missing or non-finite values in the response and x
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  
  # Define starting values and control parameters for the exponential model (nls)
  start_list <- list(a = 5, b = 3, c = 0.001)
  control_params <- nls.control(maxiter = 200, minFactor = 1e-4)
  
  # Initialize an empty list to store AIC results for each species and model
  aic_results <- list()
  
  # Loop through each species to fit models and compute AIC
  for (sp in levels(data_clean$species)) {
    sp_data <- data_clean %>% filter(species == sp)
    
    # Fit a linear model: avg_transpiration_deficit ~ x
    lm_linear <- lm(avg_transpiration_deficit ~ x, data = sp_data)
    aic_linear <- AIC(lm_linear)
    
    # Fit a quadratic (Poly2) model: avg_transpiration_deficit ~ x + I(x^2)
    lm_poly2 <- lm(avg_transpiration_deficit ~ x + I(x^2), data = sp_data)
    aic_poly2 <- AIC(lm_poly2)
    
    # Fit a cubic (Poly3) model: avg_transpiration_deficit ~ x + I(x^2) + I(x^3)
    lm_poly3 <- lm(avg_transpiration_deficit ~ x + I(x^2) + I(x^3), data = sp_data)
    aic_poly3 <- AIC(lm_poly3)
    
    # Fit an exponential model: avg_transpiration_deficit ~ a + b * exp(-c * x)
    aic_exp <- NA
    nls_exp <- tryCatch({
      nls(avg_transpiration_deficit ~ a + b * exp(-c * x),
          data = sp_data,
          start = start_list,
          control = control_params)
    }, error = function(e) {
      NULL
    })
    if (!is.null(nls_exp)) {
      aic_exp <- AIC(nls_exp)
    }
    
    # Combine the results for this species
    sp_results <- data.frame(species = sp,
                             Model = c("Linear", "Poly2", "Poly3", "Exponential"),
                             AIC = c(aic_linear, aic_poly2, aic_poly3, aic_exp))
    aic_results[[sp]] <- sp_results
  }
  
  # Combine AIC results for all species into one data frame
  aic_df <- do.call(rbind, aic_results)
  aic_df$species <- factor(aic_df$species, levels = species_order)
  
  # Define the model color palette as specified.
  model_palette <- c("Linear"      = "#E69F00",   # Orange
                     "Poly2"       = "#0072B2",   # Deep blue
                     "Poly3"       = "#009E73",   # Bluish-green
                     "Exponential" = "#F0E442")   # Yellow
  
  # Create a grouped bar plot of AIC values per species for each model
  p_aic <- ggplot(aic_df, aes(x = species, y = AIC, fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
    labs(x = "Species", y = "AIC",
         title = "Model Comparison via AIC for TDiff PSI Data") +
    scale_fill_manual(values = model_palette) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    )
  
  print(p_aic)
  
  # Save the AIC comparison plot to file
  dir.create(dirname(save_aic_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = save_aic_fig, plot = p_aic, width = 10, height = 8, dpi = 300)
}

