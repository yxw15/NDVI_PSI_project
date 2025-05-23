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
  cb_palette <- c("Oak"   = "#E69F00",
                  "Beech" = "#0072B2",
                  "Spruce"= "#009E73",
                  "Pine"  = "#F0E442")
  
  # Create a positive soil water potential variable (x = -bin_median)
  data <- data %>% mutate(x = -bin_median)
  
  # Clean data: remove rows with missing or non-finite values
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  
  # Define starting values and control parameters for the exponential model (nls)
  start_list <- list(a = 5, b = 3, c = 0.001)
  control_params <- nls.control(maxiter = 200, minFactor = 1e-4)
  
  # Initialize an empty list to store AIC results for each species and model
  aic_results <- list()
  
  # Loop through each species to fit the models and compute AIC and annotations
  for (sp in levels(data_clean$species)) {
    sp_data <- data_clean %>% filter(species == sp)
    
    # Fit models
    lm_linear <- lm(avg_value ~ x, data = sp_data)
    aic_linear <- AIC(lm_linear)
    r2_linear <- summary(lm_linear)$r.squared
    
    lm_poly2 <- lm(avg_value ~ x + I(x^2), data = sp_data)
    aic_poly2 <- AIC(lm_poly2)
    r2_poly2 <- summary(lm_poly2)$r.squared
    
    lm_poly3 <- lm(avg_value ~ x + I(x^2) + I(x^3), data = sp_data)
    aic_poly3 <- AIC(lm_poly3)
    r2_poly3 <- summary(lm_poly3)$r.squared
    
    # Exponential model
    aic_exp <- NA
    r2_exp <- NA
    nls_exp <- tryCatch({
      nls(avg_value ~ a + b * exp(-c * x), 
          data = sp_data, 
          start = start_list, 
          control = control_params)
    }, error = function(e) NULL)
    
    if (!is.null(nls_exp)) {
      aic_exp <- AIC(nls_exp)
      res <- resid(nls_exp)
      ss_res <- sum(res^2)
      ss_tot <- sum((sp_data[[value_col]] - mean(sp_data[[value_col]]))^2)
      r2_exp <- 1 - ss_res / ss_tot
    }
    
    # Create annotation labels: show only R²
    label_fun <- function(r2) {
      if (is.na(r2)) return("NA")
      return(as.character(round(r2, 2)))
    }
    
    labels <- c(
      label_fun(r2_linear),
      label_fun(r2_poly2),
      label_fun(r2_poly3),
      label_fun(r2_exp)
    )
    
    # Combine AICs and labels
    sp_results <- data.frame(
      species = sp,
      Model = c("Linear", "Poly2", "Poly3", "Exponential"),
      AIC = c(aic_linear, aic_poly2, aic_poly3, aic_exp),
      Label = labels
    )
    
    # Set Y position to center of bar
    sp_results$y_label_pos <- sp_results$AIC / 2
    
    aic_results[[sp]] <- sp_results
  }
  
  # Combine results
  aic_df <- do.call(rbind, aic_results)
  aic_df$species <- factor(aic_df$species, levels = species_order)
  
  # Print AIC table
  print(aic_df)
  
  # Model color palette
  model_palette <- c("Linear"      = "#E69F00",
                     "Poly2"       = "#0072B2",
                     "Poly3"       = "#009E73",
                     "Exponential" = "#F0E442")
  
  # Plot
  p_aic <- ggplot(aic_df, aes(x = species, y = AIC, fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(aes(label = Label, y = y_label_pos),
              position = position_dodge(width = 0.9),
              vjust = 0.5, size = 3, color = "black") +
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
  
  # Save
  dir.create(dirname(save_aic_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = save_aic_fig, plot = p_aic, width = 10, height = 8, dpi = 300)
}

plot_NDVI_Q_PSIbin_AIC_exp_linear <- function(data, save_aic_fig) {
  
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
  cb_palette <- c("Oak"   = "#E69F00",
                  "Beech" = "#0072B2",
                  "Spruce"= "#009E73",
                  "Pine"  = "#F0E442")
  
  # Create a positive soil water potential variable (x = -bin_median)
  data <- data %>% mutate(x = -bin_median)
  
  # Clean data: remove rows with missing or non-finite values
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  
  # Define starting values and control parameters for the exponential model (nls)
  start_list <- list(a = 5, b = 3, c = 0.001)
  control_params <- nls.control(maxiter = 200, minFactor = 1e-4)
  
  # Initialize an empty list to store AIC results for each species and model
  aic_results <- list()
  
  # Loop through each species to fit the models and compute AIC and annotations
  for (sp in levels(data_clean$species)) {
    sp_data <- data_clean %>% filter(species == sp)
    
    # Fit the linear model
    lm_linear <- lm(avg_value ~ x, data = sp_data)
    aic_linear <- AIC(lm_linear)
    r2_linear <- summary(lm_linear)$r.squared
    
    # Fit the exponential model using nls
    aic_exp <- NA
    r2_exp <- NA
    nls_exp <- tryCatch({
      nls(avg_value ~ a + b * exp(-c * x), 
          data = sp_data, 
          start = start_list, 
          control = control_params)
    }, error = function(e) NULL)
    
    if (!is.null(nls_exp)) {
      aic_exp <- AIC(nls_exp)
      res <- resid(nls_exp)
      ss_res <- sum(res^2)
      ss_tot <- sum((sp_data[[value_col]] - mean(sp_data[[value_col]]))^2)
      r2_exp <- 1 - ss_res / ss_tot
    }
    
    # Create annotation labels: show only R² values
    label_fun <- function(r2) {
      if (is.na(r2)) return("NA")
      return(as.character(round(r2, 2)))
    }
    
    labels <- c(
      label_fun(r2_linear),
      label_fun(r2_exp)
    )
    
    # Combine AICs and labels for the two models
    sp_results <- data.frame(
      species = sp,
      Model = c("Linear", "Exponential"),
      AIC = c(aic_linear, aic_exp),
      Label = labels
    )
    
    # Set Y position to the center of the bar for text annotation
    sp_results$y_label_pos <- sp_results$AIC / 2
    
    aic_results[[sp]] <- sp_results
  }
  
  # Combine results from all species
  aic_df <- do.call(rbind, aic_results)
  aic_df$species <- factor(aic_df$species, levels = species_order)
  
  # Print AIC table for inspection
  print(aic_df)
  
  # Define the model color palette for the two models
  model_palette <- c("Linear" = "orange",
                     "Exponential" = "dodgerblue")
  
  # Plot the AIC comparisons
  p_aic <- ggplot(aic_df, aes(x = species, y = AIC, fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(aes(label = Label, y = y_label_pos),
              position = position_dodge(width = 0.9),
              vjust = 0.5, size = 6, color = "black") +
    labs(x = "Species", y = "AIC", 
         title = "Model Comparison via AIC for NDVI Quantiles - PSI") +
    scale_fill_manual(values = model_palette) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 22, face = "bold", color = "black"),
      plot.subtitle = element_text(hjust = 0.5, size = 18),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(color = "black", size = 14),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.text = element_text(size = 14),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 16)
    )
  
  print(p_aic)
  
  # Save the resulting plot to the specified file
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
  
  # Define starting values and control parameters for the exponential model 
  # (updated to match the exponential function used in your first function)
  start_list <- list(a = 5, b = 7, c = 0.04)
  control_params <- nls.control(maxiter = 1200, minFactor = 1e-09)
  
  # Initialize an empty list to store AIC results for each species and model
  aic_results <- list()
  
  # Loop through each species to fit the models and compute AIC and R² values
  for (sp in levels(data_clean$species)) {
    sp_data <- data_clean %>% filter(species == sp)
    
    # Fit linear model: avg_value ~ x
    lm_linear <- lm(avg_value ~ x, data = sp_data)
    aic_linear <- AIC(lm_linear)
    r2_linear <- summary(lm_linear)$r.squared
    
    # Fit quadratic (poly2) model: avg_value ~ x + I(x^2)
    lm_poly2 <- lm(avg_value ~ x + I(x^2), data = sp_data)
    aic_poly2 <- AIC(lm_poly2)
    r2_poly2 <- summary(lm_poly2)$r.squared
    
    # Fit cubic (poly3) model: avg_value ~ x + I(x^2) + I(x^3)
    lm_poly3 <- lm(avg_value ~ x + I(x^2) + I(x^3), data = sp_data)
    aic_poly3 <- AIC(lm_poly3)
    r2_poly3 <- summary(lm_poly3)$r.squared
    
    # Fit exponential model: avg_value ~ a + b * exp(-c * x)
    aic_exp <- NA
    r2_exp <- NA
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
      res <- resid(nls_exp)
      ss_res <- sum(res^2)
      ss_tot <- sum((sp_data[[value_col]] - mean(sp_data[[value_col]]))^2)
      r2_exp <- 1 - ss_res / ss_tot
    }
    
    # Create annotation labels (rounded R² values)
    label_fun <- function(r2) {
      if (is.na(r2)) return("NA")
      return(as.character(round(r2, 2)))
    }
    
    labels <- c(
      label_fun(r2_linear),
      label_fun(r2_poly2),
      label_fun(r2_poly3),
      label_fun(r2_exp)
    )
    
    # Combine results for this species
    sp_results <- data.frame(
      species = sp,
      Model = c("Linear", "Poly2", "Poly3", "Exponential"),
      AIC = c(aic_linear, aic_poly2, aic_poly3, aic_exp),
      Label = labels
    )
    sp_results$y_label_pos <- sp_results$AIC / 2
    
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
  
  # Create a grouped bar plot of AIC values per species for each model with R² annotations
  p_aic <- ggplot(aic_df, aes(x = species, y = AIC, fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
    geom_text(aes(label = Label, y = y_label_pos),
              position = position_dodge(width = 0.8),
              vjust = 0.5, size = 3, color = "black") +
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

plot_NDVI_Q_TDiffbin_AIC_exp_linear <- function(data, save_aic_fig) {
  
  # Process data with NDVI_TDiffbin and remove missing rows
  data <- NDVI_TDiffbin(data)
  data <- na.omit(data)
  
  # Load required libraries
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(purrr)
  
  # Identify the value column and order species
  value_col <- "avg_value"
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  data$species <- factor(data$species, levels = species_order)
  
  # For TDiff, use bin_median as x
  data <- data %>% mutate(x = bin_median)
  
  # Clean data: remove rows with missing or non-finite values in avg_value and x
  data_clean <- data %>% filter(!is.na(.data[[value_col]]), is.finite(x))
  
  # Define starting values and control parameters for the exponential model 
  start_list <- list(a = 5, b = 7, c = 0.04)
  control_params <- nls.control(maxiter = 1200, minFactor = 1e-09)
  
  # Initialize an empty list to store AIC results for each species and model
  aic_results <- list()
  
  # Loop through each species
  for (sp in levels(data_clean$species)) {
    sp_data <- data_clean %>% filter(species == sp)
    
    # Fit linear model: avg_value ~ x
    lm_linear <- lm(avg_value ~ x, data = sp_data)
    aic_linear <- AIC(lm_linear)
    r2_linear <- summary(lm_linear)$r.squared
    
    # Fit exponential model: avg_value ~ a + b * exp(-c * x)
    aic_exp <- NA
    r2_exp <- NA
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
      res <- resid(nls_exp)
      ss_res <- sum(res^2)
      ss_tot <- sum((sp_data[[value_col]] - mean(sp_data[[value_col]]))^2)
      r2_exp <- 1 - ss_res / ss_tot
    }
    
    # Create annotation labels (rounded R² values)
    label_fun <- function(r2) {
      if (is.na(r2)) return("NA")
      return(as.character(round(r2, 2)))
    }
    
    labels <- c(
      label_fun(r2_linear),
      label_fun(r2_exp)
    )
    
    # Combine results for this species
    sp_results <- data.frame(
      species = sp,
      Model = c("Linear", "Exponential"),
      AIC = c(aic_linear, aic_exp),
      Label = labels
    )
    sp_results$y_label_pos <- sp_results$AIC / 2
    
    aic_results[[sp]] <- sp_results
  }
  
  # Combine AIC results for all species into one data frame
  aic_df <- do.call(rbind, aic_results)
  aic_df$species <- factor(aic_df$species, levels = species_order)
  
  # Define the model palette for the two models
  model_palette <- c("Linear" = "orange", "Exponential" = "dodgerblue")
  
  # Create a grouped bar plot of AIC values per species for each model with R² annotations
  p_aic <- ggplot(aic_df, aes(x = species, y = AIC, fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    geom_text(aes(label = Label, y = y_label_pos),
              position = position_dodge(width = 0.8),
              vjust = 0.5, size = 6, color = "black") +
    labs(x = "Species", y = "AIC", 
         title = "Model Comparison via AIC for NDVI Quantiles - TDiff") +
    scale_fill_manual(values = model_palette) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      plot.title = element_text(hjust = 0.5, size = 22, face = "bold", color = "black"),
      plot.subtitle = element_text(hjust = 0.5, size = 18),
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(color = "black", size = 14),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.text = element_text(size = 14),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
      strip.text = element_text(face = "bold", size = 16)
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
  
  # Initialize an empty list to store AIC and R² results for each species and model
  aic_results <- list()
  
  # Loop through each species to fit models and compute AIC and R²
  for (sp in levels(data_clean$species)) {
    sp_data <- data_clean %>% filter(species == sp)
    
    # Fit a linear model: avg_transpiration_deficit ~ x
    lm_linear <- lm(avg_transpiration_deficit ~ x, data = sp_data)
    aic_linear <- AIC(lm_linear)
    r2_linear <- summary(lm_linear)$r.squared
    
    # Fit a quadratic (Poly2) model: avg_transpiration_deficit ~ x + I(x^2)
    lm_poly2 <- lm(avg_transpiration_deficit ~ x + I(x^2), data = sp_data)
    aic_poly2 <- AIC(lm_poly2)
    r2_poly2 <- summary(lm_poly2)$r.squared
    
    # Fit a cubic (Poly3) model: avg_transpiration_deficit ~ x + I(x^2) + I(x^3)
    lm_poly3 <- lm(avg_transpiration_deficit ~ x + I(x^2) + I(x^3), data = sp_data)
    aic_poly3 <- AIC(lm_poly3)
    r2_poly3 <- summary(lm_poly3)$r.squared
    
    # Fit an exponential model: avg_transpiration_deficit ~ a + b * exp(-c * x)
    aic_exp <- NA
    r2_exp <- NA
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
      # Compute R² for the exponential model
      y <- sp_data$avg_transpiration_deficit
      y_pred <- predict(nls_exp)
      r2_exp <- 1 - sum((y - y_pred)^2) / sum((y - mean(y))^2)
    }
    
    # Combine the results for this species into a data frame
    sp_results <- data.frame(species = sp,
                             Model = c("Linear", "Poly2", "Poly3", "Exponential"),
                             AIC = c(aic_linear, aic_poly2, aic_poly3, aic_exp),
                             R2 = c(r2_linear, r2_poly2, r2_poly3, r2_exp))
    aic_results[[sp]] <- sp_results
  }
  
  # Combine AIC and R² results for all species into one data frame
  aic_df <- do.call(rbind, aic_results)
  aic_df$species <- factor(aic_df$species, levels = species_order)
  
  # Define the model color palette as specified.
  model_palette <- c("Linear"      = "#E69F00",   # Orange
                     "Poly2"       = "#0072B2",   # Deep blue
                     "Poly3"       = "#009E73",   # Bluish-green
                     "Exponential" = "#F0E442")   # Yellow
  
  # Create a grouped bar plot of AIC values per species for each model
  p_aic <- ggplot(aic_df, aes(x = species, y = AIC, fill = Model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    labs(x = "Species", y = "AIC",
         title = "Model Comparison via AIC for TDiff PSI Data") +
    scale_fill_manual(values = model_palette) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
      axis.text.y = element_text(size = 16),
      axis.title = element_text(face = "bold", size = 20),
      plot.title = element_text(hjust = 0.5, size = 24, face = "bold", color = "black"),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 18, face = "bold"),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(fill = "white", color = "white"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    ) +
    # Add R² labels at the center of each bar (using AIC/2 as the y position)
    geom_text(aes(label = paste0(round(R2, 2)), y = AIC/2),
              position = position_dodge(width = 0.9),
              color = "black",
              size = 5)
  
  print(p_aic)
  
  # Save the AIC comparison plot to file
  dir.create(dirname(save_aic_fig), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = save_aic_fig, plot = p_aic, width = 10, height = 8, dpi = 300)
}

