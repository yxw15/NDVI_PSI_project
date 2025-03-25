# Load required libraries
library(tidyverse)
library(spdep)
library(ggplot2)

# Set working directory and load the data
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")
load("results/Data/All_Species_Quantiles_PSI_TDiff.RData")

# Subset the data for species "Pine"
pine_data <- subset(all_results_df, species == "Pine")

# Estimate the parameter 'a' for Pine (using the minimum Quantiles)
a_est <- min(pine_data$Quantiles)

# Transform the dependent variable to linearize the exponential model:
# Y = log(Quantiles - a_est)
# Ensure that Quantiles - a_est > 0 for all observations
pine_data <- pine_data %>%
  mutate(Y = log(Quantiles - a_est))

# Construct the spatial weights matrix using 4-nearest neighbors based on (x, y) coordinates
coords <- cbind(pine_data$x, pine_data$y)
knn_neighbors <- knearneigh(coords, k = 4)
nb <- knn2nb(knn_neighbors)
listw <- nb2listw(nb)

# Fit the spatial lag model on the transformed variable Y
# The model is: Y = ρ * W * Y + β₀ + β₁ * soil_water_potential + ε
model_pine <- lagsarlm(Y ~ soil_water_potential, data = pine_data, listw = listw)

# Print summary of the spatial lag model
summary(model_pine)

# Generate predictions on the transformed scale over a range of soil_water_potential values
soil_seq <- seq(min(pine_data$soil_water_potential),
                max(pine_data$soil_water_potential), length.out = 100)
coefs <- model_pine$coefficients  # β₀ (intercept) and β₁ (slope)
pred_Y <- coefs[1] + coefs[2] * soil_seq

# Back-transform the predictions to the original scale:
# Fitted Quantiles = a_est + exp(predicted Y)
pred_Quantiles <- a_est + exp(pred_Y)

pred_df <- data.frame(soil_water_potential = soil_seq,
                      fitted_Quantiles = pred_Quantiles)

# Plot the original data and fitted exponential curve
ggplot(pine_data, aes(x = soil_water_potential, y = Quantiles)) +
  geom_point(alpha = 0.7) +
  geom_line(data = pred_df, aes(x = soil_water_potential, y = fitted_Quantiles),
            color = "blue", size = 1) +
  labs(title = "Fitted Exponential Function for Pine",
       x = "Soil Water Potential",
       y = "Quantiles") +
  theme_minimal()
