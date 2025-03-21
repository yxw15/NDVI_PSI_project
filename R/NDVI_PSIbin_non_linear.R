# ================================
# SET WORKING DIRECTORY & LOAD DATA
# ================================
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")
load("results/Data/NDVI_PSI_TDiff_species_df.RData")

# Assuming you have functions NDVI_PSIbin() and NDVI_TDiffbin_df() defined elsewhere
NDVI_PSIbin_df <- NDVI_PSIbin(NDVI_PSI_TDiff_species)
NDVI_PSIbin_df <- na.omit(NDVI_PSIbin_df)

# ================================
# LOAD REQUIRED LIBRARIES
# ================================
library(ggplot2)
library(nlme)      # for nonlinear modeling (nlsList)
library(dplyr)     # for data manipulation
library(tibble)    # for converting rownames to a column
library(patchwork) # for combining plots
library(purrr)     # for mapping functions

# ================================
# DATA PREPARATION
# ================================
# Order species and define a color palette
NDVI_PSIbin_df$species <- factor(NDVI_PSIbin_df$species, levels = c("Oak", "Beech", "Spruce", "Pine"))
cb_palette <- c("Oak" = "#E69F00", "Beech" = "#56B4E9", 
                "Spruce" = "#009E73", "Pine" = "#F0E442")

# Create a positive version of soil water potential (x = -bin_median)
NDVI_PSIbin_df <- NDVI_PSIbin_df %>% mutate(x = -bin_median)

# Clean data: remove rows with missing or non-finite values
data_clean <- NDVI_PSIbin_df %>% filter(!is.na(avg_quantiles), !is.na(bin_median), is.finite(x))

# ================================
# NONLINEAR MODELING PER SPECIES
# ================================
# Fit the model: avg_quantiles = a + b * exp(-c * x)
start_list <- list(a = 5, b = 3, c = 0.001)
control_params <- nls.control(maxiter = 200, minFactor = 1e-4)

nls_models <- nlsList(avg_quantiles ~ a + b * exp(-c * x) | species,
                      data = data_clean,
                      start = start_list,
                      control = control_params)
print(summary(nls_models))

# ================================
# EXTRACT COEFFICIENTS & CALCULATE DERIVED VALUES
# ================================
# Convert model coefficients to a data frame.
coef_df <- as.data.frame(coef(nls_models), optional = TRUE) %>% 
  rownames_to_column(var = "species") %>%
  filter(!is.na(a))

# For each species, compute the reference x (mean of x) and derive slope and intercept.
ref_vals <- data_clean %>%
  group_by(species) %>%
  summarize(x_ref = mean(x, na.rm = TRUE), .groups = "drop")

coef_df <- left_join(coef_df, ref_vals, by = "species") %>%
  mutate(slope_at_x_ref = -b * c * exp(-c * x_ref),
         intercept_at_x0 = a + b)

# Calculate x50 and slope50 (the x and slope where avg_quantiles reaches 50% of max)
threshold <- 0.5 * max(data_clean$avg_quantiles, na.rm = TRUE)
coef_df <- coef_df %>%
  mutate(x50 = ifelse((threshold - a) > 0 & b > 0,
                      -log((threshold - a)/b) / c,
                      NA),
         slope50 = ifelse(!is.na(x50),
                          -b * c * exp(-c * x50),
                          NA))

# ================================
# GENERATE PREDICTIONS FOR PANEL A
# ================================
pred_list <- data_clean %>%
  group_by(species) %>%
  do({
    sp <- unique(.$species)
    x_seq <- seq(min(.$x, na.rm = TRUE), max(.$x, na.rm = TRUE), length.out = 100)
    # Predict using the species-specific model
    sp_model <- nls_models[[as.character(sp)]]
    pred <- predict(sp_model, newdata = data.frame(x = x_seq))
    data.frame(x = x_seq, pred = pred)
  })
pred_all <- bind_rows(pred_list)

# ================================
# PANEL A: Observed Data and Fitted Curves
# ================================
# Reverse the x-axis to display original negative values.
x_scale <- scale_x_continuous(trans = "reverse", labels = function(x) -x)
min_quantile <- min(data_clean$avg_quantiles, na.rm = TRUE)
max_quantile <- max(data_clean$avg_quantiles, na.rm = TRUE)

p_combined <- ggplot() +
  geom_point(data = data_clean, aes(x = x, y = avg_quantiles, color = species)) +
  geom_line(data = pred_all, aes(x = x, y = pred, color = species), size = 1) +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "black", size = 1) +
  scale_color_manual(values = cb_palette) +
  x_scale +
  coord_cartesian(ylim = c(min_quantile, max_quantile)) +
  labs(x = "Soil Water Potential", y = "Quantiles") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
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

# ================================
# PANEL B: Bar Plot with p-value and R² labels
# ================================
# Compute additional statistics (slope SE, p-value, and R²) for each species.
additional_stats <- map_dfr(as.character(unique(coef_df$species)), function(sp) {
  if (!(sp %in% names(nls_models))) {
    warning("No model found for species: ", sp)
    return(NULL)
  }
  mod <- nls_models[[sp]]
  mod_sum <- summary(mod)
  se_b <- mod_sum$parameters["b", "Std. Error"]
  se_c <- mod_sum$parameters["c", "Std. Error"]
  cov_bc <- vcov(mod)["b", "c"]
  
  # Calculate R² using the species data
  df_sp <- data_clean %>% filter(species == sp) %>% drop_na(avg_quantiles, x)
  fitted_vals <- predict(mod, newdata = df_sp)
  SS_res <- sum((df_sp$avg_quantiles - fitted_vals)^2)
  SS_tot <- sum((df_sp$avg_quantiles - mean(df_sp$avg_quantiles))^2)
  r_squared <- 1 - SS_res / SS_tot
  
  # Retrieve species-specific coefficients for x50 and slope50.
  row_sp <- coef_df %>% filter(species == sp)
  x50_val <- row_sp$x50
  b_val <- row_sp$b
  c_val <- row_sp$c
  
  if (is.na(x50_val)) {
    slope_se <- NA
    p_val <- NA
  } else {
    # Compute partial derivatives of slope50 = -b * c * exp(-c * x50)
    dfdB <- c_val * exp(-c_val * x50_val)
    dfdC <- b_val * exp(-c_val * x50_val) * (c_val * x50_val - 1)
    var_slope <- (dfdB^2) * se_b^2 + (dfdC^2) * se_c^2 + 2 * dfdB * dfdC * cov_bc
    slope_se <- sqrt(var_slope)
    t_val <- row_sp$slope50 / slope_se
    p_val <- 2 * (1 - pnorm(abs(t_val)))  # normal approximation
  }
  
  tibble(species = sp, p_value = p_val, r_squared = r_squared)
})

# Merge the additional statistics with the coefficient data.
slope_data <- coef_df %>%
  left_join(additional_stats, by = "species") %>%
  mutate(slope_abs = abs(slope50)) %>%
  arrange(slope_abs) %>%
  mutate(species = factor(species, levels = unique(species)),
         label_text = ifelse(p_value < 0.005,
                             sprintf("p < 0.005\nR² = %.2f", r_squared),
                             sprintf("p = %.3f\nR² = %.2f", p_value, r_squared)))

p_bar <- ggplot(slope_data, aes(x = species, y = slope_abs, fill = species)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = cb_palette) +
  labs(x = "Species", y = "Absolute Slope at 50% Max Quantile") +
  coord_flip() +
  geom_text(aes(x = species, y = slope_abs/2, label = label_text),
            size = 4, color = "black") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    plot.background = element_rect(fill = "white", color = "white"),
    panel.background = element_rect(fill = "white"),
    legend.background = element_rect(fill = "white", color = "white"),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    panel.border = element_rect(color = NA, fill = NA, linewidth = 0),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    strip.text = element_text(face = "bold", size = 12)
  )

# -------------------------
# COMBINE PANEL A and PANEL B SIDE-BY-SIDE
# -------------------------
final_plot <- p_combined + p_bar + plot_layout(widths = c(2, 1))
print(final_plot)

# Ensure the output directory exists before saving the plot.
dir.create("results/Figures", recursive = TRUE, showWarnings = FALSE)
output_file <- "results/Figures/NDVI_PSI_non_linear.png"
ggsave(filename = output_file, plot = final_plot, width = 10, height = 6, dpi = 300)
