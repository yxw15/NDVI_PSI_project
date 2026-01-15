# ==========================================
# Complete TDiff x PSI plotting script (4-panel) with messages
# ==========================================

# --------------------------
# 0. Libraries
# --------------------------
message("Loading libraries...")
library(terra)
library(dplyr)
library(ggplot2)
library(scales)

# --------------------------
# 1. Parameters
# --------------------------
message("Setting parameters...")
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

species <- c("Oak", "Beech")
months  <- c("July")
years   <- 2018
base_dir <- "results_monthly_rootzone"

# --------------------------
# 2. Load species masks
# --------------------------
message("Loading species masks...")
species_masks <- lapply(species, function(sp) {
  message("  - Loading mask for ", sp)
  rast(file.path("species_map_MODIS", paste0(sp, ".tif")))
})
names(species_masks) <- species

# --------------------------
# 3. Build panel dataframe function
# --------------------------
message("Defining function to build panel dataframe with month/year...")

build_panel_df <- function(mask_sp, month_i, year_i) {
  message("Processing panel: ", mask_sp, " | Month: ", month_i, " | Year: ", year_i)
  
  mask <- species_masks[[mask_sp]]
  
  # --- TDiff: only panel species ---
  tdiff_path <- file.path(base_dir, month_i, mask_sp, paste0("tdiff_", year_i, ".tif"))
  tdiff <- rast(tdiff_path) |> project(mask) |> mask(mask)
  tdiff_df <- as.data.frame(tdiff, xy = TRUE, na.rm = TRUE)
  names(tdiff_df)[3] <- "transpiration_deficit"
  tdiff_df$mask_species <- mask_sp
  tdiff_df$month <- month_i
  tdiff_df$year  <- year_i
  
  # --- PSI: all species masked by panel ---
  psi_list <- lapply(species, function(psi_sp) {
    message("  - Masking PSI species ", psi_sp, " to panel ", mask_sp)
    psi_path <- file.path(base_dir, month_i, psi_sp, paste0("psi_", year_i, ".tif"))
    psi <- rast(psi_path) |> project(mask) |> mask(mask)
    df <- as.data.frame(psi, xy = TRUE, na.rm = TRUE)
    names(df)[3] <- "soil_water_potential"
    df$PSI_species <- psi_sp
    df$mask_species <- mask_sp
    df$month <- month_i
    df$year  <- year_i
    df
  })
  psi_df <- bind_rows(psi_list)
  
  # Join TDiff and PSI by x, y, mask_species, month, year
  left_join(tdiff_df, psi_df, by = c("x", "y", "mask_species", "month", "year"))
}

# --------------------------
# 4. Generate full dataset across months, years, species
# --------------------------
message("Building full dataset across all months, years, and species...")
full_df_list <- list()
for (m in months) {
  for (y in years) {
    for (sp in species) {
      key <- paste(sp, m, y, sep="_")
      full_df_list[[key]] <- build_panel_df(sp, m, y)
    }
  }
}
full_df <- bind_rows(full_df_list)

bin_summary <- full_df %>%
  filter(!is.na(transpiration_deficit), !is.na(soil_water_potential)) %>%
  mutate(PSI_bin = cut(soil_water_potential,
                       breaks = seq(floor(min(soil_water_potential)), ceiling(max(soil_water_potential)), by = 50),
                       include.lowest = TRUE, right = FALSE)) %>%
  group_by(mask_species, PSI_species, PSI_bin) %>%
  summarise(
    avg_transpiration_deficit = mean(transpiration_deficit, na.rm = TRUE),
    count    = n(),
    .groups  = 'drop'
  ) %>%
  group_by(mask_species, PSI_species) %>%
  mutate(
    total_pixels     = sum(count),
    pixel_percentage = count / total_pixels,
    bin_median       = sapply(as.character(PSI_bin), function(lbl) {
      # double-escaped backslashes so gsub() sees \[ \] \( \)
      nums <- as.numeric(strsplit(
        gsub("\\[|\\]|\\(|\\)", "", lbl),
        ","
      )[[1]])
      mean(nums)
    })
  ) %>%
  filter(pixel_percentage >= 0.001) %>%
  ungroup()

cb_palette <- c(Oak="#E69F00", Beech="#0072B2", Spruce="#009E73", Pine="#F0E442")
plot_df <- bin_summary %>%
  mutate(
    mask_species = factor(mask_species, levels = species),
    PSI_species  = factor(PSI_species,  levels = species),
    pixel_pct    = pixel_percentage * 100
  )

plot_df_clean <- plot_df %>%
  filter(!is.na(avg_transpiration_deficit), !is.na(bin_median))

fit_poly_curves <- function(df, key) {
  # Clean data: remove NAs and ensure we have enough points for Poly 3 (needs >= 4)
  df_clean <- df[is.finite(df$bin_median) & is.finite(df$avg_transpiration_deficit), ]
  vals <- df_clean$bin_median
  
  if (nrow(df_clean) < 4) {
    return(tibble(bin_median = numeric(0), fit = numeric(0), model = character(0)))
  }
  
  # Group info for logging
  m_sp <- key$mask_species
  p_sp <- key$PSI_species
  cat("\nFitting Poly 1-3 for Panel:", as.character(m_sp), "| Line:", as.character(p_sp), "\n")
  
  # Fit the three models
  m1 <- lm(avg_transpiration_deficit ~ poly(bin_median, 1, raw = TRUE), data = df_clean)
  m2 <- lm(avg_transpiration_deficit ~ poly(bin_median, 2, raw = TRUE), data = df_clean)
  m3 <- lm(avg_transpiration_deficit ~ poly(bin_median, 3, raw = TRUE), data = df_clean)
  
  # Calculate AICs
  aics <- c("poly1" = AIC(m1), "poly2" = AIC(m2), "poly3" = AIC(m3))
  best_model_name <- names(which.min(aics))
  best_model <- list("poly1" = m1, "poly2" = m2, "poly3" = m3)[[best_model_name]]
  
  cat("AICs -> P1:", round(aics[1],1), "P2:", round(aics[2],1), "P3:", round(aics[3],1), 
      "| Chosen:", best_model_name, "\n")
  
  # Create smooth prediction grid
  grid <- tibble(bin_median = seq(min(vals), max(vals), length.out = 100))
  grid$fit   <- predict(best_model, newdata = grid)
  grid$model <- best_model_name
  
  return(grid)
}

# 1. Compute fitted curves
fitted_df <- plot_df_clean %>%
  group_by(mask_species, PSI_species) %>%
  group_modify(~ fit_poly_curves(.x, .y)) %>%
  ungroup()

# 2. Build the Plot
p_tdiff_psi <- ggplot(plot_df_clean, aes(
  x     = bin_median, 
  y     = avg_transpiration_deficit, 
  color = PSI_species,
  size  = pixel_pct
)) +
  # Raw data points
  geom_point(alpha = 0.5) +
  # Best-fit polynomial lines
  geom_line(
    data = fitted_df,
    aes(
      x        = bin_median,
      y        = fit,
      color    = PSI_species,
      linetype = model, # Line type shows which poly degree was chosen
      group    = interaction(mask_species, PSI_species)
    ),
    linewidth = 1.1,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ mask_species, ncol = 2) +
  scale_color_manual(values = cb_palette) + 
  scale_size_continuous(range = c(0.5, 4)) +
  labs(
    x = "Soil Water Potential (bin median)",
    y = "Avg Transpiration Deficit",
    title = "Transpiration Deficit vs PSI: Polynomial AIC Selection"
  ) +
  theme_minimal() +
  theme(
    axis.text.x       = element_text(angle = 0, hjust = 0.5),
    plot.background   = element_rect(fill = "white", color = "white"),
    panel.background  = element_rect(fill = "white"),
    legend.background = element_rect(fill = "white", color = "white"),
    plot.title        = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title        = element_text(face = "bold", size = 14),
    axis.text         = element_text(color = "black", size = 12),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.5),
    strip.background  = element_rect(fill = "white", color = "black", linewidth = 0.5),
    strip.text        = element_text(face = "bold", size = 12),
    legend.position   = "right",
    legend.text       = element_text(size = 12)
  )

print(p_tdiff_psi)
