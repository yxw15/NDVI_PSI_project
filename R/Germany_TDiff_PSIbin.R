setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

library(terra)
library(dplyr)

# -----------------------------
# USER SETTINGS
# -----------------------------
months   <- c("July", "August")
species  <- c("Oak", "Beech", "Spruce", "Pine")
years    <- 2003:2022

base_dir <- "results_monthly_rootzone"

# -----------------------------
# FUNCTION TO READ & MERGE
# -----------------------------
read_month_species_year <- function(month, species, year) {
  
  psi_path   <- file.path(base_dir, month, species, paste0("psi_", year, ".tif"))
  tdiff_path <- file.path(base_dir, month, species, paste0("tdiff_", year, ".tif"))
  
  # Skip if files missing
  if (!file.exists(psi_path) | !file.exists(tdiff_path)) {
    message("Missing files: ", month, " ", species, " ", year)
    return(NULL)
  }
  
  psi   <- rast(psi_path)
  tdiff <- rast(tdiff_path)
  
  psi_df   <- as.data.frame(psi, xy = TRUE)
  tdiff_df <- as.data.frame(tdiff, xy = TRUE)
  
  names(psi_df)[3]   <- "soil_water_potential"
  names(tdiff_df)[3] <- "transpiration_deficit"
  
  df <- left_join(psi_df, tdiff_df, by = c("x", "y")) %>%
    mutate(
      month   = month,
      year    = year,
      species = species
    )
  
  return(df)
}

# -----------------------------
# LOOP OVER EVERYTHING
# -----------------------------
final_df <- bind_rows(
  lapply(months, function(m)
    lapply(species, function(s)
      lapply(years, function(y)
        read_month_species_year(m, s, y)
      )
    )
  ) |> unlist(recursive = FALSE)
)

# -----------------------------
# OPTIONAL CLEANING
# -----------------------------
final_df <- final_df %>%
  filter(
    !is.na(soil_water_potential),
    !is.na(transpiration_deficit)
  )

# -----------------------------
# CHECK RESULT
# -----------------------------
str(final_df)
head(final_df)

# -----------------------------
# SAVE RESULT
# -----------------------------
save(final_df, file = "results_rootzone/Data/Germany_PSI_TDiff_all_species_months_years.RData")

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

TDiff_PSIbin_df <- TDiff_PSIbin(final_df)

plot_TDiff_PSIbin_onlyA <- function(data, figure_output) {
  
  # -------------------------
  # Libraries & preprocessing
  # -------------------------
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(tibble)
  library(purrr)
  library(scales)
  
  # Prepare data
  df <- TDiff_PSIbin(data)
  df <- na.omit(df)
  
  # Species & palette
  cb_palette <- c("Oak"   = "#E69F00",
                  "Beech" = "#0072B2",
                  "Spruce"= "#009E73",
                  "Pine"  = "#F0E442")
  species_levels <- c("Oak", "Beech", "Spruce", "Pine")
  df$species <- factor(df$species, levels = species_levels)
  
  xcol <- "bin_median"
  ycol <- "avg_transpiration_deficit"
  
  # -------------------------
  # Fit models & choose best degree by AIC
  # -------------------------
  degrees <- c(1, 2, 3)
  fit_one_species <- function(sp) {
    sdf <- df %>% filter(species == sp) %>% select(all_of(c(xcol, ycol)))
    if (nrow(sdf) < 5) {
      return(list(species = sp, best_deg = NA, best_fit = NULL))
    }
    fits <- lapply(degrees, function(d) {
      fm <- reformulate(termlabels = paste0("poly(", xcol, ", ", d, ")"), response = ycol)
      tryCatch(lm(fm, data = sdf), error = function(e) NULL)
    })
    aics <- sapply(fits, function(m) if (is.null(m)) Inf else AIC(m))
    best_idx <- which.min(aics)
    list(species = sp, best_deg = degrees[best_idx], best_fit = fits[[best_idx]])
  }
  fits_all <- lapply(species_levels, fit_one_species)
  
  # ---------------------------------
  # Predictions
  # ---------------------------------
  pred_list <- lapply(fits_all, function(ff) {
    if (is.null(ff$best_fit)) return(NULL)
    sp <- ff$species
    sp_df <- df %>% filter(species == sp)
    xr <- range(sp_df[[xcol]], na.rm = TRUE)
    xseq <- seq(xr[1], xr[2], length.out = 200)
    nd <- data.frame(bin_median = xseq)
    nd$species <- factor(sp, levels = species_levels)
    nd$predicted <- predict(ff$best_fit, newdata = nd)
    nd
  })
  pred_df <- bind_rows(pred_list)
  
  # -------------
  # Final plot: ONLY points + fitted lines
  # -------------
  plot_clean <- ggplot() +
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
    scale_color_manual(values = cb_palette, name = "") +
    scale_shape_manual(values = c("Oak" = 16, "Beech" = 17, "Spruce" = 15, "Pine" = 18),
                       guide = "none") +
    scale_size_continuous(name = "",
                          range = c(1, 8),
                          labels = percent_format(accuracy = 1)) +
    guides(
      color = guide_legend(order = 1),
      size  = guide_legend(order = 2)
    ) +
    labs(x = "soil water potential (kPa)",
         y = "transpiration deficit") +
    theme_minimal() +
    theme(
      plot.title = element_blank(),        # remove title "(a)"
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title  = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14, face = "bold"),
      legend.position = "bottom",
      panel.background = element_rect(fill = "white"),
      plot.background  = element_rect(fill = "white", color = "white"),
      panel.grid       = element_blank()
    )
  
  print(plot_clean)
  
  if (!missing(figure_output)) {
    ggsave(figure_output, plot = plot_clean, width = 8, height = 6, dpi = 300)
  }
  
  invisible(plot_clean)
}

plot_TDiff_PSIbin_onlyA(final_df, "results_rootzone/Figures_till2022/supplementary/Germany_TDiff_PSIbin_onlyA_till2022.png")

#################################
# different species mathces
#################################

species_masks <- list(
  Oak    = rast("species_map_MODIS/Oak.tif"),
  Beech  = rast("species_map_MODIS/Beech.tif"),
  Spruce = rast("species_map_MODIS/Spruce.tif"),
  Pine   = rast("species_map_MODIS/Pine.tif")
)

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

species <- c("Oak", "Beech","Spruce","Pine")
months  <- c("July","August")
years   <- 2003:2022
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

save(full_df, file = "results_rootzone/Data/Germany_PSI_TDiff_diff_match_all_species_months_years.RData")

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

save(plot_df_clean, file = "results_rootzone/Data/Germany_PSI_TDiff_diff_match_all_species_months_years_plot_df_clean.RData")

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
  scale_color_manual(
    values = cb_palette, 
    # Use expression() for mathematical symbols
    name = expression(species~of~Psi[~~soil])
  ) + 
  scale_size_continuous(
    # Simple text for the percentage legend
    name = "pixel percentage (%)",
    range = c(0.5, 4)
  ) +
  labs(
    x = expression(bold(soil~water~potential~(kPa))), # Optional: add math to axis too
    y = "transpiration deficit (mm)",
    title = ""
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
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
    legend.position = "bottom",
    legend.box = "horizontal",        # Keeps the color legend above the others
    legend.box.just = "center",
    legend.spacing.x = unit(0.5, 'cm'), # Adds horizontal space between size and linetype
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0)
  )


print(p_tdiff_psi)
ggsave("results_rootzone/Figures_till2022/supplementary/TDiff_PSIbin_mix_fitted_4panel.png",
       plot=p_tdiff_psi, width=10, height=8, dpi=300)

message("Script complete ✅ Each panel = TDiff masked by panel species, each curve = PSI species masked by panel species.")

# ==========================
# Script complete
# Each panel = TDiff masked by panel species
# Each curve = PSI species masked by panel species
# Polynomial curves fitted using AIC selection
# ==========================

