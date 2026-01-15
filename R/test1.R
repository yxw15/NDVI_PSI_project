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
      legend.position = "bottom",
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14, face = "bold"),
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
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

species <- c("Oak", "Beech", "Spruce", "Pine")
months  <- c("July", "August")
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

# --------------------------
# 5. Save dataset
# --------------------------
message("Saving full dataset as RData...")
save(full_df, file = "results_rootzone/Data/Germany_PSI_TDiff_mix_all_species_months_years.RData")
message("✅ Full dataset saved successfully!")

# --------------------------
# 6. Bin PSI & compute avg TDiff + pixel percentage
# --------------------------
message("Defining PSI binning function...")
TDiff_PSIbin <- function(df, bin_width = 50) {
  species_totals <- df %>% group_by(PSI_species) %>% summarise(total_pixels = n(), .groups="drop")
  
  breaks <- seq(floor(min(df$soil_water_potential, na.rm=TRUE)),
                ceiling(max(df$soil_water_potential, na.rm=TRUE)),
                by = bin_width)
  
  df %>%
    mutate(PSI_bin = cut(soil_water_potential, breaks,
                         include.lowest = TRUE, right=FALSE)) %>%
    group_by(PSI_species, PSI_bin) %>%
    summarise(
      avg_TDiff = mean(transpiration_deficit, na.rm = TRUE),
      count = n(),
      .groups="drop"
    ) %>%
    mutate(bin_median = sapply(as.character(PSI_bin), function(b) {
      mean(as.numeric(strsplit(gsub("\\[|\\]|\\(|\\)", "", b), ",")[[1]]))
    })) %>%
    left_join(species_totals, by="PSI_species") %>%
    mutate(pixel_pct = count / total_pixels)
}

# --------------------------
# 7. Fit polynomial curves per PSI species
# --------------------------
message("Defining curve fitting function...")
fit_curves <- function(df_bin) {
  degrees <- 1:3
  fits <- lapply(species, function(sp) {
    sdf <- df_bin %>% filter(PSI_species == sp)
    if(nrow(sdf) < 5) return(NULL)
    
    models <- lapply(degrees, function(d) lm(avg_TDiff ~ poly(bin_median, d), data=sdf))
    best <- models[[which.min(sapply(models, AIC))]]
    
    xs <- seq(min(sdf$bin_median), max(sdf$bin_median), length.out=200)
    data.frame(
      PSI_species = sp,
      bin_median  = xs,
      fit         = predict(best, newdata = data.frame(bin_median=xs)),
      model       = paste0("poly", which.min(sapply(models, AIC)))
    )
  })
  bind_rows(fits)
}

# --------------------------
# 8. Plot function (4-panel)
# --------------------------
message("Defining plotting function...")
plot_TDiff_PSIbin <- function(plot_df, cb_palette) {
  df_clean <- TDiff_PSIbin(plot_df)
  fitted_df <- fit_curves(df_clean)
  
  ggplot(df_clean, aes(x=bin_median, y=avg_TDiff, color=PSI_species, shape=PSI_species, size=pixel_pct)) +
    geom_point(alpha=0.85) +
    geom_line(data=fitted_df,
              aes(x=bin_median, y=fit, color=PSI_species, linetype=model,
                  group=interaction(PSI_species, model)),
              linewidth=1, inherit.aes = FALSE) +
    scale_color_manual(values=cb_palette, name="PSI species") +
    scale_shape_manual(values=c(Oak=16, Beech=17, Spruce=15, Pine=18), name="PSI species") +
    scale_size_continuous(name="pixel percentage (%)", range=c(1,8)) +
    facet_wrap(~mask_species, ncol=2) +
    labs(x="soil water potential (kPa)", y="transpiration deficit") +
    theme(
      axis.text.x = element_text(angle=0, hjust=0.5),
      plot.background = element_rect(fill="white", color="white"),
      panel.background=element_rect(fill="white"),
      legend.background=element_rect(fill="white", color="white"),
      plot.title = element_text(hjust=0.5, size=16, face="bold"),
      axis.title = element_text(face="bold", size=14),
      axis.text  = element_text(color="black", size=12),
      panel.border=element_rect(color="black", fill=NA, linewidth=0.5),
      strip.background=element_rect(fill="white", color="black", linewidth=0.5),
      strip.text=element_text(face="bold", size=12),
      legend.position="right",
      legend.text=element_text(size=12)
    )
}

# --------------------------
# 9. Generate combined plot
# --------------------------
message("Generating combined panel dataframe for plotting...")
plot_df_all <- bind_rows(lapply(species, function(sp) {
  build_panel_df(sp, months[1], years[1])  # Example: first month & year, can loop as needed
}))

cb_palette <- c(Oak="#E69F00", Beech="#0072B2", Spruce="#009E73", Pine="#F0E442")

message("Creating 4-panel TDiff x PSI plot...")
p_combined_TDiff <- plot_TDiff_PSIbin(plot_df_all, cb_palette)

print(p_combined_TDiff)

message("Saving figure...")
ggsave("results_rootzone/test/TDiff_PSIbin_mix_fitted_4panel.png",
       plot=p_combined_TDiff, width=10, height=8, dpi=300)

message("Script complete ✅ Each panel = TDiff masked by panel species, each curve = PSI species masked by panel species.")

# ==========================
# Script complete
# Each panel = TDiff masked by panel species
# Each curve = PSI species masked by panel species
# Polynomial curves fitted using AIC selection
# ==========================

