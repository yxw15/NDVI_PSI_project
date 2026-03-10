setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")
load("results/Data/AllSpecies_AllMonths_rootzone.RData")

data_all <- combined %>% filter(Quantiles > 0) %>% filter(month %in% c("May", "June", "July", "August"))
data_all <- na.omit(data_all)
data <- data_all %>% filter(month %in% c("July", "August")) %>% filter(year < 2023)

# ----- theme -----
base_theme <- theme_minimal() +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.text = element_text(color = "black", size = 14),
    legend.position = "bottom",
    plot.title  = element_text(hjust = 0.5, size = 18, color = "black"),
    axis.title  = element_text(size = 16),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
    axis.text.y = element_text(angle = 0, hjust = 0.5, size = 12),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey92", linewidth = 0.25),
    panel.border = element_blank(),
    strip.text = element_text(size = 14)
  )

# ----- NDVI PSI month linear -----
plot_NDVI_PSI_month_species_linear_orig <- function(data, figure_output = NULL) {
  
  library(ggplot2)
  library(dplyr)
  
  # ---- base theme ----
  base_theme <- theme_minimal() +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      legend.text = element_text(color = "black", size = 14),
      legend.position = "bottom",
      plot.title  = element_text(hjust = 0.5, size = 18, color = "black"),
      axis.title  = element_text(size = 16),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
      axis.text.y = element_text(angle = 0, hjust = 0.5, size = 12),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
      panel.grid.minor = element_line(color = "grey92", linewidth = 0.25),
      panel.border = element_blank(),
      strip.text = element_text(size = 14)
    )
  
  # ---- data prep ----
  data <- data %>% filter(Quantiles > 0)
  
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  month_order   <- c("May", "June", "July", "August")
  
  data$species <- factor(data$species, levels = species_order)
  data$month   <- factor(data$month, levels = month_order)
  
  # ---- blue gradient: light → dark ----
  blue_palette <- c(
    "May"    = "#c6dbef",  # very light blue
    "June"   = "#6baed6",
    "July"   = "#3182bd",
    "August" = "#08519c"   # dark blue
  )
  
  # ---- plot ----
  p1 <- ggplot(data, aes(x = soil_water_potential, y = Quantiles, color = month)) +
    geom_smooth(method = "lm", se = FALSE, aes(group = month), linewidth = 1.1) +
    facet_wrap(~ species, scales = "free_x") +
    scale_color_manual(values = blue_palette) +
    labs(
      x = "soil water potential (kPa)",
      y = "NDVI quantiles (rank)",
      color = ""
    ) +
    ylim(0, 22) +
    base_theme
  
  if (!is.null(figure_output)) {
    ggsave(filename = figure_output, plot = p1, width = 12, height = 8)
  }
  
  return(p1)
}

# usage
plot_NDVI_PSI_month_species_linear_orig(
  data_all,
  "results_rootzone/Figures_till2022/SI_PSI/SI_NDVI_PSI_linear_month.png"
)


# ----- NDVI TDiff month linear -----
plot_NDVI_TDiff_month_species_linear_orig <- function(data, figure_output = NULL) {
  
  library(ggplot2)
  library(dplyr)
  
  # ---- base theme ----
  base_theme <- theme_minimal() +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      legend.text = element_text(color = "black", size = 14),
      legend.position = "bottom",
      plot.title  = element_text(hjust = 0.5, size = 18, color = "black"),
      axis.title  = element_text(size = 16),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
      axis.text.y = element_text(angle = 0, hjust = 0.5, size = 12),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
      panel.grid.minor = element_line(color = "grey92", linewidth = 0.25),
      panel.border = element_blank(),
      strip.text = element_text(size = 14)
    )
  
  # ---- data prep ----
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  month_order   <- c("May", "June", "July", "August")
  
  data$species <- factor(data$species, levels = species_order)
  data$month   <- factor(data$month, levels = month_order)
  
  # ---- blue gradient: light → dark ----
  blue_palette <- c(
    "May"    = "#c6dbef",
    "June"   = "#6baed6",
    "July"   = "#3182bd",
    "August" = "#08519c"
  )
  
  # ---- plot ----
  p1 <- ggplot(data, aes(x = transpiration_deficit, y = Quantiles, color = month)) +
    geom_smooth(method = "lm", se = FALSE, aes(group = month), linewidth = 1.1) +
    facet_wrap(~ species, scales = "free_x") +
    scale_color_manual(values = blue_palette) +
    labs(
      x = "transpiration deficit (mm)",
      y = "NDVI quantiles (rank)",
      color = ""
    ) +
    base_theme
  
  if (!is.null(figure_output)) {
    ggsave(filename = figure_output, plot = p1, width = 10, height = 8)
  }
  
  return(p1)
}

# usage
plot_NDVI_TDiff_month_species_linear_orig(
  data_all,
  "results_rootzone/Figures_till2022/SI_PSI/SI_NDVI_TDiff_linear_month.png"
)


# ----- NDVI PSI density linear -----
plot_density_NDVI_PSI_linear <- function(data, swp_bin_width = 50, output_path) {
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  # ---- base theme (consistent across figures) ----
  base_theme <- theme_minimal() +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      legend.text = element_text(color = "black", size = 14),
      legend.position = "bottom",
      plot.title  = element_text(hjust = 0.5, size = 18, color = "black"),
      axis.title  = element_text(size = 16),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
      axis.text.y = element_text(angle = 0, hjust = 0.5, size = 12),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
      panel.grid.minor = element_line(color = "grey92", linewidth = 0.25),
      panel.border = element_blank(),
      strip.text = element_text(size = 14)
    )
  
  # ---- species order ----
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  
  species_df <- data %>% 
    filter(species %in% species_order) %>%
    mutate(species = factor(species, levels = species_order))
  
  # ---- binning ----
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
  
  # ---- density calculation ----
  density_df <- species_binned %>%
    group_by(species, quantile_bin, swp_bin) %>%
    summarize(count = n(), .groups = "drop") %>%
    group_by(species, swp_bin) %>%
    mutate(density_percent = (count / sum(count)) * 100) %>%
    ungroup()
  
  species_density <- species_binned %>%
    left_join(density_df, by = c("species", "quantile_bin", "swp_bin"))
  
  max_density <- ceiling(max(species_density$density_percent, na.rm = TRUE))
  density_breaks <- seq(0, max_density, length.out = 5)
  
  # ---- plot ----
  p <- ggplot(species_density, 
              aes(x = soil_water_potential, 
                  y = Quantiles, 
                  color = density_percent)) +
    
    geom_point(size = 2.5, alpha = 0.9) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1) +
    
    # ---- light → dark blue ----
  scale_color_gradientn(
    colors = c("#deebf7", "#9ecae1", "#3182bd", "#08519c"),
    name = "density (%)",
    trans = "sqrt",
    breaks = density_breaks
  ) +
    
    labs(
      title = "",
      x = "soil water potential (kPa)",
      y = "NDVI quantiles (rank)"
    ) +
    
    facet_wrap(~species, ncol = 2) +
    
    base_theme
  
  print(p)
  ggsave(output_path, plot = p, width = 12, height = 8, dpi = 300)
}

plot_density_NDVI_PSI_linear(data,
                             swp_bin_width = 50,
                             output_path = "results_rootzone/Figures_till2022/SI_PSI/SI_NDVI_PSI_density.png")


# ----- NDVI TDiff density linear -----
plot_density_NDVI_TDiff_linear <- function(data, tdiff_bin_width = 3, output_path) {
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  # ---- base theme (consistent across figures) ----
  base_theme <- theme_minimal() +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      legend.text = element_text(color = "black", size = 14),
      legend.position = "bottom",
      plot.title  = element_text(hjust = 0.5, size = 18, color = "black"),
      axis.title  = element_text(size = 16),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
      axis.text.y = element_text(angle = 0, hjust = 0.5, size = 12),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
      panel.grid.minor = element_line(color = "grey92", linewidth = 0.25),
      panel.border = element_blank(),
      strip.text = element_text(size = 14)
    )
  
  # ---- blue gradient: light → dark ----
  blue_grad <- c("#deebf7", "#9ecae1", "#3182bd", "#08519c")
  
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  
  # 1) Clean + restrict to target species
  df <- data %>%
    mutate(
      species = as.character(species),
      Quantiles = suppressWarnings(as.numeric(Quantiles)),
      transpiration_deficit = suppressWarnings(as.numeric(transpiration_deficit))
    ) %>%
    filter(species %in% species_order) %>%
    mutate(species = factor(species, levels = species_order)) %>%
    drop_na(Quantiles, transpiration_deficit)
  
  if (nrow(df) == 0) stop("No rows left after filtering/NA removal.")
  
  # 2) Robust bin edges for transpiration_deficit
  tmin <- min(df$transpiration_deficit, na.rm = TRUE)
  tmax <- max(df$transpiration_deficit, na.rm = TRUE)
  if (!is.finite(tmin) || !is.finite(tmax)) stop("Non-finite transpiration_deficit range.")
  if (tmin == tmax) { tmin <- tmin - tdiff_bin_width/2; tmax <- tmax + tdiff_bin_width/2 }
  
  edges <- seq(
    floor(tmin / tdiff_bin_width) * tdiff_bin_width,
    ceiling(tmax / tdiff_bin_width) * tdiff_bin_width,
    by = tdiff_bin_width
  )
  if (length(edges) < 2) edges <- c(tmin, tmax)
  
  species_binned <- df %>%
    mutate(
      quantile_bin = cut(Quantiles, breaks = seq(0, 22, 1), include.lowest = TRUE, right = FALSE),
      tdiff_bin    = cut(transpiration_deficit, breaks = edges, include.lowest = TRUE, right = FALSE)
    ) %>%
    drop_na(quantile_bin, tdiff_bin)
  
  if (nrow(species_binned) == 0) stop("All rows dropped during binning; check Quantiles range and tdiff_bin_width.")
  
  # 3) Density within each (species, tdiff_bin)
  density_df <- species_binned %>%
    count(species, quantile_bin, tdiff_bin, name = "count") %>%
    group_by(species, tdiff_bin) %>%
    mutate(density_percent = count / sum(count) * 100) %>%
    ungroup()
  
  species_density <- species_binned %>%
    left_join(density_df, by = c("species", "quantile_bin", "tdiff_bin"))
  
  # 4) Legend safety
  max_density <- suppressWarnings(max(species_density$density_percent, na.rm = TRUE))
  has_color   <- is.finite(max_density) && max_density > 0
  if (!has_color) max_density <- 1
  density_breaks <- unique(pretty(c(0, max_density), n = 5))
  
  # 5) Only draw smoother where there’s enough variation
  smooth_df <- species_density %>%
    group_by(species) %>%
    filter(n() >= 3,
           dplyr::n_distinct(transpiration_deficit) > 1,
           dplyr::n_distinct(Quantiles) > 1) %>%
    ungroup()
  
  p <- ggplot(
    species_density,
    aes(x = transpiration_deficit, y = Quantiles, color = density_percent)
  ) +
    geom_point(size = 2.5, alpha = 0.9, na.rm = TRUE) +
    { if (nrow(smooth_df) > 0)
      geom_smooth(
        data = smooth_df,
        aes(x = transpiration_deficit, y = Quantiles),
        method = "lm", formula = y ~ x, se = TRUE,
        color = "black", linewidth = 1, na.rm = TRUE
      ) else NULL
    } +
    { if (has_color)
      scale_color_gradientn(
        colors = blue_grad,
        name = "density (%)",
        trans = "sqrt",
        breaks = density_breaks,
        limits = c(0, max_density),
        na.value = NA
      ) else
        scale_color_continuous(guide = "none", na.value = NA)
    } +
    labs(
      title = "",
      x = "transpiration deficit (mm)",
      y = "NDVI quantiles (rank)"
    ) +
    facet_wrap(~species, ncol = 2, drop = TRUE) +
    base_theme
  
  print(p)
  ggsave(output_path, plot = p, width = 12, height = 8, dpi = 300)
}

# usage
plot_density_NDVI_TDiff_linear(
  data,
  tdiff_bin_width = 3,
  output_path = "results_rootzone/Figures_till2022/SI_PSI/NDVI_Tdiff_density_till2022.png"
)



# ----- NDVI PSI linear denstiy -----'
# =========================================================
# 8-panel figure (4 rows = species; 2 cols = linear | density)
# Panels tagged (a) ... (h) at top-left of each subpanel
# =========================================================

plot_NDVI_PSI_linear_density_8panel <- function(data_linear,
                                                data_density = data_linear,
                                                swp_bin_width = 50,
                                                figure_output = NULL,
                                                width = 12,
                                                height = 16,
                                                dpi = 300) {
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(grid)
  
  # -------------------------
  # base theme
  # -------------------------
  axis_title_size <- 16
  
  base_theme <- theme_minimal() +
    theme(
      plot.background   = element_rect(fill = "white", color = NA),
      panel.background  = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      legend.position   = "bottom",
      legend.text       = element_text(color = "black", size = 12),
      legend.title      = element_text(color = "black", size = 13),
      plot.title        = element_text(hjust = 0.5, size = 18, color = "black"),
      plot.tag          = element_text(face = "plain", size = axis_title_size),
      plot.tag.position = c(0.02, 0.98),
      axis.title        = element_text(size = axis_title_size),
      axis.text.x       = element_text(angle = 0, hjust = 0.5, size = 12),
      axis.text.y       = element_text(angle = 0, hjust = 0.5, size = 12),
      panel.grid.major  = element_line(color = "grey85", linewidth = 0.4),
      panel.grid.minor  = element_line(color = "grey92", linewidth = 0.25),
      panel.border      = element_blank(),
      strip.text        = element_text(size = 14)
    )
  
  # -------------------------
  # orders
  # -------------------------
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  month_order   <- c("May", "June", "July", "August")
  
  # -------------------------
  # month colors
  # -------------------------
  blue_palette <- c(
    "May"    = "#c6dbef",
    "June"   = "#6baed6",
    "July"   = "#3182bd",
    "August" = "#08519c"
  )
  
  # =========================
  # 1) Linear data prep
  # =========================
  lin_df <- data_linear %>%
    filter(Quantiles > 0, species %in% species_order, month %in% month_order) %>%
    mutate(
      species = factor(species, levels = species_order),
      month   = factor(month, levels = month_order)
    )
  
  make_linear_one_species <- function(sp, tag_label, show_x = FALSE) {
    ggplot(
      filter(lin_df, species == sp),
      aes(x = soil_water_potential, y = Quantiles, color = month)
    ) +
      geom_smooth(
        method = "lm", se = FALSE,
        aes(group = month),
        linewidth = 1.1
      ) +
      scale_color_manual(values = blue_palette, name = "Month") +
      labs(
        tag   = tag_label,
        title = NULL,
        x     = if (show_x) "Soil water potential (kPa)" else NULL,
        y     = "NDVI quantiles (rank)"
      ) +
      ylim(0, 22) +
      base_theme +
      theme(
        legend.position = "bottom",
        axis.title.x = if (show_x) {
          element_text(size = axis_title_size)
        } else {
          element_blank()
        }
      )
  }
  
  # =========================
  # 2) Density data prep
  # =========================
  den_df <- data_density %>%
    filter(species %in% species_order, Quantiles > 0) %>%
    mutate(species = factor(species, levels = species_order))
  
  swp_breaks <- seq(
    floor(min(den_df$soil_water_potential, na.rm = TRUE) / swp_bin_width) * swp_bin_width,
    ceiling(max(den_df$soil_water_potential, na.rm = TRUE) / swp_bin_width) * swp_bin_width,
    by = swp_bin_width
  )
  
  den_binned <- den_df %>%
    mutate(
      quantile_bin = cut(
        Quantiles,
        breaks = seq(0, 22, 1),
        include.lowest = TRUE,
        right = FALSE
      ),
      swp_bin = cut(
        soil_water_potential,
        breaks = swp_breaks,
        include.lowest = TRUE,
        right = FALSE
      )
    ) %>%
    drop_na(quantile_bin, swp_bin)
  
  density_df <- den_binned %>%
    group_by(species, quantile_bin, swp_bin) %>%
    summarize(count = n(), .groups = "drop") %>%
    group_by(species, swp_bin) %>%
    mutate(density_percent = (count / sum(count)) * 100) %>%
    ungroup()
  
  species_density <- den_binned %>%
    left_join(density_df, by = c("species", "quantile_bin", "swp_bin"))
  
  max_density <- ceiling(max(species_density$density_percent, na.rm = TRUE))
  density_breaks <- pretty(c(0, max_density), n = 5)
  
  make_density_one_species <- function(sp, tag_label, show_x = FALSE) {
    ggplot(
      filter(species_density, species == sp),
      aes(x = soil_water_potential, y = Quantiles, color = density_percent)
    ) +
      geom_point(size = 2.5, alpha = 0.9) +
      geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1) +
      scale_color_gradientn(
        colors = c("#deebf7", "#9ecae1", "#3182bd", "#08519c"),
        name   = "Density (%)",
        trans  = "sqrt",
        breaks = density_breaks,
        guide  = guide_colorbar(
          title.position = "top",
          title.hjust = 0.5,
          barwidth = unit(14, "lines"),   # wider color bar
          barheight = unit(0.8, "lines")
        )
      ) +
      labs(
        tag   = tag_label,
        title = NULL,
        x     = if (show_x) "Soil water potential (kPa)" else NULL,
        y     = "NDVI quantiles (rank)"
      ) +
      ylim(0, 22) +
      base_theme +
      theme(
        legend.position = "bottom",
        axis.title.x = if (show_x) {
          element_text(size = axis_title_size)
        } else {
          element_blank()
        }
      )
  }
  
  # =========================
  # 3) Build plots
  # =========================
  p_lin_oak    <- make_linear_one_species("Oak",    "(a)", show_x = FALSE)
  p_den_oak    <- make_density_one_species("Oak",   "(b)", show_x = FALSE)
  p_lin_beech  <- make_linear_one_species("Beech",  "(c)", show_x = FALSE)
  p_den_beech  <- make_density_one_species("Beech", "(d)", show_x = FALSE)
  p_lin_spruce <- make_linear_one_species("Spruce", "(e)", show_x = FALSE)
  p_den_spruce <- make_density_one_species("Spruce","(f)", show_x = FALSE)
  p_lin_pine   <- make_linear_one_species("Pine",   "(g)", show_x = TRUE)
  p_den_pine   <- make_density_one_species("Pine",  "(h)", show_x = TRUE)
  
  # Left column: one shared legend for a,c,e,g
  left_column <- (p_lin_oak / p_lin_beech / p_lin_spruce / p_lin_pine) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  # Right column: one shared legend for b,d,f,h
  right_column <- (p_den_oak / p_den_beech / p_den_spruce / p_den_pine) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  
  # Final 2-column layout
  big <- left_column | right_column
  
  if (!is.null(figure_output)) {
    ggsave(figure_output, plot = big, width = width, height = height, dpi = dpi)
  }
  
  return(big)
}

# =========================
# Example usage
# =========================
big_fig <- plot_NDVI_PSI_linear_density_8panel(
  data_linear  = data_all,   # your linear dataset
  data_density = data,       # your density dataset (or same as data_all)
  swp_bin_width = 50,
  figure_output = "results_rootzone/Figures_till2022/SI_PSI/SI_NDVI_PSI_linear_density.png",
  width = 12, height = 16, dpi = 300
)

# ----- NDVI TDiff linear denstiy -----'
# =========================================================
# 8-panel figure for NDVI ~ Tdiff
# 4 rows = Oak/Beech/Spruce/Pine
# 2 cols = monthly linear fits | density
# Panels tagged (a) ... (h) at top-left of each subpanel
# =========================================================

plot_NDVI_TDiff_linear_density_8panel <- function(data_linear,
                                                  data_density = data_linear,
                                                  tdiff_bin_width = 3,
                                                  figure_output = NULL,
                                                  width = 12,
                                                  height = 16,
                                                  dpi = 300) {
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  
  # ---- base theme (same as yours) ----
  base_theme <- theme_minimal() +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      legend.text = element_text(color = "black", size = 14),
      legend.position = "bottom",
      plot.title  = element_text(hjust = 0.5, size = 18, color = "black"),
      axis.title  = element_text(size = 16),
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
      axis.text.y = element_text(angle = 0, hjust = 0.5, size = 12),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
      panel.grid.minor = element_line(color = "grey92", linewidth = 0.25),
      panel.border = element_blank(),
      strip.text = element_text(size = 14)
    )
  
  # ---- orders ----
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  month_order   <- c("May", "June", "July", "August")
  
  # ---- month colors (light -> dark blue) ----
  blue_palette <- c(
    "May"    = "#c6dbef",
    "June"   = "#6baed6",
    "July"   = "#3182bd",
    "August" = "#08519c"
  )
  
  # =========================
  # 1) Linear data prep
  # =========================
  lin_df <- data_linear %>%
    filter(species %in% species_order, month %in% month_order) %>%
    mutate(
      species = factor(species, levels = species_order),
      month   = factor(month, levels = month_order)
    )
  
  make_linear_one_species <- function(sp) {
    ggplot(filter(lin_df, species == sp),
           aes(x = transpiration_deficit, y = Quantiles, color = month)) +
      geom_smooth(method = "lm", se = FALSE, aes(group = month), linewidth = 1.1) +
      scale_color_manual(values = blue_palette) +
      labs(
        title = NULL,
        x = "transpiration deficit (mm)",
        y = "NDVI quantiles (rank)",
        color = ""
      ) +
      base_theme
  }
  
  # =========================
  # 2) Density data prep
  # =========================
  df <- data_density %>%
    mutate(
      species = as.character(species),
      Quantiles = suppressWarnings(as.numeric(Quantiles)),
      transpiration_deficit = suppressWarnings(as.numeric(transpiration_deficit))
    ) %>%
    filter(species %in% species_order) %>%
    mutate(species = factor(species, levels = species_order)) %>%
    drop_na(Quantiles, transpiration_deficit)
  
  if (nrow(df) == 0) stop("No rows left after filtering/NA removal.")
  
  # robust bin edges
  tmin <- min(df$transpiration_deficit, na.rm = TRUE)
  tmax <- max(df$transpiration_deficit, na.rm = TRUE)
  if (!is.finite(tmin) || !is.finite(tmax)) stop("Non-finite transpiration_deficit range.")
  if (tmin == tmax) { tmin <- tmin - tdiff_bin_width/2; tmax <- tmax + tdiff_bin_width/2 }
  
  edges <- seq(
    floor(tmin / tdiff_bin_width) * tdiff_bin_width,
    ceiling(tmax / tdiff_bin_width) * tdiff_bin_width,
    by = tdiff_bin_width
  )
  if (length(edges) < 2) edges <- c(tmin, tmax)
  
  species_binned <- df %>%
    mutate(
      quantile_bin = cut(Quantiles, breaks = seq(0, 22, 1), include.lowest = TRUE, right = FALSE),
      tdiff_bin    = cut(transpiration_deficit, breaks = edges, include.lowest = TRUE, right = FALSE)
    ) %>%
    drop_na(quantile_bin, tdiff_bin)
  
  if (nrow(species_binned) == 0) stop("All rows dropped during binning; check Quantiles range and tdiff_bin_width.")
  
  density_df <- species_binned %>%
    count(species, quantile_bin, tdiff_bin, name = "count") %>%
    group_by(species, tdiff_bin) %>%
    mutate(density_percent = count / sum(count) * 100) %>%
    ungroup()
  
  species_density <- species_binned %>%
    left_join(density_df, by = c("species", "quantile_bin", "tdiff_bin"))
  
  # legend safety
  max_density <- suppressWarnings(max(species_density$density_percent, na.rm = TRUE))
  has_color   <- is.finite(max_density) && max_density > 0
  if (!has_color) max_density <- 1
  density_breaks <- unique(pretty(c(0, max_density), n = 5))
  
  # smoother only when valid
  smooth_df <- species_density %>%
    group_by(species) %>%
    filter(n() >= 3,
           dplyr::n_distinct(transpiration_deficit) > 1,
           dplyr::n_distinct(Quantiles) > 1) %>%
    ungroup()
  
  blue_grad <- c("#deebf7", "#9ecae1", "#3182bd", "#08519c")
  
  make_density_one_species <- function(sp) {
    p <- ggplot(filter(species_density, species == sp),
                aes(x = transpiration_deficit, y = Quantiles, color = density_percent)) +
      geom_point(size = 2.5, alpha = 0.9, na.rm = TRUE)
    
    # add smoother if possible for this species
    sp_smooth <- filter(smooth_df, species == sp)
    if (nrow(sp_smooth) > 0) {
      p <- p + geom_smooth(
        data = sp_smooth,
        aes(x = transpiration_deficit, y = Quantiles),
        method = "lm", formula = y ~ x, se = TRUE,
        color = "black", linewidth = 1, na.rm = TRUE
      )
    }
    
    p +
      { if (has_color)
        scale_color_gradientn(
          colors = blue_grad,
          name = "density (%)",
          trans = "sqrt",
          breaks = density_breaks,
          limits = c(0, max_density),
          na.value = NA
        ) else
          scale_color_continuous(guide = "none", na.value = NA)
      } +
      labs(
        title = NULL,
        x = "transpiration deficit (mm)",
        y = "NDVI quantiles (rank)"
      ) +
      base_theme
  }
  
  # =========================
  # 3) Build 8 plots + layout
  # =========================
  p_lin_oak    <- make_linear_one_species("Oak")
  p_den_oak    <- make_density_one_species("Oak")
  p_lin_beech  <- make_linear_one_species("Beech")
  p_den_beech  <- make_density_one_species("Beech")
  p_lin_spruce <- make_linear_one_species("Spruce")
  p_den_spruce <- make_density_one_species("Spruce")
  p_lin_pine   <- make_linear_one_species("Pine")
  p_den_pine   <- make_density_one_species("Pine")
  
  big <- (p_lin_oak    | p_den_oak) /
    (p_lin_beech  | p_den_beech) /
    (p_lin_spruce | p_den_spruce) /
    (p_lin_pine   | p_den_pine)
  
  # Tag each subpanel (a)-(h) at top-left
  big <- big +
    plot_annotation(tag_levels = "a") &
    theme(
      plot.tag = element_text(face = "bold", size = 14),
      plot.tag.position = c(0.02, 0.98)
    )
  
  if (!is.null(figure_output)) {
    ggsave(figure_output, plot = big, width = width, height = height, dpi = dpi)
  }
  
  return(big)
}

# =========================
# Example usage
# =========================
big_fig_tdiff <- plot_NDVI_TDiff_linear_density_8panel(
  data_linear  = data_all,
  data_density = data,
  tdiff_bin_width = 3,
  figure_output = "results_rootzone/Figures_till2022/SI_PSI/SI_NDVI_TDiff_linear_density.png",
  width = 12, height = 16, dpi = 300
)
