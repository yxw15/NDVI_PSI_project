# ------------------------- Setup & Packages -----------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(terra)
  library(GGally)
  library(readr)
  library(patchwork)
})

# ---------------------------- Paths & Config ----------------------------------
root_dir    <- "/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project"
setwd(root_dir)

date_tag    <- "2024_07_28"
r_dir       <- "results_rootzone/Project_mean"
figures_dir <- "results_rootzone/Figures"
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# Output file names
csv_path        <- file.path(figures_dir, "PSI_sample_common_pixels_1000.csv")
pairs_ggally    <- file.path(figures_dir, "PSI_pairs_sample1000.png")
heat_path       <- file.path(figures_dir, "PSI_correlation_heatmap_sample1000.png")
pairs_custom    <- file.path(figures_dir, "PSI_selected_pairs.png")   # <-- pairs, not "paris"

message("Figures will be saved to: ", normalizePath(figures_dir))

# --------------------------- Load rasters -------------------------------------
PSI_Oak    <- rast(file.path(r_dir, paste0("PSI_Oak_",    date_tag, ".tif")))
PSI_Beech  <- rast(file.path(r_dir, paste0("PSI_Beech_",  date_tag, ".tif")))
PSI_Spruce <- rast(file.path(r_dir, paste0("PSI_Spruce_", date_tag, ".tif")))
PSI_Pine   <- rast(file.path(r_dir, paste0("PSI_Pine_",   date_tag, ".tif")))

names(PSI_Oak)    <- paste0("PSI_Oak_",    date_tag)
names(PSI_Beech)  <- paste0("PSI_Beech_",  date_tag)
names(PSI_Spruce) <- paste0("PSI_Spruce_", date_tag)
names(PSI_Pine)   <- paste0("PSI_Pine_",   date_tag)

# Combine (assumes identical grid; align beforehand if needed)
PSI_stack <- c(PSI_Oak, PSI_Beech, PSI_Spruce, PSI_Pine)

# -------------------- Sample same pixels & extract values ---------------------
set.seed(1234)

# Mask of cells where ALL four layers are non-NA
complete_mask <- terra::app(PSI_stack, function(v) as.integer(all(!is.na(v))))

available <- as.integer(terra::global(complete_mask, "sum", na.rm = TRUE)[1, 1])
message("Pixels with data for ALL species: ", available)
if (available == 0L) stop("No common non-NA pixels across species.")

# Choose sample size (up to 1000 or all if fewer)
n_samp <- min(1000L, available)

# Sample locations from the mask (identical pixels for all species)
pts <- terra::spatSample(
  x         = complete_mask,
  size      = n_samp,
  method    = "random",   # use "regular" for systematic sampling if preferred
  as.points = TRUE,
  na.rm     = TRUE,
  values    = FALSE
)

# Extract PSI at those exact pixels
samp_df <- terra::extract(PSI_stack, pts, xy = TRUE) |>
  dplyr::select(-ID) |>
  dplyr::rename(
    Oak    = !!names(PSI_Oak),
    Beech  = !!names(PSI_Beech),
    Spruce = !!names(PSI_Spruce),
    Pine   = !!names(PSI_Pine)
  )

# Save sampled values (with coordinates) for traceability
readr::write_csv(samp_df, csv_path)
message("Saved sampled values: ", normalizePath(csv_path))

# ----------------------------- Correlations -----------------------------------
num_df <- dplyr::select(samp_df, Oak, Beech, Spruce, Pine)

corr_pearson  <- cor(num_df, use = "complete.obs", method = "pearson")
corr_spearman <- cor(num_df, use = "complete.obs", method = "spearman")

print(round(corr_pearson,  3))
print(round(corr_spearman, 3))

# ----------------------------- Visualizations ---------------------------------
# A) Pairs plot over the same pixels (GGally)
p_pairs <- GGally::ggpairs(
  num_df,
  progress = FALSE,
  upper = list(continuous = GGally::wrap("cor", size = 3)),
  lower = list(continuous = "points"),
  diag  = list(continuous = "densityDiag")
)

# Robust saver for ggpairs
save_pairs <- function(p, path, width = 10, height = 10, dpi = 300) {
  tryCatch({
    if ("GGally" %in% loadedNamespaces() && "ggsave" %in% getNamespaceExports("GGally")) {
      GGally::ggsave(filename = path, plot = p, width = width, height = height, dpi = dpi)
    } else {
      ggplot2::ggsave(filename = path, plot = p, width = width, height = height, dpi = dpi)
    }
    message("Saved pairs plot: ", normalizePath(path))
  }, error = function(e) {
    png(path, width = width, height = height, units = "in", res = dpi)
    print(p)
    dev.off()
    message("Saved pairs plot via PNG device: ", normalizePath(path))
  })
}

save_pairs(p_pairs, pairs_ggally)

# B) Correlation heatmap (Pearson)
corr_long <- as.data.frame(as.table(corr_pearson))
names(corr_long) <- c("Var1", "Var2", "corr")

p_heat <- ggplot(corr_long, aes(Var1, Var2, fill = corr)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", corr)), size = 4) +
  scale_fill_gradient2(limits = c(-1, 1), midpoint = 0) +
  coord_fixed() +
  labs(title = "",
       x = NULL, y = NULL, fill = "Pearson r") +
  theme_minimal(base_size = 12) +
  theme(
    plot.background   = element_rect(fill = "white", color = "white"),
    panel.background  = element_rect(fill = "white"),
    legend.background = element_rect(fill = "white", color = "white"),
    plot.title        = element_text(hjust = 0.5, size = 18, face = "bold"),
    plot.subtitle     = element_text(hjust = 0.5, size = 14),
    axis.title        = element_text(face = "bold", size = 14),
    axis.text.y       = element_text(color = "black", size = 12),
    axis.text.x       = element_blank(),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    legend.position   = "top",
    legend.text       = element_text(size = 12),
    strip.background  = element_rect(fill = "white", color = "black", linewidth = 0.5),
    strip.text        = element_text(face = "bold", size = 12)
  )

ggsave(heat_path, p_heat, width = 6, height = 5, dpi = 300)
message("Saved heatmap: ", normalizePath(heat_path))

# C) Selected pairwise scatter plots with LM line + r label
plot_pair <- function(df, xvar, yvar) {
  r <- suppressWarnings(cor(df[[xvar]], df[[yvar]], use = "complete.obs"))
  ggplot(df, aes(x = .data[[xvar]], y = .data[[yvar]])) +
    geom_point(alpha = 0.4, size = 0.8) +
    geom_smooth(method = "lm", se = FALSE) +
    annotate("text",
             x = min(df[[xvar]], na.rm = TRUE),
             y = max(df[[yvar]], na.rm = TRUE),
             label = paste0("r = ", round(r, 2)),
             hjust = 0, vjust = 1, size = 5, fontface = "bold") +
    labs(x = xvar, y = yvar) +
    theme_minimal(base_size = 12) +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
}

pairs_to_plot <- list(
  c("Oak", "Beech"),
  c("Oak", "Spruce"),
  c("Oak", "Pine"),
  c("Beech", "Spruce"),
  c("Beech", "Pine"),
  c("Spruce", "Pine")
)

plots <- lapply(pairs_to_plot, function(p) plot_pair(num_df, p[1], p[2]))
combined <- (plots[[1]] | plots[[2]] | plots[[3]]) /
  (plots[[4]] | plots[[5]] | plots[[6]])

print(combined)
pairs_custom <- file.path(figures_dir, "PSI_pairs.png")
ggsave(pairs_custom, combined, width = 12, height = 8, dpi = 300)
message("Saved custom pairs plot: ", normalizePath(pairs_custom))

# -------------------------------- Notes ---------------------------------------
# • “Same pixel” is guaranteed by sampling from a mask where all species are non-NA.
# • Switch to method = "regular" in spatSample() for a spatially even sample.
# • Increase/decrease n_samp by changing the 1000L cap.
# • If rasters are on different grids, align before this (project/crop/resample).

list.files(figures_dir, pattern = "PSI", full.names = TRUE)
