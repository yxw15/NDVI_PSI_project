# Set working directory and load packages
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")
library(dplyr)
library(ggplot2)
library(terra)
library(GGally)
library(patchwork)

# Create directory for figures
figures_dir <- "results_rootzone/Figures"
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# Load and combine rasters
PSI_Oak_2024_07_28_TIFF <- rast("results_rootzone/Project_mean/PSI_Oak_2024_07_28.tif")
PSI_Beech_2024_07_28_TIFF <- rast("results_rootzone/Project_mean/PSI_Beech_2024_07_28.tif")
PSI_Spruce_2024_07_28_TIFF <- rast("results_rootzone/Project_mean/PSI_Spruce_2024_07_28.tif")
PSI_Pine_2024_07_28_TIFF <- rast("results_rootzone/Project_mean/PSI_Pine_2024_07_28.tif")

names(PSI_Oak_2024_07_28_TIFF) <- "PSI_Oak_2024_07_28"
names(PSI_Beech_2024_07_28_TIFF) <- "PSI_Beech_2024_07_28"
names(PSI_Spruce_2024_07_28_TIFF) <- "PSI_Spruce_2024_07_28"
names(PSI_Pine_2024_07_28_TIFF) <- "PSI_Pine_2024_07_28"

PSI_2024_07_28_TIFF <- c(PSI_Oak_2024_07_28_TIFF,
                         PSI_Beech_2024_07_28_TIFF,
                         PSI_Spruce_2024_07_28_TIFF,
                         PSI_Pine_2024_07_28_TIFF)

# Randomly sample 1000 points (excluding NA values)
set.seed(123) # for reproducibility
sample_points <- spatSample(PSI_2024_07_28_TIFF, size = 1000, method = "random", 
                            na.rm = TRUE, as.df = TRUE)

# Check the structure of sampled data
cat("Sample data structure:\n")
print(str(sample_points))
cat("\nNumber of non-NA samples:", nrow(sample_points), "\n")

# Rename columns for better readability
colnames(sample_points) <- c("Oak", "Beech", "Spruce", "Pine")

# Option 1: Comprehensive scatter plot matrix using GGally
cat("Creating scatter plot matrix...\n")
scatter_matrix <- ggpairs(sample_points,
                          title = "Scatter Plot Matrix of PSI Values (1000 random points)",
                          lower = list(continuous = wrap("points", alpha = 0.6, size = 1)),
                          diag = list(continuous = wrap("barDiag", bins = 20)),
                          upper = list(continuous = wrap("cor", size = 4))) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Save scatter plot matrix
ggsave(file.path(figures_dir, "PSI_scatter_matrix.png"), scatter_matrix, 
       width = 12, height = 10, dpi = 300)
cat("Saved: PSI_scatter_matrix.png\n")

# Option 2: Individual scatter plots for each pair
cat("Creating individual scatter plots...\n")
layer_names <- c("Oak", "Beech", "Spruce", "Pine")
pairs <- combn(layer_names, 2, simplify = FALSE)

# Create individual scatter plots
plots <- list()

for (pair in pairs) {
  p <- ggplot(sample_points, aes(x = .data[[pair[1]]], y = .data[[pair[2]]])) +
    geom_point(alpha = 0.6, size = 2, color = "steelblue") +
    geom_smooth(method = "lm", se = TRUE, color = "red", fill = "pink", alpha = 0.3) +
    labs(x = paste0(pair[1], " PSI (kPa)"),
         y = paste0(pair[2], " PSI (kPa)"),
         title = paste0(pair[1], " vs ", pair[2])) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          panel.grid.major = element_line(color = "grey90"),
          panel.grid.minor = element_blank())
  
  plots[[paste0(pair[1], "_vs_", pair[2])]] <- p
  
  # Save each individual plot
  ggsave(file.path(figures_dir, paste0("PSI_", pair[1], "_vs_", pair[2], ".png")), p, 
         width = 8, height = 6, dpi = 300)
  cat(paste0("Saved: PSI_", pair[1], "_vs_", pair[2], ".png\n"))
}

# Arrange individual plots in a grid
cat("Arranging plots in grid...\n")
grid_plot <- wrap_plots(plots, ncol = 2) + 
  plot_annotation(title = "Pairwise Scatter Plots of PSI Values",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16)))

# Save grid of individual plots
ggsave(file.path(figures_dir, "PSI_pairwise_scatters_grid.png"), grid_plot, 
       width = 14, height = 12, dpi = 300)
cat("Saved: PSI_pairwise_scatters_grid.png\n")

# Option 3: Simple base R scatter matrix (saved as PNG)
cat("Creating base R scatter matrix...\n")
png(file.path(figures_dir, "PSI_baseR_scatter_matrix.png"), 
    width = 10, height = 10, units = "in", res = 300)
pairs(sample_points, 
      main = "Scatter Plot Matrix of PSI Values (1000 random points)",
      pch = 16, 
      cex = 0.8,
      col = rgb(0.2, 0.4, 0.6, 0.6),
      gap = 0.5)
dev.off()
cat("Saved: PSI_baseR_scatter_matrix.png\n")

# Additional analysis: Correlation matrix
cat("\nCorrelation Analysis:\n")
cor_matrix <- cor(sample_points, use = "complete.obs")
print("Correlation Matrix:")
print(round(cor_matrix, 3))

# Save correlation matrix as CSV
write.csv(round(cor_matrix, 3), file.path(figures_dir, "PSI_correlation_matrix.csv"))
cat("Saved: PSI_correlation_matrix.csv\n")

# Summary statistics
cat("\nSummary Statistics:\n")
summary_stats <- summary(sample_points)
print(summary_stats)

# Save summary statistics as CSV
summary_df <- as.data.frame(do.call(cbind, lapply(sample_points, summary)))
write.csv(summary_df, file.path(figures_dir, "PSI_summary_statistics.csv"))
cat("Saved: PSI_summary_statistics.csv\n")

# Create a comprehensive report of all files saved
saved_files <- list.files(figures_dir, pattern = "PSI_")
cat("\n=== FILES SAVED IN results_rootzone/Figures/ ===\n")
for (file in saved_files) {
  cat(paste0("- ", file, "\n"))
}

cat(paste0("\nAnalysis complete! All files saved in: ", figures_dir, "\n"))
cat(paste0("Total files created: ", length(saved_files), "\n"))