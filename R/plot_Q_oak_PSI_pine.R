library(terra)
library(ggplot2)
library(dplyr)
library(stringr)

# Set working directory and folders
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")
PSI.pine.folder <- "results/Pine"
NDVI.oak.folder <- "results/Oak/NDVI_mask"

# List files
psi_files <- list.files(PSI.pine.folder, pattern = "psi_\\d{4}\\.tif$", full.names = TRUE)
ndvi_mask_files <- list.files(NDVI.oak.folder, pattern = "NDVI_\\d{4}_mask\\.tif$", full.names = TRUE)

# Extract years
get_year <- function(filename) {
  str_extract(basename(filename), "\\d{4}")
}

# Initialize data frame for results
results <- data.frame(Year = character(), Mean_PSI = numeric(), Mean_NDVI = numeric(), stringsAsFactors = FALSE)

# Loop through years that exist in both datasets
common_years <- intersect(
  sapply(psi_files, get_year),
  sapply(ndvi_mask_files, get_year)
)

for (year in common_years) {
  psi_path <- psi_files[grepl(year, psi_files)]
  ndvi_path <- ndvi_mask_files[grepl(year, ndvi_mask_files)]
  
  psi_rast <- rast(psi_path)
  ndvi_mask_rast <- rast(ndvi_path)
  
  # Reproject PSI to NDVI
  psi_projected <- project(psi_rast, ndvi_mask_rast, method = "bilinear")
  
  # Mask both rasters
  psi_masked <- mask(psi_projected, ndvi_mask_rast)
  ndvi_masked <- mask(ndvi_mask_rast, ndvi_mask_rast)
  
  # Calculate means
  mean_psi <- global(psi_masked, fun = "mean", na.rm = TRUE)[1, 1]
  mean_ndvi <- global(ndvi_masked, fun = "mean", na.rm = TRUE)[1, 1]
  
  # Append to results
  results <- rbind(results, data.frame(Year = as.numeric(year), Mean_PSI = mean_psi, Mean_NDVI = mean_ndvi))
}

# Calculate correlation
cor_val <- cor(results$Mean_PSI, results$Mean_NDVI, use = "complete.obs")

# Plot
p <- ggplot(results, aes(x = Mean_PSI, y = Mean_NDVI)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "blue") +
  labs(
    x = "Yearly mean soil water potential for pine",
    y = "Yearly mean NDVI quanteils for oak",
    title = "Correlation between mean soil water potential and NDVI quanteils",
    subtitle = paste0("Pearson Correlation: ", round(cor_val, 3))
  ) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.background = element_rect(fill = "white", color = "white"),
    panel.background = element_rect(fill = "white"),
    legend.background = element_rect(fill = "white", color = "white"),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
    plot.subtitle = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_text(color = "black", size = 14),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 14),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    strip.text = element_text(face = "bold", size = 12)
  )

print(p)
ggsave(filename = "results/key_displays/Q_oak_PSI_pine.png", plot = p, device = "png", width = 10, height = 8, dpi = 300)
