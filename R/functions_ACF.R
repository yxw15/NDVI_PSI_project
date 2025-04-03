# Set working directory
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

# Load required package
library(terra)

# Species and variable settings
species_list <- c("Oak", "Beech", "Spruce", "Pine")
vars <- c("soil_water_potential", "transpiration_deficit", "Quantiles")
var_names <- c("PSI", "TDiff", "NDVI")
names(var_names) <- vars

# Color palette for species
cb_palette <- c(
  "Oak"    = "#E69F00",  # Orange
  "Beech"  = "#0072B2",  # Deep blue
  "Spruce" = "#009E73",  # Bluish-green
  "Pine"   = "#F0E442"   # Yellow
)

# Function to compute and plot ACF boxplot across all pixels
plot_acf_for_raster <- function(rast_obj, species, var, output_dir, suffix = "original", cb_palette) {
  vals <- values(rast_obj)
  out_box <- file.path(output_dir, paste0(species, "_", var, "_ACF_all_pixels_boxplot_", suffix, ".png"))
  png(filename = out_box, width = 1200, height = 800)
  
  valid_pixels <- apply(vals, 1, function(x) all(!is.na(x)))
  valid_vals <- vals[valid_pixels, ]
  
  if (nrow(valid_vals) == 0) {
    plot.new()
    title(main = paste("No valid pixels:", species, var, suffix), col.main = cb_palette[species])
    dev.off()
    return()
  }
  
  n_years <- ncol(valid_vals)
  lag_max <- n_years - 1
  
  acf_matrix <- t(apply(valid_vals, 1, function(ts) {
    acf_res <- acf(ts, plot = FALSE, lag.max = lag_max)
    as.vector(acf_res$acf)
  }))
  
  lag_names <- paste0("lag_", 0:lag_max)
  acf_df <- as.data.frame(acf_matrix)
  colnames(acf_df) <- lag_names
  
  ci_bound <- 1.96 / sqrt(n_years)
  
  boxplot(acf_df,
          main = paste("Distribution of", species, gsub("_", " ", var), "ACF (", suffix, ")"),
          xlab = "Lag", ylab = "ACF",
          col = cb_palette[species])
  abline(h = ci_bound, col = "blue", lty = 2)
  abline(h = -ci_bound, col = "blue", lty = 2)
  abline(h = 0, col = "red", lty = 2)
  dev.off()
  
  message("✅ Saved ACF boxplot for ", species, " - ", var, " (", suffix, ")")
}

# Main loop: species × variables
for (species in species_list) {
  for (var in vars) {
    var_short <- var_names[var]
    message("Processing: ", species, " - ", var_short)
    
    # Define folder containing raster masks
    var_folder <- file.path("results", species, paste0(var_short, "_mask"))
    
    # Find all *_<year>_mask.tif files
    tif_files <- list.files(var_folder, pattern = "^.*_\\d{4}_mask\\.tif$", full.names = TRUE)
    
    if (length(tif_files) == 0) {
      message("⚠️ No files found for ", species, " - ", var_short)
      next
    }
    
    # Extract and sort by year
    tif_years <- as.numeric(gsub("^.*_(\\d{4})_mask\\.tif$", "\\1", tif_files))
    sorted_idx <- order(tif_years)
    tif_files <- tif_files[sorted_idx]
    
    # Stack rasters
    var_stack <- rast(tif_files)
    
    # Output directory
    output_dir <- file.path("results_acf")
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    # Plot ACF boxplot
    plot_acf_for_raster(
      rast_obj   = var_stack,
      species    = species,
      var        = var_short,
      output_dir = output_dir,
      suffix     = "original",
      cb_palette = cb_palette
    )
  }
}
