# Load required libraries
load_libraries <- function() {
  cat("Loading libraries...\n")
  library(dplyr)
  library(sf)
  library(spdep)
  library(ggplot2)
  cat("Libraries loaded.\n\n")
}

# Set working directory and load data
load_data <- function(path, data_file) {
  cat("Setting working directory and loading data...\n")
  setwd(path)
  load(data_file)
  cat("Data loaded.\n\n")
}

# Define species and color palette
get_species_palette <- function() {
  cat("Defining species list and color palette for Pine...\n")
  species_list <- "Pine"
  cb_palette <- c("Pine" = "#009E73")
  cat("Species list and color palette defined.\n\n")
  list(species_list = species_list, cb_palette = cb_palette)
}

# Define ggplot2 theme
create_custom_theme <- function() {
  cat("Defining custom theme...\n")
  theme(
    plot.background = element_rect(fill = "white", color = "white"),
    panel.background = element_rect(fill = "white"),
    legend.background = element_rect(fill = "white", color = "white"),
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold", color = "black"),
    axis.title = element_text(size = 16, face = "bold", color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    axis.text.x = element_text(hjust = 1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13),
    legend.position = "top"
  )
}

# Create required directories
create_result_dirs <- function() {
  cat("Creating results directories if they do not exist...\n")
  dir.create("results_moran", showWarnings = FALSE)
  dir.create("results_moran/data", showWarnings = FALSE, recursive = TRUE)
  cat("Results directories created.\n\n")
}

# Define Moran's I calculation using permutation test
calculate_moran <- function(distance_threshold, values, coords) {
  nb <- dnearneigh(coords, 0, distance_threshold)
  lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
  
  # Monte Carlo permutation test with 999 simulations
  moran <- moran.mc(values, lw, nsim = 999, zero.policy = TRUE)
  
  list(
    distance = distance_threshold,
    observed_I = moran$statistic,
    expected_I = mean(moran$res),
    p_value = moran$p.value
  )
}

# Process species for a given variable
process_species <- function(sp, variable, distances, all_results_df) {
  cat(paste0("  [", sp, "] Processing started...\n"))
  
  # Filter data
  cat("    Filtering species data...\n")
  species_data <- all_results_df %>%
    filter(species == sp, year == "2018")
  cat("    Species data filtered.\n")
  
  # Convert to spatial
  cat("    Converting data to spatial object...\n")
  sp_sf <- st_as_sf(species_data, coords = c("x", "y"), crs = 4326)
  sp_utm <- st_transform(sp_sf, crs = 31467)
  coords <- st_coordinates(sp_utm)
  values <- sp_utm[[variable]]
  cat("    Spatial conversion done.\n")
  
  # Calculate Moran's I
  cat("    Calculating Moran's I for different distances...\n")
  results <- lapply(distances, function(d) {
    calculate_moran(d, values, coords)
  })
  cat("    Moran's I calculation done.\n")
  
  # Build results dataframe
  cat("    Building intermediate results data frame...\n")
  moran_df <- do.call(rbind, lapply(results, function(res) {
    data.frame(
      Distance = res$distance,
      Moran_I = res$observed_I,
      Expected_I = res$expected_I,
      P_value = res$p_value,
      Significance = ifelse(res$p_value < 0.05, "*", ""),
      Species = sp,
      stringsAsFactors = FALSE
    )
  }))
  cat("    Data frame built.\n")
  
  return(moran_df)
}

# Generate and save plot
generate_plot <- function(data, variable, palette, theme) {
  cat("  Creating plot for variable: ", variable, "\n", sep = "")
  plot <- ggplot(data, aes(x = Distance, y = Moran_I, color = Species)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    geom_text(aes(label = Significance), vjust = 2, size = 5) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_color_manual(values = palette) +
    labs(
      title = paste("Moran's I vs Distance for", variable, "in 2018"),
      x = "Distance (meters)",
      y = "Moran's I",
      color = "Species"
    ) +
    theme
  
  print(plot)
  output_file <- paste0("results_moran/Pine_", variable, "_moran_I_vs_distance.png")
  cat("  Saving plot to ", output_file, "...\n", sep = "")
  ggsave(output_file, plot, width = 10, height = 6, dpi = 300)
  cat("  Plot saved.\n\n")
}

# Save results to CSV
save_results <- function(data, variable) {
  file <- paste0("results_moran/data/Pine_", variable, "_moran_data.csv")
  cat("  Saving intermediate data to ", file, "...\n", sep = "")
  write.csv(data, file = file, row.names = FALSE)
  cat("  Data saved.\n")
}

# ---------------- Main Execution ---------------- #

main <- function() {
  # Config
  setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")
  load("results/Data/All_Species_Quantiles_PSI_TDiff.RData")
  variable <- "transpiration_deficit"
  distances <- seq(1000, 20000, by = 4000)
  
  # Run steps
  load_libraries()
  create_result_dirs()
  
  params <- get_species_palette()
  species_list <- params$species_list
  cb_palette <- params$cb_palette
  custom_theme <- create_custom_theme()
  
  cat("Processing variable: ", variable, "\n")
  results_list <- list()
  
  for (sp in species_list) {
    results_list[[sp]] <- process_species(sp, variable, distances, all_results_df)
  }
  
  combined_df <- do.call(rbind, results_list)
  save_results(combined_df, variable)
  generate_plot(combined_df, variable, cb_palette, custom_theme)
  
  cat("All processing completed.\n")
}

# Run the main function
main()
