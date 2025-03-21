library(terra)

create_species_maps <- function(input_file, output_folder) {
  # Load the raster stack
  Species_Maps <- rast(input_file)
  
  # Ensure the output folder exists
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  
  # Create binary masks for each species
  Oak_mask <- ifel(Species_Maps[[4]] > 0.7, 1, NA)
  Beech_mask <- ifel(Species_Maps[[2]] > 0.7, 2, NA)
  Spruce_mask <- ifel(Species_Maps[[6]] > 0.7, 3, NA)
  Pine_mask <- ifel(Species_Maps[[7]] > 0.7, 4, NA)
  
  # Save the masks to the output folder
  writeRaster(Oak_mask, file.path(output_folder, "Oak_mask.tif"), overwrite = TRUE)
  writeRaster(Beech_mask, file.path(output_folder, "Beech_mask.tif"), overwrite = TRUE)
  writeRaster(Spruce_mask, file.path(output_folder, "Spruce_mask.tif"), overwrite = TRUE)
  writeRaster(Pine_mask, file.path(output_folder, "Pine_mask.tif"), overwrite = TRUE)
  
  cat("Species masks saved successfully in", output_folder, "\n")
}


extract_species_maps <- function(input_raster_path, output_dir, species_values) {
  # Load the raster file
  species_map <- rast(input_raster_path)
  
  # Loop through each species
  for (species_name in names(species_values)) {
    species_value <- species_values[[species_name]]
    
    # Extract cells for the current species
    species_raster <- species_map == species_value
    
    # Define the output file path
    species_dir <- paste0(output_dir, "/", species_name)
    output_file <- paste0(species_dir, "/", species_name, ".tif")
    
    # Create the directory if it doesn't exist
    dir.create(species_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Save the raster with the species name
    writeRaster(species_raster, output_file, overwrite = TRUE)
    
    # Print a success message
    print(paste(species_name, "raster saved as", output_file))
  }
}


aggregate_species <- function(input_files, output_dir, agg_factor = 100) {
  # Define a custom function to calculate the percentage of area with value 1
  percent_cover <- function(values, ...) {
    valid <- !is.na(values)  # Exclude NA values
    sum(values[valid] == 1) / length(values[valid]) * 100
  }
  
  # Loop through each species file
  for (species_name in names(input_files)) {
    input_file <- input_files[[species_name]]
    
    # Load the raster
    species_raster <- rast(input_file)
    
    # Aggregate the raster
    aggregated_raster <- aggregate(species_raster, fact = agg_factor, fun = percent_cover, na.rm = TRUE)
    
    # Define the output file path
    output_file <- paste0(output_dir, "/", species_name, "/", species_name, "_1km_percentage.tif")
    
    # Save the aggregated raster
    writeRaster(aggregated_raster, output_file, overwrite = TRUE)
    
    # Print a success message
    print(paste(species_name, "aggregated raster saved as", output_file))
  }
}

mask_and_reclassify <- function(input_files, threshold = 70) {
  # Loop through each species file
  for (species_name in names(input_files)) {
    input_file <- input_files[[species_name]]
    
    # Load the raster
    species_raster <- rast(input_file)
    
    # Mask and reclassify: keep values > threshold and assign them to 1, others to NA
    masked_raster <- classify(species_raster, rcl = matrix(c(-Inf, threshold, NA, threshold, Inf, 1), ncol = 3, byrow = TRUE))
    
    # Define the output file path (same folder as input file)
    output_file <- gsub("_1km_percentage\\.tif$", "_masked.tif", input_file)
    
    # Save the masked raster
    writeRaster(masked_raster, output_file, overwrite = TRUE)
    
    # Print a success message
    print(paste(species_name, "masked raster saved as", output_file))
  }
}

