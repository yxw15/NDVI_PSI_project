# Load required packages
library(terra)           # Raster data handling
library(dplyr)           # Data manipulation
library(ggplot2)         # Plotting
library(sf)              # Vector spatial data
library(tidyr)           # Data reshaping
library(rnaturalearth)   # Country boundaries

# Set working directory
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

# Species configuration (with both psi and tdiff)
species_config <- list(
  Beech  = list(
    psi   = "../ALLAN_PIA_SoilMoisture/PSIM/PSImean_rootzone_AMJJA_8days_Bu_bfv_20032024_compressed.nc",
    tdiff = "../ALLAN_PIA_SoilMoisture/TDIFF/TDiffsum_AMJJA_8days_Bu_bfv20032024_compressed.nc"
  ),
  Oak    = list(
    psi   = "../ALLAN_PIA_SoilMoisture/PSIM/PSImean_rootzone_AMJJA_8days_Ei_bfv_20032024_compressed.nc",
    tdiff = "../ALLAN_PIA_SoilMoisture/TDIFF/TDiffsum_AMJJA_8days_Ei_bfv20032024_compressed.nc"
  ),
  Spruce = list(
    psi   = "../ALLAN_PIA_SoilMoisture/PSIM/PSImean_rootzone_AMJJA_8days_Fi_bfv_20032024_compressed.nc",
    tdiff = "../ALLAN_PIA_SoilMoisture/TDIFF/TDiffsum_AMJJA_8days_Fi_bfv20032024_compressed.nc"
  ),
  Pine   = list(
    psi   = "../ALLAN_PIA_SoilMoisture/PSIM/PSImean_rootzone_AMJJA_8days_Ki_bfv_20032024_compressed.nc",
    tdiff = "../ALLAN_PIA_SoilMoisture/TDIFF/TDiffsum_AMJJA_8days_Ki_bfv20032024_compressed.nc"
  )
)

# Color-blind–friendly palette for species
cb_palette <- c(
  Oak    = "#E69F00",
  Beech  = "#0072B2",
  Spruce = "#009E73",
  Pine   = "#F0E442"
)

# Load NDVI template (first layer)
NDVI_template <- rast("../WZMAllDOYs/Quantiles_241.nc")[[1]]

# Fetch Germany boundary from Natural Earth and reproject
germany_border_gk <- ne_countries(
  country     = "Germany",
  scale       = "medium",
  returnclass = "sf"
) %>%
  st_transform(crs(NDVI_template))

# Create output directories
dir.create("results_rootzone/Project_mean", recursive=TRUE, showWarnings=FALSE)
dir.create("results_rootzone/Figures",      recursive=TRUE, showWarnings=FALSE)

# Set target date and year
target_date <- as.Date("2024-07-28")
target_year <- 2024

# Function 1: Process data and save to files
process_variable_data <- function(var_name) {
  message(sprintf("--- Processing data for %s ---", toupper(var_name)))
  dfs <- list()
  
  # 1. Extract & project each species' layer
  for (sp in names(species_config)) {
    nc      <- species_config[[sp]][[var_name]]
    rast_all<- rast(nc)
    idx     <- which(time(rast_all) == target_date)
    if (length(idx) != 1) {
      warning(sprintf("Date %s not found for %s", target_date, sp))
      next
    }
    single <- rast_all[[idx]]
    proj   <- project(single, NDVI_template)
    
    # Save projected raster
    out_tif <- file.path(
      "results_rootzone/Project_mean",
      sprintf("%s_%s_%s.tif", toupper(var_name), sp, format(target_date, "%Y_%m_%d"))
    )
    writeRaster(proj, out_tif, overwrite=TRUE)
    
    # Convert to data.frame
    df <- as.data.frame(proj, xy=TRUE, na.rm=TRUE)
    colnames(df)[3] <- var_name
    df$species <- sp
    df$year    <- target_year
    dfs[[sp]]  <- df
  }
  
  # 2. Combine & save CSV
  all_df <- bind_rows(dfs)
  write.csv(
    all_df,
    file.path("results_rootzone/Project_mean",
              sprintf("%s_points_%s.csv", toupper(var_name), format(target_date, "%Y_%m_%d"))
    ),
    row.names=FALSE
  )
  
  message(sprintf("Data processing completed for %s", toupper(var_name)))
  return(invisible(TRUE))
}

# Function 2: Create value distribution map
plot_distribution <- function(var_name, legend_title) {
  message(sprintf("--- Creating mean distribution plot for %s ---", toupper(var_name)))
  
  # Load existing CSV file
  csv_file <- file.path(
    "results_rootzone/Project_mean",
    sprintf("%s_points_%s.csv", toupper(var_name), format(target_date, "%Y_%m_%d"))
  )
  
  if (!file.exists(csv_file)) {
    warning(sprintf("CSV file not found: %s", csv_file))
    return(NULL)
  }
  
  all_df <- read.csv(csv_file)
  
  # Prepare data for plotting
  dfm <- all_df %>%
    group_by(x, y, species, year) %>%
    summarise(mean_val = mean(.data[[var_name]], na.rm=TRUE), .groups="drop") %>%
    mutate(species = factor(species, levels=c("Oak","Beech","Spruce","Pine")))
  
  # Define color palette
  pal <- if (var_name=="psi") rev(c("blue","dodgerblue","cyan","yellow","orange","red"))
  else c("blue","dodgerblue","cyan","yellow","orange","red")
  
  # Create plot
  p_dist <- ggplot(dfm, aes(x=x, y=y, color=mean_val)) +
    geom_point(size=0.5) +
    geom_sf(data=germany_border_gk, fill=NA, color="black", inherit.aes=FALSE) +
    scale_color_gradientn(colours=pal, name=legend_title) +
    facet_grid(species ~ year) +
    coord_sf(crs=crs(NDVI_template), expand=FALSE) +
    labs(x="longitude", y="latitude") +
    theme_minimal() +
    theme(
      legend.position  = "top",
      axis.title       = element_text(face="bold", size=16),
      axis.text        = element_text(size=14, color="black"),
      strip.text       = element_text(face="bold", size=12),
      panel.border     = element_rect(color="black", fill=NA),
      panel.grid       = element_blank()
    )
  
  # Save plot
  ggsave(
    file.path("results_rootzone/Figures",
              sprintf("mean_%s_all_species_%s.png", var_name, format(target_date, "%Y_%m_%d"))
    ),
    p_dist, width=14, height=6, dpi=300
  )
  
  message(sprintf("Mean distribution plot completed for %s", toupper(var_name)))
  return(invisible(TRUE))
}

# Function 3: Create faceted scatter plot (Oak vs others)
plot_faceted_scatter <- function(var_name) {
  message(sprintf("--- Creating faceted scatter plot for %s ---", toupper(var_name)))
  
  # Load TIFF files
  tifs <- list.files(
    "results_rootzone/Project_mean",
    pattern    = sprintf("^%s_.*_%s\\.tif$", toupper(var_name), format(target_date, "%Y_%m_%d")),
    full.names = TRUE
  )
  
  if (length(tifs) == 0) {
    warning(sprintf("No TIFF files found for %s", var_name))
    return(NULL)
  }
  
  names(tifs) <- sub(
    sprintf("^%s_(.*)_%s\\.tif$", toupper(var_name), format(target_date, "%Y_%m_%d")),
    "\\1", basename(tifs)
  )
  
  # Sample random points
  st <- rast(tifs)
  set.seed(42)
  samp <- sampleRandom(st, 10000, xy=TRUE, na.rm=TRUE)
  df_s <- as.data.frame(samp)
  
  # Prepare data for plotting
  df_long <- df_s %>%
    pivot_longer(
      cols = setdiff(names(df_s), c("x","y","Oak")),
      names_to = "species",
      values_to = "value"
    )
  
  # Create faceted scatter plot
  p_facet <- ggplot(df_long, aes(x=value, y=Oak, color=species)) +
    geom_point(alpha=0.5, size=1) +
    facet_wrap(~species, ncol=1) +
    scale_color_manual(values=cb_palette) +
    labs(x=sprintf("Other %s", toupper(var_name)),
         y=sprintf("Oak %s", toupper(var_name))) +
    theme_minimal() +
    theme(
      legend.position  = "top",
      axis.title       = element_text(face="bold", size=16),
      axis.text        = element_text(size=14, color="black"),
      strip.text       = element_text(face="bold", size=12),
      panel.border     = element_rect(color="black", fill=NA),
      panel.grid       = element_blank()
    )
  
  # Save plot
  ggsave(
    file.path("results_rootzone/Figures",
              sprintf("scatter_facet_Oak_vs_others_%s.png", var_name)
    ),
    p_facet, width=5, height=15, dpi=300
  )
  
  message(sprintf("Faceted scatter plot completed for %s", toupper(var_name)))
  return(invisible(TRUE))
}

# Function 4: Create overlay scatter plot (Oak vs others)
plot_overlay_scatter <- function(var_name) {
  message(sprintf("--- Creating overlay scatter plot for %s ---", toupper(var_name)))
  
  # Load TIFF files
  tifs <- list.files(
    "results_rootzone/Project_mean",
    pattern    = sprintf("^%s_.*_%s\\.tif$", toupper(var_name), format(target_date, "%Y_%m_%d")),
    full.names = TRUE
  )
  
  if (length(tifs) == 0) {
    warning(sprintf("No TIFF files found for %s", var_name))
    return(NULL)
  }
  
  names(tifs) <- sub(
    sprintf("^%s_(.*)_%s\\.tif$", toupper(var_name), format(target_date, "%Y_%m_%d")),
    "\\1", basename(tifs)
  )
  
  # Sample random points
  st <- rast(tifs)
  set.seed(42)
  samp <- sampleRandom(st, 10000, xy=TRUE, na.rm=TRUE)
  df_s <- as.data.frame(samp)
  
  # Prepare data for plotting
  df_long <- df_s %>%
    pivot_longer(
      cols = setdiff(names(df_s), c("x","y","Oak")),
      names_to = "species",
      values_to = "value"
    )
  
  # Create overlay scatter plot
  p_overlay <- ggplot(df_long, aes(x=value, y=Oak, color=species)) +
    geom_point(alpha=0.5, size=1) +
    scale_color_manual(values=cb_palette) +
    labs(x=sprintf("Other %s", toupper(var_name)),
         y=sprintf("Oak %s", toupper(var_name))) +
    theme_minimal() +
    theme(
      legend.position  = "top",
      axis.title       = element_text(face="bold", size=16),
      axis.text        = element_text(size=14, color="black"),
      panel.border     = element_rect(color="black", fill=NA),
      panel.grid       = element_blank()
    )
  
  # Save plot
  ggsave(
    file.path("results_rootzone/Figures",
              sprintf("scatter_overlay_Oak_vs_others_%s.png", var_name)
    ),
    p_overlay, width=6, height=6, dpi=300
  )
  
  message(sprintf("Overlay scatter plot completed for %s", toupper(var_name)))
  return(invisible(TRUE))
}

# Function 5: Master plotting function that calls all individual plot functions
plot_variable <- function(var_name, legend_title) {
  message(sprintf("--- Creating all plots for %s ---", toupper(var_name)))
  
  # Create all plot types
  plot_mean_distribution(var_name, legend_title)
  plot_faceted_scatter(var_name)
  plot_overlay_scatter(var_name)
  
  message(sprintf("All plots completed for %s", toupper(var_name)))
  return(invisible(TRUE))
}

# OPTION 1: Process data first (if files don't exist yet)
process_variable_data("psi")
process_variable_data("tdiff")

# OPTION 2: Create all plots
plot_variable("psi", "soil water potential (kPa)")
plot_variable("tdiff", "transpiration deficit (mm)")

# OPTION 3: Create individual plots (if you want to run them separately)
plot_mean_distribution("psi", "soil water potential (kPa)")
plot_faceted_scatter("psi")
plot_overlay_scatter("psi")
plot_mean_distribution("tdiff", "transpiration deficit (mm)")
plot_faceted_scatter("tdiff")
plot_overlay_scatter("tdiff")