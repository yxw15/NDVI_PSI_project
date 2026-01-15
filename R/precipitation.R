# R script to download DWD gridded monthly precipitation for Jul & Aug (2003-2022)
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

library(terra)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)
library(tidyr)

# base URL for monthly precipitation grids
base_url <- "https://opendata.dwd.de/climate_environment/CDC/grids_germany/monthly/precipitation"

# function to create full urls for given years and month codes
create_urls <- function(years, month_code) {
  sprintf("%s/%s/grids_germany_monthly_precipitation_%04d%02d.asc.gz",
          base_url, month_code, years, as.integer(sub(".*_(\\d+)$", "\\1", month_code)))
}

# years and months of interest
years <- 2003:2022
months <- c("07_Jul", "08_Aug")

# create a data frame with all combinations
download_list <- expand.grid(year = years, month = months, stringsAsFactors = FALSE)

# generate the URLs
download_list$url <- with(download_list, 
                          sprintf("%s/%s/grids_germany_monthly_precipitation_%04d%s.asc.gz",
                                  base_url, month,
                                  year, 
                                  ifelse(month == "07_Jul", "07", "08")))

# local download folder
out_dir <- "../DWD_precipitation"
dir.create(out_dir, showWarnings = FALSE)

# download each file if not already present
for (i in seq_len(nrow(download_list))) {
  url <- download_list$url[i]
  destfile <- file.path(out_dir, basename(url))
  if (!file.exists(destfile)) {
    message("Downloading: ", basename(url))
    try(download.file(url, destfile, mode = "wb"), silent = TRUE)
  } else {
    message("Already exists, skipping: ", basename(url))
  }
}

# 1) Decompress Precipitation Files
library(R.utils)
folder_precip <- "../DWD_precipitation"
gz_precip <- list.files(folder_precip, pattern = "\\.gz$", full.names = TRUE)

for (f in gz_precip) {
  gunzip(f, remove = TRUE, overwrite = TRUE)
  cat("Unzipped:", basename(f), "\n")
}

# 2) Parameters & Setup
species_list <- c("Oak", "Beech", "Spruce", "Pine")
years <- 2003:2022
mean_precip_maps <- list()

# 3) Loop over species to process rasters
for (sp in species_list) {
  # Load mask (using your MODIS species maps)
  mask_r <- rast(sprintf("species_map_MODIS/%s.tif", sp))
  
  precip_yr_list <- list()
  
  for (yr in years) {
    # Read July and August Precipitation
    p_j <- rast(sprintf("../DWD_precipitation/grids_germany_monthly_precipitation_%d07.asc", yr))
    p_a <- rast(sprintf("../DWD_precipitation/grids_germany_monthly_precipitation_%d08.asc", yr))
    
    # Calculate mean of Jul/Aug
    # Note: DWD precip is usually in 1/10 mm or mm. Check your metadata. 
    # If it's 1/10 mm like temp, keep the /10. If it's mm, remove it.
    p_mean_month <- app(c(p_j, p_a), mean) 
    
    crs(p_mean_month) <- "epsg:31467"
    p_mean_month <- mask(project(p_mean_month, mask_r), mask_r)
    precip_yr_list[[as.character(yr)]] <- p_mean_month
  }
  
  # Mean across the period 2003-2022
  mean_precip_maps[[sp]] <- app(rast(precip_yr_list), mean, na.rm=TRUE)
}

# 4) Convert to Data Frame
df_list_precip <- lapply(species_list, function(sp) {
  df <- as.data.frame(mean_precip_maps[[sp]], xy=TRUE)
  names(df)[3] <- "precip"
  df <- na.omit(df)
  df$species <- sp
  df
})

precip_df <- bind_rows(df_list_precip) %>%
  mutate(species = factor(species, levels=species_list))

save(precip_df, file = "results/Data/mean_precip_2003_2022.RData")

boundary_germany <- ne_countries(scale="medium",
                                 country="Germany",
                                 returnclass="sf")


# 5) Plotting (Same style as Temperature)
p_precip <- ggplot(precip_df) +
  geom_point(aes(x = x, y = y, color = precip), size = 0.5) +
  geom_sf(data = boundary_germany, fill = NA, color = "black", inherit.aes = FALSE) +
  scale_color_gradientn(
    colours = c("red", "yellow", "cyan", "dodgerblue", "blue"),
    name    = "precipitation (mm)"
  ) +
  facet_wrap(~species, nrow = 1) +
  coord_sf(expand = FALSE) +
  labs(
    x = "longitude",
    y = "latitude"
  ) +
  guides(
    color = guide_colorbar(
      title.position  = "top",
      title.hjust     = 0.5,
      barwidth        = unit(8, "cm"),
      barheight       = unit(0.4, "cm"),
      ticks           = TRUE
    )
  ) +
  theme_minimal() +
  theme(
    axis.text.x        = element_text(hjust = 0.5),
    plot.background    = element_rect(fill = "white", color = NA),
    panel.background   = element_rect(fill = "white"),
    legend.background  = element_rect(fill = "white", color = NA),
    legend.position    = "top",
    legend.text        = element_text(size = 14),
    legend.title       = element_text(size = 14, face = "bold"),
    axis.title         = element_text(face = "bold", size = 16),
    axis.text          = element_text(color = "black", size = 14),
    panel.border       = element_rect(color = "black", fill = NA),
    panel.grid         = element_blank(),
    strip.background   = element_rect(fill = "white", color = "black"),
    strip.text         = element_text(face = "bold", size = 12)
  )

# 6) Save Result
print(p_precip)
output_path_precip <- "results_rootzone/Figures_till2022/supplementary/mean_precipitation_July_August.png"
ggsave(filename = output_path_precip, plot = p_precip, width = 14, height = 6, dpi = 300)

load("results/Data/mean_Temp_VPD.RData")
library(dplyr)

combined_df <- precip_df %>%
  left_join(all_df, by = c("x", "y", "species"))

head(combined_df)

plot_bar_mean_Temp_VPD_Precip <- function(df, output_file) {
  # packages
  require(dplyr)
  require(tidyr)
  require(ggplot2)
  
  # enforce species order & color palette
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  cb_palette <- c(
    "Oak"    = "#E69F00",
    "Beech"  = "#0072B2",
    "Spruce" = "#009E73",
    "Pine"   = "#F0E442"
  )
  
  # compute mean values per species (now incl. precipitation)
  summary_df <- df %>%
    group_by(species) %>%
    summarise(
      Mean_Temp   = mean(temp,   na.rm = TRUE),
      Mean_VPD    = mean(vpd,    na.rm = TRUE),
      Mean_Precip = mean(precip, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_longer(
      cols = c(Mean_Temp, Mean_VPD, Mean_Precip),
      names_to  = "Parameter",
      values_to = "Mean_Value"
    ) %>%
    mutate(
      species   = factor(species, levels = species_order),
      Parameter = dplyr::recode(Parameter,
                                Mean_Temp   = "temperature (°C)",
                                Mean_VPD    = "VPD (kPa)",
                                Mean_Precip = "precipitation (mm)"),
      Parameter = factor(
        Parameter,
        levels = c("temperature (°C)", "VPD (kPa)", "precipitation (mm)")
      )
    )
  
  # build three-panel bar plot
  p <- ggplot(summary_df, aes(x = species, y = Mean_Value, fill = species)) +
    geom_col(color = NA, width = 0.7) +
    facet_wrap(~ Parameter, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = cb_palette) +
    labs(
      x     = "",
      y     = "mean value",
      title = NULL,
      fill  = NULL
    ) +
    theme_minimal() +
    theme(
      plot.background   = element_rect(fill = "white", color = NA),
      panel.background  = element_rect(fill = "white"),
      panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.5),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      axis.title        = element_text(face = "bold", size = 16),
      axis.text         = element_text(color = "black", size = 14),
      legend.position   = "none",
      strip.background  = element_rect(fill = "white", color = "black", size = 0.5),
      strip.text        = element_text(face = "bold", size = 14)
    )
  
  # save & print
  ggsave(filename = output_file, plot = p, width = 14, height = 5, dpi = 300)
  print(p)
  invisible(p)
}

plot_bar_mean_Temp_VPD_Precip(
  combined_df,
  "results_rootzone/Figures_till2022/supplementary/bar_mean_temp_vpd_precip_till2022.png"
)
