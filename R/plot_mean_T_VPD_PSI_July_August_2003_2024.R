setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

# Required package
library(rvest)
library(terra)

# URLs
website_Temp_07 <- "https://opendata.dwd.de/climate_environment/CDC/grids_germany/monthly/air_temperature_mean/07_Jul/"
website_Temp_08 <- "https://opendata.dwd.de/climate_environment/CDC/grids_germany/monthly/air_temperature_mean/08_Aug/"

# Output folder
output_dir <- "../DWD_temp"
dir.create(output_dir, showWarnings = FALSE)

# Function to download all files from a DWD directory listing
download_from_dwd <- function(base_url, output_dir) {
  # Read HTML content from the URL
  page <- read_html(base_url)
  
  # Extract all hrefs from anchor tags
  links <- page %>% html_nodes("a") %>% html_attr("href")
  
  # Filter for files only (e.g., .gz, .zip, .txt)
  file_links <- links[grepl("\\.(zip|gz|txt)$", links)]
  
  # Construct full download URLs
  full_urls <- paste0(base_url, file_links)
  
  # Download each file
  for (i in seq_along(full_urls)) {
    destfile <- file.path(output_dir, file_links[i])
    cat("Downloading", file_links[i], "...\n")
    download.file(full_urls[i], destfile, mode = "wb")
  }
}

# Download from both months
download_from_dwd(website_Temp_07, output_dir)
download_from_dwd(website_Temp_08, output_dir)

library(R.utils)

# Define folder
folder <- "../DWD_temp"

# List all .gz files
gz_files <- list.files(folder, pattern = "\\.gz$", full.names = TRUE)

# Decompress and remove each .gz file
for (gz_file in gz_files) {
  gunzip(gz_file, remove = TRUE, overwrite = TRUE)
  cat("Unzipped and deleted:", basename(gz_file), "\n")
}

asc_files <- list.files(folder, pattern = "\\.asc$", full.names = TRUE)

for (file in asc_files) {
  # Extract year from filename (assumes format like ..._YYYYMM.asc)
  filename <- basename(file)
  year_match <- regmatches(filename, regexpr("\\d{6}", filename))
  year <- as.integer(substr(year_match, 1, 4))
  
  if (is.na(year) || year < 2003 || year > 2024) {
    file.remove(file)
    cat("Deleted file outside range:", filename, "\n")
  }
}

library(terra)

# temperature
Beech.NDVI.2003.mask <- rast("results/Beech/NDVI_mask/NDVI_2003_mask.tif")
Beech.temp.2003.07 <- rast("../DWD_temp/grids_germany_monthly_air_temp_mean_200307.asc")
Beech.temp.2003.08 <- rast("../DWD_temp/grids_germany_monthly_air_temp_mean_200308.asc")
Beech.temp.2003 <- c(Beech.temp.2003.07, Beech.temp.2003.08)
Beech.temp.2003 <- app(Beech.temp.2003, mean)
Beech.temp.2003 <- Beech.temp.2003/10
crs(Beech.temp.2003) <- "epsg:31467"
Beech.temp.2003 <- project(Beech.temp.2003, Beech.NDVI.2003.mask)
Beech.temp.2003 <- mask(Beech.temp.2003, Beech.NDVI.2003.mask)
names(Beech.temp.2003) <- c("Beech.temp.2003")

es.func <- function(tmean)
{
  es <- 6.11 * exp((2.5e6 / 461) * (1 / 273 - 1 / (273 + tmean)))
  return(es)
}

Beech.es.2003 <- app(Beech.temp.2003, es.func)

hy <- rast("../E_OBS_HU/hu_ens_mean_0.1deg_reg_v31.0e.nc")
all_dates <- time(hy)
all_dates <- as.Date(all_dates)
target_dates <- as.Date(c("2003-07-28", "2003-08-29"))
layer_indices <- which(all_dates %in% target_dates)
hy.2003 <- hy[[layer_indices]]
Beech.hy.2003 <- app(hy.2003, mean)
Beech.hy.2003 <- project(Beech.hy.2003, Beech.NDVI.2003.mask)
Beech.hy.2003 <- mask(Beech.hy.2003, Beech.NDVI.2003.mask)

Beech.es.hy.2003 <- c(Beech.es.2003, Beech.hy.2003)

# VPD
vpd.func <- function(es.hy){
  ## calculate saturation vapor pressure
  es <- es.hy[1]
  hy <- es.hy[2]
  ## calculate vapor pressure deficit
  vpd <- ((100 - hy) / 100) * es
  return(vpd)
}

Beech.vpd.2003 <- app(Beech.es.hy.2003, vpd.func)
names(Beech.vpd.2003) <- c("Beech.vpd.2003")
Beech.vpd.2003 <- Beech.vpd.2003/10 # unit: kPa

Beech.temp.vpd.2003 <- c(Beech.temp.2003, Beech.vpd.2003)
Beech.temp.vpd.2003.df <- as.data.frame(Beech.temp.vpd.2003, xy = T)
Beech.temp.vpd.2003.df <- na.omit(Beech.temp.vpd.2003.df)
names(Beech.temp.vpd.2003.df) <- c("x", "y", "temp", "vpd")
Beech.temp.vpd.2003.df$year <- 2003
Beech.temp.vpd.2003.df$species <- "Beech"


library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)

boundary_germany <- ne_countries(scale = "medium", 
                                 country = "Germany", 
                                 returnclass = "sf")

ggplot() +
  # VPD points
  geom_point(data = Beech.temp.vpd.2003.df,
             aes(x = x, y = y, color = vpd),
             size = 0.5) +
  # Germany boundary
  geom_sf(data = boundary_germany,
          fill = NA,
          color = "black",
          inherit.aes = FALSE) +
  # switch to coord_sf
  coord_sf() +
  scale_color_gradientn(
    colours = c("blue", "dodgerblue", "yellow", "orange", "red"),
    name    = "VPD (kPa)"
  ) +
  theme(
    axis.text.x       = element_text(angle = 0, hjust = 0.5),
    plot.background   = element_rect(fill = "white", color = "white"),
    panel.background  = element_rect(fill = "white"),
    legend.background = element_rect(fill = "white", color = "white"),
    plot.title        = element_text(hjust = 0.5, size = 18, face = "bold"),
    plot.subtitle     = element_text(hjust = 0.5, size = 14),
    axis.title        = element_text(face = "bold", size = 16),
    axis.text         = element_text(color = "black", size = 14),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    legend.position   = "top",
    legend.text       = element_text(size = 14),
    strip.background  = element_rect(fill = "white", color = "black", linewidth = 0.5),
    strip.text        = element_text(face = "bold", size = 12)
  )


#–––––––––––––––––––––––––––––––––––––––––
# Compute & Plot Mean Temp & VPD (2003–2024)
#–––––––––––––––––––––––––––––––––––––––––

# 1) Libraries
library(terra)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)
library(tidyr)

# 2) Helper functions
es.func <- function(tmean) {
  6.11 * exp((2.5e6 / 461) * (1/273 - 1/(273 + tmean)))
}
vpd.func <- function(es.hy) {
  es <- es.hy[1]; hy <- es.hy[2]
  ((100 - hy) / 100) * es
}

# 3) Parameters
species_list <- c("Oak","Beech","Spruce","Pine")    # desired order
years <- 2003:2024
hu_file      <- "../E_OBS_HU/hu_ens_mean_0.1deg_reg_v31.0e.nc"

# 4) Containers for each species
mean_temp_maps <- list()
mean_vpd_maps  <- list()

# 5) Loop over species
for (sp in species_list) {
  # read one mask (assume same extent for both years)
  mask_r <- rast(sprintf("species_map_MODIS/%s.tif", sp, years[1]))
  
  tmp_list  <- list()
  vpd_list  <- list()
  
  for (yr in years) {
    # — Temperature (July & August)
    t_j <- rast(sprintf("../DWD_temp/grids_germany_monthly_air_temp_mean_%d07.asc", yr))
    t_a <- rast(sprintf("../DWD_temp/grids_germany_monthly_air_temp_mean_%d08.asc", yr))
    tmp <- app(c(t_j, t_a), mean) / 10
    crs(tmp) <- "epsg:31467"
    tmp <- mask(project(tmp, mask_r), mask_r)
    tmp_list[[as.character(yr)]] <- tmp
    
    # — Saturation vapour pressure
    es  <- app(tmp, es.func)
    
    # — Relative humidity
    hy_all  <- rast(hu_file)
    dates   <- as.Date(time(hy_all))
    tgt     <- as.Date(c(sprintf("%d-07-28",yr),
                         sprintf("%d-08-29",yr)))
    idx     <- which(dates %in% tgt)
    hy      <- app(hy_all[[idx]], mean)
    hy      <- mask(project(hy, mask_r), mask_r)
    
    # — VPD raster
    vpd_r <- app(c(es, hy), vpd.func) / 10
    vpd_list[[as.character(yr)]] <- vpd_r
  }
  
  # — mean across years
  mean_temp_maps[[sp]] <- app(rast(tmp_list), mean, na.rm=TRUE)
  mean_vpd_maps[[sp]]  <- app(rast(vpd_list), mean, na.rm=TRUE)
}

# 6) Convert all species’ rasters into one data frame
df_list <- lapply(species_list, function(sp) {
  # temp
  tdf <- as.data.frame(mean_temp_maps[[sp]], xy=TRUE)
  names(tdf)[3] <- "temp"
  # vpd
  vdf <- as.data.frame(mean_vpd_maps[[sp]], xy=TRUE)
  names(vdf)[3] <- "vpd"
  # merge
  df <- inner_join(tdf, vdf, by=c("x","y"))
  df <- na.omit(df)
  df$species <- sp
  df
})

all_df <- bind_rows(df_list) %>%
  mutate(species = factor(species, levels=species_list))

save(all_df, file = "results/Data/mean_Temp_VPD.RData")

# 7) Germany boundary
boundary_germany <- ne_countries(scale="medium",
                                 country="Germany",
                                 returnclass="sf")

# 8) Plotting

## 8a) Mean Temperature map
library(grid)   # for unit()

p_temp <- ggplot(all_df) +
  geom_point(aes(x = x, y = y, color = temp), size = 0.5) +
  geom_sf(data = boundary_germany, fill = NA, color = "black", inherit.aes = FALSE) +
  scale_color_gradientn(
    colours = c("blue", "dodgerblue", "cyan", "yellow", "orange", "red"),
    name    = "temperature (°C)"
  ) +
  facet_wrap(~species, nrow = 1) +
  coord_sf(expand = FALSE) +
  labs(
    # title    = "Mean Temperature (°C), Jul–Aug 2003–2004",
    # subtitle = "Faceted by species",
    x        = "longitude",
    y        = "latitude"
  ) +
  guides(
    color = guide_colorbar(
      title.position  = "top",        # move title above bar
      title.hjust     = 0.5,          # center the title
      barwidth        = unit(8, "cm"),# make the bar wider
      barheight       = unit(0.4, "cm"),# thin it out a bit
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
    plot.title         = element_text(hjust = 0.5, size = 18, face = "bold"),
    plot.subtitle      = element_text(hjust = 0.5, size = 14),
    strip.background   = element_rect(fill = "white", color = "black"),
    strip.text         = element_text(face = "bold", size = 12)
  )

print(p_temp)
output_path_temp <- "results/key_displays_July_August/mean_temperature_July_August.png"
ggsave(filename = output_path_temp, plot = p_temp, width = 14, height = 6, dpi = 300)

## 8b) Mean VPD map
library(grid)

p_vpd <- ggplot(all_df) +
  geom_point(aes(x = x, y = y, color = vpd), size = 0.5) +
  geom_sf(data = boundary_germany, fill = NA, color = "black", inherit.aes = FALSE) +
  scale_color_gradientn(
    colours = c("blue", "dodgerblue", "cyan", "yellow", "orange", "red"),
    name    = "VPD (kPa)"
  ) +
  facet_wrap(~species, nrow = 1) +
  coord_sf(expand = FALSE) +
  labs(
    x = "longitude",
    y = "latitude"
  ) +
  guides(
    color = guide_colorbar(
      title.position = "top",
      title.hjust    = 0.5,
      barwidth       = unit(8, "cm"),
      barheight      = unit(0.4, "cm"),
      ticks          = TRUE
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
    plot.title         = element_text(hjust = 0.5, size = 18, face = "bold"),
    plot.subtitle      = element_text(hjust = 0.5, size = 14),
    strip.background   = element_rect(fill = "white", color = "black"),
    strip.text         = element_text(face = "bold", size = 12)
  )

print(p_vpd)
output_path_vpd <- "results/key_displays_July_August/mean_vpd_July_August.png"
ggsave(filename = output_path_vpd, plot = p_vpd, width = 14, height = 6, dpi = 300)


## speices map
# 1. enforce the order of species
species_order <- c("Oak", "Beech", "Spruce", "Pine")
cb_palette <- c("Oak"   = "#E69F00",  # Orange
                "Beech" = "#0072B2",  # Deep blue
                "Spruce"= "#009E73",  # Bluish-green
                "Pine"  = "#F0E442")  # Yellow
all_df$species <- factor(all_df$species, levels = species_order)
p_species <- ggplot(all_df, aes(x = x, y = y, color = species)) +
  # 2. points coloured by species
  geom_point(size = 0.5) +
  # overlay Germany boundary
  geom_sf(data = boundary_germany, fill = NA, color = "black", inherit.aes = FALSE) +
  # 3. manual palette
  scale_color_manual(
    values = cb_palette,
    name   = ""
  ) +
  # keep one‐row facets in the same order as the factor
  facet_wrap(~ species, nrow = 1) +
  coord_sf(expand = FALSE) +
  labs(
    x = "longitude",
    y = "latitude"
  ) +
  # 4. custom legend bar
  guides(
    color = guide_legend(
      title.position = "top",
      title.hjust    = 0.5,
      # for discrete, barwidth/barheight get ignored; use keywidth/keyheight
      keywidth  = unit(2, "cm"),
      keyheight = unit(0.8, "cm"),
      nrow      = 1,
      byrow     = TRUE,
      override.aes = list(size = 4)
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
    plot.title         = element_text(hjust = 0.5, size = 18, face = "bold"),
    plot.subtitle      = element_text(hjust = 0.5, size = 14),
    strip.background   = element_rect(fill = "white", color = "black"),
    strip.text         = element_text(face = "bold", size = 12)
  )

# Draw it
print(p_species)
output_path_vpd <- "results/key_displays_July_August/species_map.png"
ggsave(filename = output_path_vpd, plot = p_species, width = 14, height = 6, dpi = 300)

load("results/Data/All_species_month_year_Quantiles_PSI_TDiff.RData")
library(dplyr)
library(ggplot2)
library(grid)     # for unit()
library(sf)       # if boundary_germany is an sf object

# 1. Summarise over all years for July & August
all_df <- final_df %>%
  filter(month %in% c("July", "August")) %>%
  group_by(x, y, species) %>%
  summarise(
    mean_swp = mean(soil_water_potential, na.rm = TRUE),
    mean_td  = mean(transpiration_deficit, na.rm = TRUE)
  ) %>%
  ungroup()

species_order <- c("Oak", "Beech", "Spruce", "Pine")
all_df$species <- factor(all_df$species, levels = species_order)

# 2. Define a helper to build each plot with identical theming
make_map <- function(data, var, legend_title) {
  ggplot(data) +
    geom_point(aes_string(x = "x", y = "y", color = var), size = 0.5) +
    geom_sf(data = boundary_germany, fill = NA, color = "black", inherit.aes = FALSE) +
    scale_color_gradientn(
      colours = c("blue", "dodgerblue", "cyan", "yellow", "orange", "red"),
      name    = legend_title
    ) +
    facet_wrap(~species, nrow = 1) +
    coord_sf(expand = FALSE) +
    labs(x = "longitude", y = "latitude") +
    guides(
      color = guide_colorbar(
        title.position = "top",
        title.hjust    = 0.5,
        barwidth       = unit(8, "cm"),
        barheight      = unit(0.4, "cm"),
        ticks          = TRUE
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
      plot.title         = element_text(hjust = 0.5, size = 18, face = "bold"),
      plot.subtitle      = element_text(hjust = 0.5, size = 14),
      strip.background   = element_rect(fill = "white", color = "black"),
      strip.text         = element_text(face = "bold", size = 12)
    )
}

# 3a. Mean soil‐water potential map
rev_pal <- rev(c("blue", "dodgerblue", "cyan", "yellow", "orange", "red"))
p_swp <- make_map(all_df, "mean_swp", "soil water potential (kPa)") +
  scale_color_gradientn(
    colours = rev_pal,
    limits  = range(all_df$mean_swp, na.rm = TRUE),
    name    = "soil water potential (kPa)"
  )
# p_swp <- make_map(all_df, "mean_swp", "soil water potential (kPa)")
print(p_swp)
ggsave(
  filename = "results/key_displays_July_August/mean_soil_water_potential_July_August.png",
  plot     = p_swp,
  width    = 14,
  height   = 6,
  dpi      = 300
)

# 3b. Mean transpiration deficit map
p_td <- make_map(all_df, "mean_td", "transpiration deficit (mm)")
print(p_td)
ggsave(
  filename = "results/key_displays_July_August/mean_transpiration_deficit_July_August.png",
  plot     = p_td,
  width    = 14,
  height   = 6,
  dpi      = 300
)


library(terra)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)
library(tidyr)

# 2) Helper functions
es.func <- function(tmean) {
  6.11 * exp((2.5e6 / 461) * (1/273 - 1/(273 + tmean)))
}
vpd.func <- function(es.hy) {
  es  <- es.hy[1]
  hy  <- es.hy[2]
  ((100 - hy) / 100) * es
}

# 3) Parameters
species_list <- c("Oak", "Beech", "Spruce", "Pine")  # order matters

hu_file      <- "../E_OBS_HU/hu_ens_mean_0.1deg_reg_v31.0e.nc"

# 4) Containers
mean_temp_maps <- list()
mean_vpd_maps  <- list()
ts_list        <- list()  # for per-year species means

# 5) Loop over species
for (sp in species_list) {
  
  # read mask (same for all years)
  mask_r <- rast(sprintf("results/%s/NDVI_mask/NDVI_%d_mask.tif", sp, years[1]))
  
  # temporary holders
  tmp_list <- list()
  vpd_list <- list()
  recs     <- tibble(year = integer(), mean_temp = numeric(), mean_vpd = numeric())
  
  for (yr in years) {
    # — Temperature (July & August)
    t_j <- rast(sprintf("../DWD_temp/grids_germany_monthly_air_temp_mean_%d07.asc", yr))
    t_a <- rast(sprintf("../DWD_temp/grids_germany_monthly_air_temp_mean_%d08.asc", yr))
    tmp <- app(c(t_j, t_a), mean) / 10
    crs(tmp) <- "epsg:31467"
    tmp <- mask(project(tmp, mask_r), mask_r)
    
    # — Saturation vapour pressure
    es <- app(tmp, es.func)
    
    # — Relative humidity (select the two days)
    hy_all <- rast(hu_file)
    dates  <- as.Date(time(hy_all))
    tgt    <- as.Date(c(sprintf("%d-07-28", yr), sprintf("%d-08-29", yr)))
    idx    <- which(dates %in% tgt)
    hy     <- app(hy_all[[idx]], mean)
    hy     <- mask(project(hy, mask_r), mask_r)
    
    # — VPD raster (kPa)
    vpd_r <- app(c(es, hy), vpd.func) / 10
    
    # store rasters for annual mean‐across‐years maps
    tmp_list[[as.character(yr)]] <- tmp
    vpd_list[[as.character(yr)]] <- vpd_r
    
    # compute global mean for this year & species
    mean_tmp <- terra::global(tmp, "mean", na.rm = TRUE)[1,1]
    mean_vpd <- terra::global(vpd_r, "mean", na.rm = TRUE)[1,1]
    
    recs <- bind_rows(
      recs,
      tibble(year = yr,
             mean_temp = mean_tmp,
             mean_vpd  = mean_vpd)
    )
  }
  
  # mean across years for each pixel
  mean_temp_maps[[sp]] <- app(rast(tmp_list), mean, na.rm = TRUE)
  mean_vpd_maps[[sp]]  <- app(rast(vpd_list), mean, na.rm = TRUE)
  
  # store time series for this species
  recs$species <- sp
  ts_list[[sp]] <- recs
}

# 6) Convert pixel‐level mean rasters into one data frame & save

df_list <- lapply(species_list, function(sp) {
  tdf <- as.data.frame(mean_temp_maps[[sp]], xy = TRUE)
  names(tdf)[3] <- "temp"
  vdf <- as.data.frame(mean_vpd_maps[[sp]], xy = TRUE)
  names(vdf)[3] <- "vpd"
  df  <- inner_join(tdf, vdf, by = c("x","y"))
  na.omit(df) %>%
    mutate(species = sp)
})
all_df <- bind_rows(df_list) %>%
  mutate(species = factor(species, levels = species_list))
# Save spatial means (no year)
save(all_df, file = "results/Data/mean_Temp_VPD.RData")

# Also save species‐level time series (with year)
ts_df <- bind_rows(ts_list) %>%
  mutate(species = factor(species, levels = species_list))
save(ts_df, file = "results/Data/species_Temp_VPD_timeseries.RData")

# 7) (Optional) reshape for plotting if needed
# ts_long <- ts_df %>% pivot_longer(cols = c(mean_temp, mean_vpd), names_to = "variable", values_to = "value")

# 8) Define colour palette and custom theme (as before)
cb_palette <- c(
  "Oak"    = "#E69F00",
  "Beech"  = "#0072B2",
  "Spruce" = "#009E73",
  "Pine"   = "#F0E442"
)

custom_theme <- theme(
  axis.text.x        = element_text(angle = 0, hjust = 0.5),
  plot.background    = element_rect(fill = "white", color = "white"),
  panel.background   = element_rect(fill = "white"),
  legend.background  = element_rect(fill = "white", color = "white"),
  plot.title         = element_text(hjust = 0.5, size = 18, face = "bold", color = "black"),
  plot.subtitle      = element_text(hjust = 0.5, size = 14),
  axis.title         = element_text(face = "bold", size = 16),
  axis.text          = element_text(color = "black", size = 14),
  panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.5),
  panel.grid.major   = element_blank(),
  panel.grid.minor   = element_blank(),
  legend.position    = "top",
  legend.text        = element_text(size = 14),
  strip.background   = element_rect(fill = "white", color = "black", linewidth = 0.5),
  strip.text         = element_text(face = "bold", size = 12)
)

# 9) Plot and save separate time series for Temperature and VPD

# Temperature plot
p_temp <- ggplot(ts_df, aes(x = year, y = mean_temp, color = species)) +
  geom_line(size = 1.2) +
  scale_x_continuous(breaks = seq(from = min(years, na.rm = TRUE),
                                  to   = max(years, na.rm = TRUE),
                                  by   = 5)) +
  scale_color_manual(values = cb_palette) +
  labs(
    x     = "",
    y     = "temperature (°C)",
    color = "",
    title = ""
  ) +
  custom_theme

ggsave(
  filename = "results/key_displays_July_August/mean_temp_time_series.png",
  plot     = p_temp,
  width    = 10,
  height   = 6,
  dpi      = 300
)

# VPD plot
p_vpd <- ggplot(ts_df, aes(x = year, y = mean_vpd, color = species)) +
  geom_line(size = 1.2) +
  scale_x_continuous(breaks = seq(from = min(years, na.rm = TRUE),
                                  to   = max(years, na.rm = TRUE),
                                  by   = 5)) +
  scale_color_manual(values = cb_palette) +
  labs(
    x     = "",
    y     = "VPD (kPa)",
    color = "",
    title = ""
  ) +
  custom_theme

ggsave(
  filename = "results/key_displays_July_August/mean_vpd_time_series.png",
  plot     = p_vpd,
  width    = 10,
  height   = 6,
  dpi      = 300
)

