# Improved NDVI & PSI Temperature & VPD Processing Script
# Author: Yixuan Wang
# Date: 2025-05-29

# 0) Setup -------------------------------------------------------------------
# Load libraries
library(terra)
library(dplyr)
library(purrr)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(furrr)
library(progressr)

# Global options
default_crs      <- "EPSG:31467"
plan(multisession)
handlers(global = TRUE)

# Paths & parameters --------------------------------------------------------
base_dir         <- "/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project"
setwd(base_dir)
output_dir       <- file.path("..", "DWD_temp")    # assume .asc files present here
years            <- 2003:2015                        # adjust as needed
species_list     <- c("Oak","Beech","Spruce","Pine")
hu_file          <- file.path("..","E_OBS_HU","hu_ens_mean_0.1deg_reg_v31.0e.nc")
boundary_germany <- ne_countries(scale = "medium", country = "Germany", returnclass = "sf")

# 1) Decompress & Filter ---------------------------------------------------
# gz_files <- list.files(output_dir, pattern = "\\.gz$", full.names = TRUE)
# purrr::walk(gz_files, ~ R.utils::gunzip(.x, remove = TRUE, overwrite = TRUE))
# asc_files <- list.files(output_dir, pattern = "\\.asc$", full.names = TRUE)
# keep_files <- asc_files %>%
#   keep(~{
#     ym <- regmatches(basename(.x), regexpr("\\d{6}", basename(.x)))
#     as.integer(substr(ym,1,4)) %in% years
#   })
# unlink(setdiff(asc_files, keep_files))

# 2) Helper functions -------------------------------------------------------
es_func <- function(tmean) {
  6.11 * exp((2.5e6/461) * (1/273 - 1/(273 + tmean)))
}
vpd_func <- function(es, rh) {
  ((100 - rh)/100) * es/10  # convert to kPa
}

# 3) Per-species processing ------------------------------------------------
process_species <- function(sp) {
  mask_file <- sprintf("results/%s/NDVI_mask/NDVI_%d_mask.tif", sp, years[1])
  mask_r    <- rast(mask_file)
  
  ts_recs <- tibble(year=integer(), mean_temp=double(), mean_vpd=double())
  tmp_list <- list(); vpd_list <- list()
  
  for (yr in years) {
    # combine July & August
    files <- c(
      sprintf("%s/grids_germany_monthly_air_temp_mean_%d07.asc", output_dir, yr),
      sprintf("%s/grids_germany_monthly_air_temp_mean_%d08.asc", output_dir, yr)
    )
    tmp    <- rast(files) %>% app(mean, na.rm=TRUE) / 10
    crs(tmp) <- default_crs
    tmp_m  <- mask(project(tmp, mask_r), mask_r)
    es_r   <- app(tmp_m, es_func)
    
    # relative humidity at two dates
    hy_all <- rast(hu_file)
    dts    <- as.Date(time(hy_all))
    idx    <- which(dts %in% as.Date(c(sprintf("%d-07-28",yr), sprintf("%d-08-29",yr))))
    hy_r   <- app(hy_all[[idx]], mean)
    hy_m   <- mask(project(hy_r, mask_r), mask_r)
    
    vpd_r  <- app(c(es_r, hy_m), function(x) vpd_func(x[1], x[2]))
    
    tmp_list[[as.character(yr)]] <- tmp_m
    vpd_list[[as.character(yr)]] <- vpd_r
    
    ts_recs <- ts_recs %>% add_row(
      year      = yr,
      mean_temp = terra::global(tmp_m, "mean", na.rm=TRUE)[1],
      mean_vpd  = terra::global(vpd_r, "mean", na.rm=TRUE)[1]
    )
  }
  
  mean_tmp <- app(rast(tmp_list), mean, na.rm=TRUE)
  mean_vpd <- app(rast(vpd_list), mean, na.rm=TRUE)
  
  list(
    ts     = ts_recs %>% mutate(species = sp),
    raster = list(mean_tmp=mean_tmp, mean_vpd=mean_vpd, species=sp)
  )
}

# 4) Execute in parallel ----------------------------------------------------
results <- future_map(species_list, process_species, .progress=TRUE)
names(results) <- species_list

# 5) Compile time series ----------------------------------------------------
ts_df <- map(results, "ts") %>% bind_rows() %>%
  mutate(species = factor(species, levels=species_list))
save(ts_df, file = "results/Data/species_Temp_VPD_timeseries.RData")

# 6) Compile pixel-level means ----------------------------------------------
pixel_df <- map(results, "raster") %>%
  map_df(~{
    tmp <- .x$mean_tmp
    vpd <- .x$mean_vpd
    sp  <- .x$species
    df  <- as.data.frame(c(tmp, vpd), xy=TRUE)
    names(df)[3:4] <- c("temp","vpd")
    df$species <- sp
    df
  }) %>% drop_na()
save(pixel_df, file = "results/Data/mean_Temp_VPD_global.RData")

# 7) Plotting ----------------------------------------------------------------
# custom_theme unchanged
custom_theme <- theme(
  axis.text.x       = element_text(angle=0,hjust=0.5),
  plot.background   = element_rect(fill="white",color="white"),
  panel.background  = element_rect(fill="white"),
  legend.background = element_rect(fill="white",color="white"),
  plot.title        = element_text(hjust=0.5,size=18,face="bold"),
  plot.subtitle     = element_text(hjust=0.5,size=14),
  axis.title        = element_text(face="bold",size=16),
  axis.text         = element_text(color="black",size=14),
  panel.border      = element_rect(color="black",fill=NA,linewidth=0.5),
  panel.grid.major  = element_blank(),
  panel.grid.minor  = element_blank(),
  legend.position   = "top",
  legend.text       = element_text(size=14),
  strip.background  = element_rect(fill="white",color="black",linewidth=0.5),
  strip.text        = element_text(face="bold",size=12)
)
cb_palette <- c(Oak="#E69F00",Beech="#0072B2",Spruce="#009E73",Pine="#F0E442")

plot_map <- function(df,var,title) {
  ggplot(df)+
    geom_point(aes(x=x,y=y,color=.data[[var]]),size=0.5)+
    geom_sf(data=boundary_germany,fill=NA,color="black",inherit.aes=FALSE)+
    scale_color_gradientn(colours=c("blue","dodgerblue","cyan","yellow","orange","red"),name=title)+
    facet_wrap(~species,nrow=1)+coord_sf(expand=FALSE)+
    labs(x="longitude",y="latitude")+
    guides(color=guide_colorbar(title.position="top",title.hjust=0.5))+custom_theme
}

ggsave("results/key_displays_July_August/mean_temperature.png",
       plot_map(pixel_df,"temp","temperature (°C)"),width=14,height=6)

ggsave("results/key_displays_July_August/mean_vpd.png",
       plot_map(pixel_df,"vpd","VPD (kPa)"),width=14,height=6)

p1 <- ggplot(ts_df,aes(x=year,y=mean_temp,color=species))+geom_line(size=1.2)+
  scale_color_manual(values=cb_palette)+
  labs(x="Year",y="temperature (°C)",color="Species")+
  custom_theme
p2 <- ggplot(ts_df,aes(x=year,y=mean_vpd,color=species))+geom_line(size=1.2)+
  scale_color_manual(values=cb_palette)+
  labs(x="Year",y="VPD (kPa)",color="Species")+
  custom_theme

ggsave("results/key_displays_July_August/mean_temp_time_series.png",p1,width=10,height=6)

ggsave("results/key_displays_July_August/mean_vpd_time_series.png",p2,width=10,height=6)
