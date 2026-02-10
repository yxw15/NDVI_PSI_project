library(terra)
library(data.table)

years   <- 2003:2022
months  <- c("August")
species <- c("Oak")

base_dir  <- "results_monthly_rootzone"
out_dir   <- "results_rootzone/Data/chunks"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_rdata <- "results_rootzone/Data/NDVI_PSI_TDiff_species_month_year_wholeGermany.RData"
dir.create(dirname(out_rdata), showWarnings = FALSE, recursive = TRUE)

NDVI_month <- list(
  July   = "../WZMAllDOYs/Quantiles_209.nc",
  August = "../WZMAllDOYs/Quantiles_241.nc"
)

# ---------------------------
# Load NDVI stacks (subset 2003-2022)
# ---------------------------
load_ndvi <- function(nc, years) {
  r <- rast(nc)
  tt <- time(r)
  if (is.null(tt)) stop("No time() found in NDVI NetCDF: ", nc)
  keep <- format(tt, "%Y") %in% as.character(years)
  r <- r[[keep]]
  names(r) <- as.character(time(r))
  r
}

NDVI <- list(
  July   = load_ndvi(NDVI_month$July, years),
  August = load_ndvi(NDVI_month$August, years)
)

tmpl <- NDVI$July[[1]]
stopifnot(crs(NDVI$August[[1]]) == crs(tmpl))
stopifnot(all(res(NDVI$August[[1]]) == res(tmpl)))
stopifnot(all(ext(NDVI$August[[1]]) == ext(tmpl)))

get_ndvi_layer <- function(ndvi_stack, year) {
  d <- as.Date(names(ndvi_stack))
  idx <- which(format(d, "%Y") == as.character(year))
  if (length(idx) != 1) stop("NDVI layer not found or duplicated for year: ", year)
  ndvi_stack[[idx]]
}

align_to_template <- function(r, tmpl, method = "bilinear") {
  if (!is.na(crs(r)) && crs(r) == crs(tmpl)) {
    resample(r, tmpl, method = method)
  } else {
    project(r, tmpl, method = method)
  }
}

mask_complete_cases <- function(s) {
  ok <- !is.na(sum(s, na.rm = FALSE))
  mask(s, ifel(ok, 1, NA))
}

missing_log <- data.table(month=character(), species=character(), year=integer(),
                          missing_psi=logical(), missing_tdiff=logical())

extract_combo <- function(month, sp, year) {
  ndvi_y <- get_ndvi_layer(NDVI[[month]], year)
  
  psi_path   <- file.path(base_dir, month, sp, sprintf("psi_%d.tif", year))
  tdiff_path <- file.path(base_dir, month, sp, sprintf("tdiff_%d.tif", year))
  
  miss_psi   <- !file.exists(psi_path)
  miss_tdiff <- !file.exists(tdiff_path)
  
  if (miss_psi || miss_tdiff) {
    missing_log <<- rbind(
      missing_log,
      data.table(month=month, species=sp, year=year,
                 missing_psi=miss_psi, missing_tdiff=miss_tdiff)
    )
    return(NULL)
  }
  
  psi   <- rast(psi_path)
  tdiff <- rast(tdiff_path)
  
  psi_a   <- align_to_template(psi, tmpl, method="bilinear")
  tdiff_a <- align_to_template(tdiff, tmpl, method="bilinear")
  
  s <- c(ndvi_y, psi_a, tdiff_a)
  names(s) <- c("NDVI", "soil_water_potential", "transpiration_deficit")
  
  s2 <- mask_complete_cases(s)
  
  df <- as.data.frame(s2, xy=TRUE, na.rm=TRUE)
  df$species <- sp
  df$month   <- month
  df$year    <- year
  
  as.data.table(df)
}

# ---------------------------
# Process + checkpoint per (month, species)
# ---------------------------
for (m in months) {
  for (sp in species) {
    
    chunk_file <- file.path(out_dir, sprintf("%s_%s_%d_%d.rds", m, sp, min(years), max(years)))
    
    # Skip if already done
    if (file.exists(chunk_file)) {
      cat("✅ Chunk exists, skipping: ", chunk_file, "\n", sep="")
      next
    }
    
    cat("🚀 Processing chunk: ", m, " | ", sp, "\n", sep="")
    
    DT_years <- vector("list", length(years))
    for (i in seq_along(years)) {
      yr <- years[i]
      cat(sprintf("  ⏳ %s | %s | %d\n", m, sp, yr))
      DT_years[[i]] <- extract_combo(m, sp, yr)
    }
    
    chunk_dt <- rbindlist(DT_years, use.names=TRUE, fill=TRUE)
    
    # Save checkpoint
    saveRDS(chunk_dt, chunk_file, compress = "xz")
    cat("💾 Saved chunk: ", chunk_file, " (rows: ", nrow(chunk_dt), ")\n\n", sep="")
  }
}

# Save missing log (optional)
if (nrow(missing_log) > 0) {
  fwrite(missing_log, file.path(out_dir, "missing_files_log.csv"))
}

cat("✅ Chunking complete.\n")


out_dir <- "results_rootzone/Data/chunks"

# list all chunk files (e.g., August_Beech_2003_2022.rds, July_Oak_2003_2022.rds, ...)
files <- list.files(out_dir, pattern = "\\.rds$", full.names = TRUE)

# fast row-bind
DT_all <- rbindlist(lapply(files, readRDS), use.names = TRUE, fill = TRUE)

# (optional) ensure types are consistent
DT_all[, year := as.integer(year)]
DT_all[, month := as.character(month)]
DT_all[, species := as.character(species)]

DT_all$area <- "Germany"

# save combined
save(DT_all, file = "results_rootzone/Data/NDVI_PSI_TDiff_species_month_year_wholeGermany.RData")
# # or as RDS:
# saveRDS(DT_all, file = "results_rootzone/Data/NDVI_PSI_TDiff_species_month_year_wholeGermany.rds", compress="xz")


load("results/Data/AllSpecies_AllMonths_rootzone.RData")

library(tidyverse)

data <- combined %>% 
  filter(month %in% c("July", "August"),
         Quantiles > 0,
         year <= 2022) %>%
  drop_na()

data$area <- "species"

DT_combined <- bind_rows(
  DT_all,
  data %>%
    rename(NDVI = Quantiles) %>%     # match column name
    select(
      x, y, NDVI, soil_water_potential, transpiration_deficit,
      species, month, year, area
    )                               # drop rootzone + reorder
)

save(DT_combined, file = "results_rootzone/Data/NDVI_PSI_TDiff_species_month_year_area.RData")