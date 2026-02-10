# setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0001/yixuan/NDVI_PSI_project")

library(terra)
library(data.table)

# ---------------------------
# Config
# ---------------------------
years   <- 2003:2022
months  <- c("July", "August")
species <- c("Oak", "Beech", "Spruce", "Pine")

base_dir  <- "results_monthly_rootzone"
out_rdata <- "results_rootzone/Data/NDVI_PSI_TDiff_species_month_year_wholeGermany.RData"
dir.create(dirname(out_rdata), showWarnings = FALSE, recursive = TRUE)

NDVI_month <- list(
  July   = "../WZMAllDOYs/Quantiles_209.nc",
  August = "../WZMAllDOYs/Quantiles_241.nc"
)

cat("🚀 Starting pipeline: NDVI + PSI + TDIFF (2003–2022) for July & August\n")

# ---------------------------
# Load NDVI stacks (subset 2003-2022)
# ---------------------------
load_ndvi <- function(nc, years) {
  cat("📥 Loading NDVI from: ", nc, "\n", sep = "")
  r <- rast(nc)
  tt <- time(r)
  if (is.null(tt)) stop("❌ No time() found in NDVI NetCDF: ", nc)
  
  keep <- format(tt, "%Y") %in% as.character(years)
  r <- r[[keep]]
  names(r) <- as.character(time(r))
  
  cat("✅ NDVI loaded. Layers kept: ", nlyr(r), " (", min(years), "–", max(years), ")\n\n", sep = "")
  r
}

NDVI <- list(
  July   = load_ndvi(NDVI_month$July, years),
  August = load_ndvi(NDVI_month$August, years)
)

# ---------------------------
# One common template (NDVI July & August already match)
# ---------------------------
tmpl <- NDVI$July[[1]]

cat("🧭 Checking NDVI July/August grid consistency...\n")
stopifnot(crs(NDVI$August[[1]]) == crs(tmpl))
stopifnot(all(res(NDVI$August[[1]]) == res(tmpl)))
stopifnot(all(ext(NDVI$August[[1]]) == ext(tmpl)))
cat("✅ July & August NDVI share the same CRS/resolution/extent. Using one template.\n\n")

# ---------------------------
# Helpers
# ---------------------------
get_ndvi_layer <- function(ndvi_stack, year) {
  d <- as.Date(names(ndvi_stack))
  idx <- which(format(d, "%Y") == as.character(year))
  if (length(idx) != 1) stop("❌ NDVI layer not found or duplicated for year: ", year)
  ndvi_stack[[idx]]
}

align_to_template <- function(r, tmpl, method = "bilinear") {
  if (!is.na(crs(r)) && crs(r) == crs(tmpl)) {
    resample(r, tmpl, method = method)
  } else {
    project(r, tmpl, method = method)
  }
}

# Remove pixels where ANY layer is NA (complete-case mask)
mask_complete_cases <- function(s) {
  ok <- !is.na(sum(s, na.rm = FALSE))   # NA if any NA
  mask(s, ifel(ok, 1, NA))
}

# ---------------------------
# Extraction for one combo
# ---------------------------
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
    cat("⚠️  Missing files -> ", month, " | ", sp, " | ", year,
        " (psi:", !miss_psi, ", tdiff:", !miss_tdiff, ")\n", sep="")
    return(NULL)
  }
  
  psi   <- rast(psi_path)
  tdiff <- rast(tdiff_path)
  
  # Align PSI & TDIFF to NDVI template
  psi_a   <- align_to_template(psi, tmpl, method = "bilinear")
  tdiff_a <- align_to_template(tdiff, tmpl, method = "bilinear")
  
  # Stack
  s <- c(ndvi_y, psi_a, tdiff_a)
  names(s) <- c("NDVI", "soil_water_potential", "transpiration_deficit")
  
  # Remove pixels with any NA across the 3 layers
  s2 <- mask_complete_cases(s)
  
  # Convert to table
  df <- as.data.frame(s2, xy = TRUE, na.rm = TRUE)
  df$species <- sp
  df$month   <- month
  df$year    <- year
  
  as.data.table(df)
}

# ---------------------------
# Run all combos with progress
# ---------------------------
total <- length(months) * length(species) * length(years)
cat("🧮 Total combinations to process: ", total, "\n\n", sep="")

DT_list <- vector("list", total)
k <- 0
t0_all <- Sys.time()

for (m in months) {
  cat("📅 Month: ", m, " 🌙\n", sep="")
  for (sp in species) {
    cat("🌳 Species: ", sp, "\n", sep="")
    for (yr in years) {
      
      k <- k + 1
      t0 <- Sys.time()
      
      cat(sprintf("  ⏳ [%d/%d] %s | %s | %d ... ", k, total, m, sp, yr))
      
      out <- extract_combo(m, sp, yr)
      
      if (is.null(out)) {
        cat("SKIPPED ⚠️\n")
        DT_list[[k]] <- NULL
      } else {
        DT_list[[k]] <- out
        dt <- round(as.numeric(difftime(Sys.time(), t0, units="secs")), 1)
        cat("DONE ✅ (", dt, "s)\n", sep="")
      }
    }
    cat("✨ Finished species: ", sp, "\n\n", sep="")
  }
  cat("🎉 Finished month: ", m, "\n\n", sep="")
}

cat("🧩 Binding all data tables...\n")
NDVI_PSI_TDiff_df <- rbindlist(DT_list, use.names = TRUE, fill = TRUE)

cat("💾 Saving final RData...\n")
save(NDVI_PSI_TDiff_df, file = out_rdata)

t_all <- round(as.numeric(difftime(Sys.time(), t0_all, units="mins")), 2)
cat("✅ All done! Total time: ", t_all, " minutes 🏁\n", sep="")
cat("📦 Saved here: ", out_rdata, "\n\n", sep="")

if (nrow(missing_log) > 0) {
  miss_path <- file.path(dirname(out_rdata), "missing_files_log.csv")
  fwrite(missing_log, miss_path)
  cat("🧾 Missing file log written to: ", miss_path, "\n", sep="")
} else {
  cat("🥳 No missing PSI/TDIFF files detected.\n")
}

