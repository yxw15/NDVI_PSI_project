### Updated log-binned functions with absolute-value transform and renamed plot function

#' Bin NDVI by log-transformed soil water potential magnitude
NDVI_PSIbin_log <- function(df, bin_width_log = 0.5) {
  library(dplyr)
  
  # Identify correct value column
  value_column <- if ("Quantiles" %in% names(df)) "Quantiles" else "Proportions"
  
  # Filter to valid, non-zero potentials
  df <- df %>% filter(!is.na(soil_water_potential), soil_water_potential != 0)
  
  # Log-transform magnitude of soil water potential
  df <- df %>% mutate(log_PSI = log(abs(soil_water_potential)))
  
  # Total pixel count per species
  species_totals <- df %>%
    group_by(species) %>%
    summarise(total_pixels = n(), .groups = "drop")
  
  # Define log-scale bin breaks
  psi_min <- floor(min(df$log_PSI, na.rm = TRUE))
  psi_max <- ceiling(max(df$log_PSI, na.rm = TRUE))
  bin_breaks <- seq(psi_min, psi_max, by = bin_width_log)
  
  # Bin log(soil_water_potential magnitude)
  df <- df %>%
    mutate(PSI_bin = cut(log_PSI, breaks = bin_breaks, include.lowest = TRUE, right = FALSE))
  
  # Compute stats per species and bin
  result <- df %>%
    group_by(species, PSI_bin) %>%
    summarise(
      avg_value = mean(.data[[value_column]], na.rm = TRUE),
      count = n(),
      .groups = 'drop'
    ) %>%
    mutate(
      bin_median = sapply(as.character(PSI_bin), function(bin_label) {
        edges <- as.numeric(strsplit(gsub("\\[|\\]|\\(|\\)", "", bin_label), ",")[[1]])
        mean(edges)
      })
    ) %>%
    left_join(species_totals, by = "species") %>%
    mutate(percentage = count / total_pixels) %>%
    filter(percentage >= 0.01) %>%
    select(species, PSI_bin, bin_median, avg_value, count, total_pixels, percentage)
  
  return(result)
}

#' Bin NDVI by log-transformed transpiration deficit
NDVI_TDiffbin <- function(df, bin_width_log = 0.2) {
  library(dplyr)
  
  # Identify correct value column
  value_column <- if ("Quantiles" %in% names(df)) "Quantiles" else "Proportions"
  
  # Filter to positive deficits
  df <- df %>% filter(!is.na(transpiration_deficit), transpiration_deficit > 0)
  
  # Log-transform transpiration deficit
  df <- df %>% mutate(log_TDiff = log(transpiration_deficit))
  
  # Total pixel count per species
  species_totals <- df %>%
    group_by(species) %>%
    summarise(total_pixels = n(), .groups = "drop")
  
  # Define log-scale bin breaks
  tdiff_min <- floor(min(df$log_TDiff, na.rm = TRUE))
  tdiff_max <- ceiling(max(df$log_TDiff, na.rm = TRUE))
  bin_breaks <- seq(tdiff_min, tdiff_max, by = bin_width_log)
  
  # Bin
  df <- df %>%
    mutate(TDiff_bin = cut(log_TDiff, breaks = bin_breaks, include.lowest = TRUE, right = FALSE))
  
  # Compute stats
  result <- df %>%
    group_by(species, TDiff_bin) %>%
    summarise(
      avg_value = mean(.data[[value_column]], na.rm = TRUE),
      count = n(),
      .groups = 'drop'
    ) %>%
    mutate(
      bin_median = sapply(as.character(TDiff_bin), function(bin_label) {
        edges <- as.numeric(strsplit(gsub("\\[|\\]|\\(|\\)", "", bin_label), ",")[[1]])
        mean(edges)
      })
    ) %>%
    left_join(species_totals, by = "species") %>%
    mutate(percentage = count / total_pixels) %>%
    filter(percentage >= 0.0001) %>%
    select(species, TDiff_bin, bin_median, avg_value, count, total_pixels, percentage)
  
  return(result)
}

#' Bin log-transformed transpiration deficit by log-transformed soil water potential magnitude
TDiff_PSIbin <- function(df, bin_width_log = 0.5) {
  library(dplyr)
  
  # Filter valid rows
  df <- df %>%
    filter(!is.na(soil_water_potential), soil_water_potential != 0,
           !is.na(transpiration_deficit), transpiration_deficit > 0)
  
  # Log-transform both
  df <- df %>%
    mutate(
      log_PSI = log(abs(soil_water_potential)),
      log_TDiff = log(transpiration_deficit)
    )
  
  # Total pixel count per species
  species_totals <- df %>%
    group_by(species) %>%
    summarise(total_pixels = n(), .groups = "drop")
  
  # Define PSI bin breaks
  psi_min <- floor(min(df$log_PSI, na.rm = TRUE))
  psi_max <- ceiling(max(df$log_PSI, na.rm = TRUE))
  bin_breaks <- seq(psi_min, psi_max, by = bin_width_log)
  
  df <- df %>%
    mutate(PSI_bin = cut(log_PSI, breaks = bin_breaks, include.lowest = TRUE, right = FALSE))
  
  # Compute stats
  result <- df %>%
    group_by(species, PSI_bin) %>%
    summarise(
      avg_transpiration_deficit = mean(log_TDiff, na.rm = TRUE),
      count = n(),
      .groups = 'drop'
    ) %>%
    mutate(
      bin_median = sapply(as.character(PSI_bin), function(bin_label) {
        edges <- as.numeric(strsplit(gsub("\\[|\\]|\\(|\\)", "", bin_label), ",")[[1]])
        mean(edges)
      })
    ) %>%
    left_join(species_totals, by = "species") %>%
    mutate(percentage = count / total_pixels) %>%
    filter(percentage > 0.001) %>%
    select(species, PSI_bin, bin_median, avg_transpiration_deficit, count, total_pixels, percentage)
  
  return(result)
}

#' Plot combined AIC and R2 using log-binned variables
plot_log_AIC_R2 <- function(data, save_combined_fig) {
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  
  species_order <- c("Oak", "Beech", "Spruce", "Pine")
  
  # Panel A: NDVI ~ log(PSI)
  data_a <- NDVI_PSIbin(data)
  data_a <- na.omit(data_a)
  data_a$species <- factor(data_a$species, levels = species_order)
  data_a <- data_a %>% mutate(x = -bin_median)
  
  aic_results_a <- lapply(levels(data_a$species), function(sp) {
    sp_data <- data_a[data_a$species == sp,]
    lm_l <- lm(avg_value ~ x, data = sp_data)
    nls_e <- tryCatch(nls(avg_value ~ a + b * exp(-c * x), data = sp_data,
                          start = list(a=5, b=3, c=0.001),
                          control = nls.control(maxiter=200)), error = function(e) NULL)
    aic_vals <- c(AIC(lm_l), if(!is.null(nls_e)) AIC(nls_e) else NA)
    r2_vals  <- c(summary(lm_l)$r.squared,
                  if(!is.null(nls_e)) {res <- resid(nls_e); 1 - sum(res^2)/sum((sp_data$avg_value - mean(sp_data$avg_value))^2)} else NA)
    data.frame(species=sp, Model=c("linear","exponential"), AIC=aic_vals, R2=r2_vals)
  })
  aic_df_a <- bind_rows(aic_results_a)
  
  p_a <- ggplot(aic_df_a, aes(x=species, y=AIC, fill=Model)) +
    geom_bar(stat="identity", position=position_dodge(0.9)) +
    geom_text(aes(label=round(R2,2), y=AIC/2), position=position_dodge(0.9), vjust=0.5, size=8) +
    labs(title="NDVI ~ log(soil water potential)") + theme_minimal() + theme(legend.position="top")
  
  # Panel B: NDVI ~ log(TDiff)
  data_b <- NDVI_TDiffbin(data)
  data_b <- na.omit(data_b)
  data_b$species <- factor(data_b$species, levels = species_order)
  data_b <- data_b %>% mutate(x = bin_median)
  
  aic_results_b <- lapply(levels(data_b$species), function(sp) {
    sp_data <- data_b[data_b$species == sp,]
    lm_l <- lm(avg_value ~ x, data = sp_data)
    nls_e <- tryCatch(nls(avg_value ~ a + b * exp(-c * x), data = sp_data,
                          start = list(a=5,b=7,c=0.04),
                          control = nls.control(maxiter=1200)), error = function(e) NULL)
    aic_vals <- c(AIC(lm_l), if(!is.null(nls_e)) AIC(nls_e) else NA)
    r2_vals  <- c(summary(lm_l)$r.squared,
                  if(!is.null(nls_e)) {res <- resid(nls_e); 1 - sum(res^2)/sum((sp_data$avg_value - mean(sp_data$avg_value))^2)} else NA)
    data.frame(species=sp, Model=c("linear","exponential"), AIC=aic_vals, R2=r2_vals)
  })
  aic_df_b <- bind_rows(aic_results_b)
  
  p_b <- ggplot(aic_df_b, aes(x=species, y=AIC, fill=Model)) +
    geom_bar(stat="identity", position=position_dodge(0.9)) +
    geom_text(aes(label=round(R2,2), y=AIC/2), position=position_dodge(0.9), vjust=0.5, size=8) +
    labs(title="NDVI ~ log(transpiration deficit)") + theme_minimal() + theme(legend.position="top")
  
  # Panel C: log(TDiff) ~ log(PSI)
  data_c <- TDiff_PSIbin(data)
  data_c <- na.omit(data_c)
  data_c$species <- factor(data_c$species, levels = species_order)
  data_c <- data_c %>% mutate(x=bin_median)
  
  aic_results_c <- lapply(levels(data_c$species), function(sp) {
    sp_data <- data_c[data_c$species == sp,]
    models <- list(
      linear = lm(avg_transpiration_deficit ~ x, data = sp_data),
      poly2  = lm(avg_transpiration_deficit ~ x + I(x^2), data = sp_data),
      poly3  = lm(avg_transpiration_deficit ~ x + I(x^2) + I(x^3), data = sp_data)
    )
    aic_vals <- sapply(models, AIC)
    r2_vals  <- sapply(models, function(m) summary(m)$r.squared)
    data.frame(species=sp, Model=names(models), AIC=aic_vals, R2=r2_vals)
  })
  aic_df_c <- bind_rows(aic_results_c)
  
  p_c <- ggplot(aic_df_c, aes(x=species, y=AIC, fill=Model)) +
    geom_bar(stat="identity", position=position_dodge(0.9)) +
    geom_text(aes(label=round(R2,2), y=AIC/2), position=position_dodge(0.9), vjust=0.5, size=8) +
    labs(title="log(transpiration deficit) ~ log(soil water potential)") + theme_minimal() + theme(legend.position="top")
  
  # Combine and save
  combined_top <- (p_a | p_b) + patchwork::plot_layout(guides="collect") & theme(legend.position="top")
  combined_plot <- combined_top / p_c
  dir.create(dirname(save_combined_fig), recursive=TRUE, showWarnings=FALSE)
  ggsave(filename=save_combined_fig, plot=combined_plot, width=14, height=12, dpi=300)
  print(combined_plot)
}

plot_log_AIC_R2(df,  "results/key_displays_July_August/AIC_log.png")

setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")

library(dplyr)
load("results/Data/All_species_month_year_Quantiles_PSI_TDiff.RData")
df <- final_df %>% filter(month %in% c("July", "August"))

NDVI_PSIbin_df <- NDVI_PSIbin(df)
NDVI_PSIbin_df <- na.omit(NDVI_PSIbin_df)

# assign one color per species
cols <- rainbow(length(unique(NDVI_PSIbin_df$species)))
names(cols) <- unique(NDVI_PSIbin_df$species)

with(NDVI_PSIbin_df, 
     plot(-bin_median+8, avg_value,
          col   = cols[species],
          pch   = 19,
          xlab  = "Bin Median",
          ylab  = "Average Value",
          main  = "Average Value vs. Bin Median"))

legend("topleft", legend = names(cols), col = cols, pch = 19, title = "Species")

library(gamm4)

g4 <- gamm4(
  avg_value ~ s(bin_median, k = 10),
  random   = ~(1|species),
  weights  = count,
  data     = NDVI_PSIbin_df
)
summary(g4$gam)      # smooth component
summary(g4$mer)      # mixed‐model component

library(ggplot2)

ggplot(NDVI_PSIbin_df, 
       aes(x = bin_median, 
           y = avg_value, 
           size = count)) +
  geom_point(alpha = 0.6) +
  # fit a separate GAM for each species, weighting by count
  geom_smooth(
    aes(weight = count),
    method  = "gam",
    formula = y ~ s(x, k = 10),
    se      = FALSE
  ) +
  facet_wrap(~ species, scales = "free_x") +
  scale_size_area(max_size = 4) +
  labs(
    x     = "Bin median PSI",
    y     = "Average NDVI",
    title = "NDVI vs PSI by Species with GAM fit"
  ) +
  theme_minimal()


# 1. Load required libraries
library(gamm4)     # for GAMM fitting
library(dplyr)     # for data manipulation
library(tidyr)     # for expand_grid()
library(ggplot2)   # for plotting
library(tibble)    # for rownames_to_column()

# 2. Fit the GAMM
#    - avg_value ~ smooth(bin_median)
#    - random intercept by species
#    - weights = count
g4 <- gamm4(
  avg_value ~ s(bin_median, k = 10),
  random   = ~(1 | species),
  weights  = count,
  data     = NDVI_PSIbin_df
)

# 3. Inspect both components
summary(g4$gam)  # smooth component
summary(g4$mer)  # mixed‐model (random effects) component


# 4. Specify species order & custom palette
species_order <- c("Oak", "Beech", "Spruce", "Pine")
cb_palette <- c(
  "Oak"    = "#E69F00",
  "Beech"  = "#0072B2",
  "Spruce" = "#009E73",
  "Pine"   = "#F0E442"
)

# 5. Re‐factor your original data
NDVI_PSIbin_df <- NDVI_PSIbin_df %>%
  mutate(species = factor(species, levels = species_order))

# 6. Build a fine prediction grid
newd <- expand_grid(
  species    = species_order,
  bin_median = seq(
    min(NDVI_PSIbin_df$bin_median),
    max(NDVI_PSIbin_df$bin_median),
    length = 200
  )
) %>%
  mutate(species = factor(species, levels = species_order)) %>%
  # 7. Add the fixed‐effect (smooth) prediction
  mutate(fit_fixed = predict(g4$gam, newdata = .))

# 8. Extract species random intercepts
ran_int <- ranef(g4$mer)$species %>%
  as.data.frame() %>%
  rownames_to_column(var = "species") %>%
  rename(rand_int = `(Intercept)`) %>%
  mutate(species = factor(species, levels = species_order))

# 9. Combine fixed + random to get full prediction
newd <- newd %>%
  left_join(ran_int, by = "species") %>%
  mutate(fit = fit_fixed + rand_int)


# 10. Define your custom ggplot2 theme
my_theme <- theme(
  axis.text.x           = element_text(angle = 0, hjust = 0.5),
  plot.background       = element_rect(fill = "white", color = "white"),
  panel.background      = element_rect(fill = "white"),
  legend.background     = element_rect(fill = "white", color = "white"),
  plot.title            = element_text(hjust = 0.5, size = 18, face = "bold"),
  plot.caption          = element_text(face = "bold", size = 16, hjust = 0),
  plot.caption.position = "plot",
  axis.title            = element_text(face = "bold", size = 16),
  axis.text             = element_text(color = "black", size = 14),
  panel.border          = element_blank(),
  panel.grid.major      = element_blank(),
  panel.grid.minor      = element_blank(),
  legend.position       = "none",
  legend.text           = element_text(size = 14)
)

# 11. Plot raw points + GAMM fits
ggplot(NDVI_PSIbin_df, aes(x = bin_median, y = avg_value, color = species)) +
  geom_point(aes(size = count), alpha = 0.6) +
  geom_line(
    data = newd,
    aes(y = fit, color = species),
    size = 1
  ) +
  scale_color_manual(values = cb_palette) +
  scale_size_area(max_size = 4) +
  labs(
    x     = "Bin-median PSI",
    y     = "Average NDVI",
    title = "NDVI vs PSI by Species (GAMM fit)"
  ) +
  my_theme

# 1. Load libraries
library(lme4)      # for lmer()
library(dplyr)     # for data manipulation
library(tidyr)     # for expand_grid()
library(ggplot2)   # for plotting
library(tibble)    # for rownames_to_column()

# 2. Prepare data: positive X and log-transform

df <- NDVI_PSIbin_df %>%
  # make X positive
  mutate(
    Xpos   = -bin_median,
    logX   = log(Xpos),
    species = factor(species, levels = c("Oak","Beech","Spruce","Pine"))
  )

# 3. Fit a linear mixed model
#    avg_value ~ logX + (1 | species), weighted by count
lmm <- lmer(
  avg_value ~ logX + (1 | species),
  weights = count,
  data    = df
)
summary(lmm)

# 4. Build prediction grid for each species
pred_grid <- expand_grid(
  species = levels(df$species),
  Xpos    = seq(min(df$Xpos), max(df$Xpos), length = 200)
) %>%
  mutate(
    logX = log(Xpos)
  )

# 5. Get fitted values including random intercepts
#    re.form = NULL includes both fixed and random effects
pred_grid$fit <- predict(lmm, newdata = pred_grid, re.form = NULL)

# 6. Convert back to original x-axis for plotting
pred_grid <- pred_grid %>%
  mutate(bin_median = -Xpos)

# 7. Custom palette & theme
cb_palette <- c(
  "Oak"    = "#E69F00",
  "Beech"  = "#0072B2",
  "Spruce" = "#009E73",
  "Pine"   = "#F0E442"
)
my_theme <- theme(
  axis.text.x           = element_text(angle = 0, hjust = 0.5),
  plot.background       = element_rect(fill = "white", color = "white"),
  panel.background      = element_rect(fill = "white"),
  legend.background     = element_rect(fill = "white", color = "white"),
  plot.title            = element_text(hjust = 0.5, size = 18, face = "bold"),
  plot.caption          = element_text(face = "bold", size = 16, hjust = 0),
  plot.caption.position = "plot",
  axis.title            = element_text(face = "bold", size = 16),
  axis.text             = element_text(color = "black", size = 14),
  panel.border          = element_blank(),
  panel.grid.major      = element_blank(),
  panel.grid.minor      = element_blank(),
  legend.position       = "none",
  legend.text           = element_text(size = 14)
)

# 8. Plot raw points + LMM fits
ggplot(df, aes(x = bin_median, y = avg_value, color = species)) +
  geom_point(aes(size = count), alpha = 0.6) +
  geom_line(
    data = pred_grid,
    aes(x = bin_median, y = fit, color = species),
    size = 1
  ) +
  scale_color_manual(values = cb_palette) +
  scale_size_area(max_size = 4) +
  labs(
    x     = "Bin-median PSI",
    y     = "Average NDVI",
    title = "Linear Mixed Model: NDVI vs log(–PSI) by Species"
  ) +
  my_theme


# 1. Load required libraries
library(lme4)      # for lmer()
library(dplyr)     # for data manipulation
library(tidyr)     # for expand_grid()
library(ggplot2)   # for plotting
library(tibble)    # for rownames_to_column()

# 2. Specify species order and palette
species_order <- c("Oak", "Beech", "Spruce", "Pine")
cb_palette <- c(
  "Oak"    = "#E69F00",
  "Beech"  = "#0072B2",
  "Spruce" = "#009E73",
  "Pine"   = "#F0E442"
)

# 3. Prepare the data
#    - Shift: Xpos = –bin_median + 8
#    - Log transform: logX = log(Xpos)
#    - Factor species in the desired order
df <- NDVI_PSIbin_df %>%
  mutate(
    species = factor(species, levels = species_order),
    Xpos    = -bin_median + 8,
    logX    = log(Xpos)
  )

# (Optional) check that Xpos is strictly positive:
if(any(df$Xpos <= 0)) {
  stop("Some Xpos values are ≤ 0. Choose a larger shift than +8.")
}

# 4. Fit the linear mixed model
#    avg_value ~ logX + (1 | species), weighted by count
lmm2 <- lmer(
  avg_value ~ logX + (1 | species),
  weights = count,
  data    = df
)
summary(lmm2)


# 5. Build a prediction grid
pred_grid2 <- expand_grid(
  species    = species_order,
  bin_median = seq(
    min(df$bin_median),
    max(df$bin_median),
    length = 200
  )
) %>%
  mutate(
    # same transformations as for df:
    Xpos = -bin_median + 8,
    logX = log(Xpos)
  )

# 6. Get model‐based predictions (including random intercepts)
#    re.form = NULL means “include both fixed and random effects”
pred_grid2$fit <- predict(lmm2, newdata = pred_grid2, re.form = NULL)

# 7. Define your custom ggplot2 theme
my_theme <- theme(
  axis.text.x           = element_text(angle = 0, hjust = 0.5),
  plot.background       = element_rect(fill = "white", color = "white"),
  panel.background      = element_rect(fill = "white"),
  legend.background     = element_rect(fill = "white", color = "white"),
  plot.title            = element_text(hjust = 0.5, size = 18, face = "bold"),
  plot.caption          = element_text(face = "bold", size = 16, hjust = 0),
  plot.caption.position = "plot",
  axis.title            = element_text(face = "bold", size = 16),
  axis.text             = element_text(color = "black", size = 14),
  panel.border          = element_blank(),
  panel.grid.major      = element_blank(),
  panel.grid.minor      = element_blank(),
  legend.position       = "none",
  legend.text           = element_text(size = 14)
)

# 8. Plot raw data and fitted lines
ggplot(df, aes(x = bin_median, y = avg_value, color = species)) +
  geom_point(aes(size = count), alpha = 0.6) +
  geom_line(
    data = pred_grid2,
    aes(x = bin_median, y = fit, color = species),
    size = 1
  ) +
  scale_color_manual(values = cb_palette) +
  scale_size_area(max_size = 4) +
  labs(
    x     = "Bin-median PSI",
    y     = "Average NDVI",
    title = "Linear Mixed Model: NDVI vs log(–PSI + 8) by Species"
  ) +
  my_theme


# 1. Prepare data ---------------------------------------------------------

# bin & NA‐cleanup (you already have this)
NDVI_PSIbin_df <- NDVI_PSIbin(df)
NDVI_PSIbin_df <- na.omit(NDVI_PSIbin_df)

# enforce species order
species_order <- c("Oak","Beech","Spruce","Pine")
NDVI_PSIbin_df$species <- factor(NDVI_PSIbin_df$species, levels = species_order)

# compute shifted predictor and its log
NDVI_PSIbin_df$Xpos <- -NDVI_PSIbin_df$bin_median + 8
if(any(NDVI_PSIbin_df$Xpos <= 0)) stop("Some Xpos ≤ 0; increase the +8 shift.")
NDVI_PSIbin_df$logX <- log(NDVI_PSIbin_df$Xpos)


# 2. Assign species colours ------------------------------------------------

cols <- rainbow(length(species_order))
names(cols) <- species_order


# 3. Fit the linear mixed‐effects model -------------------------------------

library(lme4)
lmm2 <- lmer(
  avg_value ~ logX + (1 | species),
  weights = count,
  data    = NDVI_PSIbin_df
)
summary(lmm2)


# 4. Base‐R scatter plot ----------------------------------------------------

with(
  NDVI_PSIbin_df,
  plot(
    Xpos, avg_value,
    col   = cols[species],
    pch   = 19,
    xlab  = expression(-bin_median + 8),
    ylab  = "Average Value",
    main  = expression("LMM fit: " * avg_value ~ "~" ~ log(-bin_median+8))
  )
)
legend(
  "topleft",
  legend = species_order,
  col    = cols,
  pch    = 19,
  title  = "Species"
)


# 5. Overlay fitted lines ---------------------------------------------------

# generate a grid of original bin_medians
bin_seq <- seq(
  min(NDVI_PSIbin_df$bin_median),
  max(NDVI_PSIbin_df$bin_median),
  length = 200
)
# transform to Xpos/logX
grid <- data.frame(
  bin_median = bin_seq,
  Xpos       = -bin_seq + 8,
  logX       = log(-bin_seq + 8)
)

# one line per species
for(sp in species_order) {
  newdat <- transform(grid, species = sp, count = NA)
  # predict includes both fixed + random intercepts by default
  pred <- predict(lmm2, newdata = newdat, re.form = NULL)
  lines(newdat$Xpos, pred, col = cols[sp], lwd = 2)
}


# 1. Libraries --------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(mgcv)    # needed for gam()

# 2. Species order & palette ------------------------------------------------
species_order <- c("Oak", "Beech", "Spruce", "Pine")
cb_palette <- c(
  "Oak"    = "#E69F00",
  "Beech"  = "#0072B2",
  "Spruce" = "#009E73",
  "Pine"   = "#F0E442"
)

# 3. Prepare your df ---------------------------------------------------------
df <- NDVI_PSIbin_df %>%
  # ensure species in the right order
  mutate(species = factor(species, levels = species_order)) %>%
  # (you said you already have these; just here for completeness)
  mutate(
    Xpos = bin_median,
    logX = log(Xpos)
  )

# 4. Custom ggplot2 theme ----------------------------------------------------
my_theme <- theme(
  axis.text.x           = element_text(angle = 0, hjust = 0.5),
  plot.background       = element_rect(fill = "white", color = "white"),
  panel.background      = element_rect(fill = "white"),
  legend.background     = element_rect(fill = "white", color = "white"),
  plot.title            = element_text(hjust = 0.5, size = 18, face = "bold"),
  plot.caption.position = "plot",
  axis.title            = element_text(face = "bold", size = 16),
  axis.text             = element_text(color = "black", size = 14),
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
  panel.grid.major      = element_blank(),
  panel.grid.minor      = element_blank(),
  legend.position       = "top",
  legend.text           = element_text(size = 14)
)

# 5. Plot: raw points + per-species GAM fits -------------------------------
ggplot(df, aes(x = logX, y = avg_value, color = species)) +
  # raw binned averages, sized by count
  geom_point(aes(size = count), alpha = 0.6) +
  # separate GAM per species, weighted by count
  geom_smooth(
    aes(weight = count),
    method  = "gam",
    formula = y ~ s(x, k = 10),
    se      = FALSE
  ) +
  # apply your colours and sizing
  scale_color_manual(values = cb_palette) +
  scale_size_area(max_size = 4) +
  # labels
  labs(
    x     = expression(log( - bin_median + 8 )),
    y     = "Average NDVI",
    title = "GAM fits of NDVI vs log(–PSI + 8) by Species"
  ) +
  # your theme
  my_theme


library(dplyr)
library(ggplot2)

# 1. Specify species order & palette
species_order <- c("Oak", "Beech", "Spruce", "Pine")
cb_palette <- c(
  "Oak"    = "#E69F00",
  "Beech"  = "#0072B2",
  "Spruce" = "#009E73",
  "Pine"   = "#F0E442"
)

# 2. Your custom theme
my_theme <- theme(
  axis.text.x           = element_text(angle = 0, hjust = 0.5),
  plot.background       = element_rect(fill = "white", color = "white"),
  panel.background      = element_rect(fill = "white"),
  legend.background     = element_rect(fill = "white", color = "white"),
  plot.title            = element_text(hjust = 0.5, size = 18, face = "bold"),
  plot.caption.position = "plot",
  axis.title            = element_text(face = "bold", size = 16),
  axis.text             = element_text(color = "black", size = 14),
  panel.border          = element_rect(color = "black", fill = NA, linewidth = 0.5),
  panel.grid.major      = element_blank(),
  panel.grid.minor      = element_blank(),
  legend.position       = "top",
  legend.text           = element_text(size = 14)
)

# 3. Factor your data so ggplot respects that order
df <- NDVI_PSIbin_df %>%
  mutate(species = factor(species, levels = species_order))

# 4. Plot everything in one figure
ggplot(df, aes(
  x     = bin_median,
  y     = avg_value,
  color = species,
  size  = count
)) +
  # raw points
  geom_point(alpha = 0.6) +
  # separate GAM per species
  geom_smooth(
    aes(weight = count),
    method  = "gam",
    formula = y ~ s(x, k = 10),
    se      = FALSE
  ) +
  # apply your manual colours
  scale_color_manual(values = cb_palette) +
  # size scaling for counts
  scale_size_area(max_size = 4) +
  # axis labels + title
  labs(
    x     = "Bin-median PSI",
    y     = "Average NDVI",
    title = "NDVI vs PSI by Species with GAM fits"
  ) +
  # your custom theme (legend at top by default)
  my_theme



# Load necessary libraries
library(nlme)
library(ggplot2)
library(dplyr)

# Your data (replace with your actual data frame)
# NDVI_PSIbin_df <- read.csv("path_to_your_data.csv")

# Ensure your data is structured as shown:
# species (factor), PSI_bin (factor or character), bin_median (numeric), avg_value (numeric), count (integer)

# Define non-linear function
# ndvi_model <- function(x, a, b, q) {
#   a * b^x + q
# }

ndvi_model <- function(x, a, b, q) {
    a * log(x) + q
  }

nlme_fit <- nlme(
  avg_value ~ a * log(bin_median) + q,
  data = NDVI_PSIbin_df,
  fixed = a + q ~ 1,
  random = a + q ~ 1 | species,
  start = c(a = -1, q = 12),
  control = nlmeControl(opt = "nlminb", maxIter = 1000, msMaxIter = 500, pnlsMaxIter = 50)
)

# Check model summary
summary(nlme_fit)


# Prepare predictions for plotting
pred_df <- NDVI_PSIbin_df %>%
  group_by(species) %>%
  summarize(bin_median = seq(min(bin_median), max(bin_median), length.out = 100)) %>%
  ungroup()

# Add predicted values
pred_df$predicted <- predict(nlme_fit, newdata = pred_df)

# Plot original data and model predictions
ggplot() +
  geom_point(data = NDVI_PSIbin_df,
             aes(x = bin_median, y = avg_value, color = species), alpha = 0.6) +
  geom_line(data = pred_df,
            aes(x = bin_median, y = predicted, color = species), size = 1.2) +
  theme_minimal() +
  labs(x = "Bin Median (X)",
       y = "Avg Value (Y)",
       title = "Mixed-Effects Non-linear Model: NDVI ~ PSI Bin Median",
       color = "Species")


library(nlme)
library(ggplot2)
library(dplyr)

species_order <- c("Oak", "Beech", "Spruce", "Pine")
cb_palette <- c(
  "Oak"    = "#E69F00",
  "Beech"  = "#0072B2",
  "Spruce" = "#009E73",
  "Pine"   = "#F0E442"
)

# 2. Your custom theme
my_theme <- theme(
  axis.text.x           = element_text(angle = 0, hjust = 0.5),
  plot.background       = element_rect(fill = "white", color = "white"),
  panel.background      = element_rect(fill = "white"),
  legend.background     = element_rect(fill = "white", color = "white"),
  plot.title            = element_text(hjust = 0.5, size = 18, face = "bold"),
  plot.caption.position = "plot",
  axis.title            = element_text(face = "bold", size = 16),
  axis.text             = element_text(color = "black", size = 14),
  panel.border          = element_rect(color = "black", fill = NA, linewidth = 0.5),
  panel.grid.major      = element_blank(),
  panel.grid.minor      = element_blank(),
  legend.position       = "top",
  legend.text           = element_text(size = 14)
)


# Ensure no bin_median <= 0 (log undefined for 0 or negative)
NDVI_PSIbin_df <- NDVI_PSIbin_df %>% filter(bin_median > 0)

# Define the nonlinear model explicitly
ndvi_model <- function(x, a, b) {
  a * log(x) + b
}

# Nonlinear mixed-effects model using 'nlme' correctly
nlme_fit <- nlme(
  avg_value ~ ndvi_model(-bin_median+8, a, b),
  data = NDVI_PSIbin_df,
  fixed = a + b ~ 1,
  random = a + b ~ 1 | species,
  start = c(a = -1, b = 12),
  control = nlmeControl(opt = "nlminb", maxIter = 1000, msMaxIter = 500, pnlsMaxIter = 50)
)

# Summary of model results
summary(nlme_fit)

# Prepare predictions for plotting
pred_df <- NDVI_PSIbin_df %>%
  group_by(species) %>%
  summarize(bin_median = seq(min(bin_median), max(bin_median), length.out = 100)) %>%
  ungroup() %>%
  mutate(predicted = predict(nlme_fit, newdata = .))

# Plot data and fitted curves
ggplot() +
  geom_point(data = NDVI_PSIbin_df,
             aes(x = bin_median, y = avg_value, color = species), alpha = 0.6) +
  geom_line(data = pred_df,
            aes(x = bin_median, y = predicted, color = species), linewidth = 1.2) +
  theme_minimal() +
  labs(x = "Bin Median (X)",
       y = "Avg Value (Y)",
       title = "Corrected Nonlinear Mixed-effects Model Fit",
       color = "Species")
