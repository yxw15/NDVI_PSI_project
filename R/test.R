# -----------------------------
# Setup
# -----------------------------
setwd("/dss/dssfs02/lwp-dss-0001/pr48va/pr48va-dss-0000/yixuan/NDVI_PSI_project")
load("results/Data/AllSpecies_AllMonths_rootzone.RData")

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(tibble)
library(purrr)
library(scales)

# -----------------------------
# Prepare data
# -----------------------------
data_all <- combined %>%
  filter(Quantiles > 0) %>%
  filter(month %in% c("May","June","July","August")) %>%
  na.omit()

data <- data_all %>%
  filter(month %in% c("July","August")) %>%
  filter(year < 2023)

# -----------------------------
# Function: Bin transpiration deficit by soil water potential
# -----------------------------
TDiff_PSIbin <- function(df, bin_width = 50) {
  value_column <- "transpiration_deficit"
  
  species_totals <- df %>%
    group_by(species) %>%
    summarise(total_pixels = n(), .groups = "drop")
  
  psi_min <- floor(min(df$soil_water_potential, na.rm = TRUE))
  psi_max <- ceiling(max(df$soil_water_potential, na.rm = TRUE))
  bin_breaks <- seq(psi_min, psi_max, by = bin_width)
  
  df <- df %>%
    mutate(PSI_bin = cut(soil_water_potential, breaks = bin_breaks,
                         include.lowest = TRUE, right = FALSE))
  
  meanTDiff_PSIbin_species <- df %>%
    group_by(species, PSI_bin) %>%
    summarise(
      avg_transpiration_deficit = mean(.data[[value_column]], na.rm = TRUE),
      count = n(),
      .groups = "drop"
    ) %>%
    mutate(
      bin_median = sapply(as.character(PSI_bin), function(bin_label) {
        nums <- as.numeric(strsplit(gsub("\\[|\\]|\\(|\\)", "", bin_label), ",")[[1]])
        mean(nums)
      })
    ) %>%
    left_join(species_totals, by = "species") %>%
    mutate(percentage = count / total_pixels) %>%
    filter(percentage > 0.001) %>%
    select(species, PSI_bin, bin_median, avg_transpiration_deficit,
           count, total_pixels, percentage)
  
  return(meanTDiff_PSIbin_species)
}

# -----------------------------
# Main plotting function
# -----------------------------
plot_TDiff_PSIbin_slope <- function(data,
                                    coef_output,
                                    figure_output,
                                    aic_barplot_fig) {
  # -----------------------------
  # Prepare data & species info
  # -----------------------------
  df <- TDiff_PSIbin(data) %>% na.omit()
  
  species_levels <- c("Oak","Beech","Spruce","Pine")
  cb_palette <- c("Oak" = "#E69F00", "Beech" = "#0072B2",
                  "Spruce" = "#009E73", "Pine" = "#F0E442")
  
  df$species <- factor(df$species, levels = species_levels)
  xcol <- "bin_median"
  ycol <- "avg_transpiration_deficit"
  
  # -----------------------------
  # Fit polynomial models & choose by AIC
  # -----------------------------
  degrees <- c(1,2,3)
  fit_one_species <- function(sp) {
    sdf <- df %>% filter(species == sp) %>% select(all_of(c(xcol, ycol)))
    
    if (nrow(sdf) < 5 || any(!is.finite(sdf[[xcol]])) || any(!is.finite(sdf[[ycol]]))) {
      return(list(species = sp,
                  aic_tbl = tibble(model = c("linear","poly2","poly3"), AIC = NA_real_),
                  best_deg = NA_integer_,
                  best_fit = NULL))
    }
    
    fits <- lapply(degrees, function(d) {
      fm <- reformulate(termlabels = paste0("poly(", xcol, ", ", d, ")"), response = ycol)
      tryCatch(lm(fm, data = sdf), error = function(e) NULL)
    })
    
    aics <- sapply(fits, function(m) if (is.null(m)) NA_real_ else AIC(m))
    aic_tbl <- tibble(model = c("linear","poly2","poly3"), AIC = aics)
    
    best_idx <- suppressWarnings(which.min(aics))
    best_deg <- if (length(best_idx) == 0 || all(is.na(aics))) NA_integer_ else degrees[best_idx]
    best_fit <- if (is.na(best_deg)) NULL else fits[[best_idx]]
    
    list(species = sp, aic_tbl = aic_tbl, best_deg = best_deg, best_fit = best_fit)
  }
  
  fits_all <- lapply(species_levels, fit_one_species)
  
  # -----------------------------
  # AIC barplot
  # -----------------------------
  aic_df <- bind_rows(lapply(fits_all, function(x) {
    x$aic_tbl %>% mutate(species = x$species)
  }))
  aic_df$species <- factor(aic_df$species, levels = species_levels)
  aic_df$model <- factor(aic_df$model, levels = c("linear","poly2","poly3"))
  
  p_aic <- ggplot(aic_df, aes(x = species, y = AIC, fill = model)) +
    geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.6) +
    scale_fill_manual(values = c("linear"="#0072B2","poly2"="#E69F00","poly3"="#009E73")) +
    labs(x="", y="AIC") +
    theme_minimal() +
    theme(legend.position = "top",
          axis.text = element_text(size=14),
          axis.title = element_text(size=16, face="bold"))
  
  ggsave(aic_barplot_fig, p_aic, width=8, height=6, dpi=300)
  
  # -----------------------------
  # Predictions for plotting
  # -----------------------------
  pred_list <- lapply(fits_all, function(ff) {
    sp <- ff$species
    if (is.null(ff$best_fit) || is.na(ff$best_deg)) return(NULL)
    sp_df <- df %>% filter(species == sp)
    xr <- range(sp_df[[xcol]], na.rm = TRUE)
    if (!is.finite(xr[1]) || !is.finite(xr[2]) || xr[1] == xr[2]) return(NULL)
    xseq <- seq(xr[1], xr[2], length.out = 200)
    nd <- data.frame(bin_median = xseq, species = factor(sp, levels=species_levels))
    nd$predicted <- predict(ff$best_fit, newdata=nd)
    nd$degree <- ff$best_deg
    nd
  })
  pred_df <- bind_rows(pred_list)
  
  # line_median <- median(df[[ycol]], na.rm = TRUE)
  line_mean <- mean(df[[ycol]], na.rm = TRUE)
  
  # -----------------------------
  # Panel (a)
  # -----------------------------
  plot_mixed <- ggplot() +
    geom_point(data = df,
               aes(x = .data[[xcol]], y = .data[[ycol]],
                   color = species, shape = species, size = percentage),
               alpha = 0.7) +
    geom_line(data = pred_df, aes(x = bin_median, y = predicted, color = species), linewidth=1) +
    geom_hline(yintercept = line_mean, linetype="dashed", color="black") +
    annotate("text", x=-2000, y=line_mean, label="mean", vjust=-0.3, fontface="italic", size=6) +
    scale_color_manual(values=cb_palette) +
    scale_shape_manual(values=c("Oak"=16,"Beech"=17,"Spruce"=15,"Pine"=18)) +
    scale_size_continuous(name="Pixels per bin (%)", range=c(1,8), labels=percent_format(accuracy=1)) +
    labs(x="soil water potential (kPa)", y="transpiration deficit") +
    ggtitle("(a)") +
    theme_minimal() +
    theme(legend.position="bottom", legend.text=element_text(size=14))
  
  # -----------------------------
  # Panel (b) - root at mean
  # -----------------------------
  solve_x_at_y <- function(fit, y_target, x_range, n_grid=400) {
    f <- function(x) predict(fit, newdata=data.frame(bin_median=x)) - y_target
    xs <- seq(x_range[1], x_range[2], length.out=n_grid)
    fs <- suppressWarnings(f(xs))
    sc_idx <- which(diff(sign(fs)) != 0 & is.finite(fs[-length(fs)]) & is.finite(fs[-1]))
    if (length(sc_idx)==0) return(NA_real_)
    pick <- sc_idx[which.min(pmin(abs(fs[sc_idx]), abs(fs[sc_idx+1])))]
    lo <- xs[pick]; hi <- xs[pick+1]
    tryCatch(uniroot(f, interval=c(lo, hi))$root, error=function(e) NA_real_)
  }
  
  x_at_mean <- lapply(fits_all, function(ff) {
    sp <- ff$species
    fit <- ff$best_fit
    if (is.null(fit)) return(tibble(species=sp, x_at_mean=NA_real_))
    sp_df <- df %>% filter(species==sp)
    xr <- range(sp_df[[xcol]], na.rm=TRUE)
    xm <- solve_x_at_y(fit, line_mean, xr, n_grid=600)
    tibble(species=sp, x_at_mean=xm)
  }) %>% bind_rows() %>% mutate(species=factor(species, levels=species_levels))
  
  p_bar_x <- ggplot(x_at_mean, aes(x=species, y=x_at_mean, fill=species)) +
    geom_bar(stat="identity", width=0.7, na.rm=TRUE) +
    scale_fill_manual(values=cb_palette, guide="none") +
    labs(x="", y="soil water potential (kPa)") +
    ggtitle("(b)") +
    theme_minimal()
  
  # -----------------------------
  # Panel (c) - slope at median
  # -----------------------------
  get_poly_derivative <- function(fit) {
    function(x0) {
      h <- 1e-6
      (predict(fit, newdata=data.frame(bin_median=x0+h)) -
          predict(fit, newdata=data.frame(bin_median=x0-h))) / (2*h)
    }
  }
  
  local_slope_data <- x_at_mean %>%
    left_join(bind_rows(lapply(fits_all, function(ff) tibble(species=ff$species, fit=list(ff$best_fit), degree=ff$best_deg))),
              by="species") %>%
    group_by(species) %>%
    group_modify(~{
      fit <- .x$fit[[1]]; deg <- .x$degree[1]; x_m <- .x$x_at_mean[1]
      if(is.null(fit) || is.na(x_m)) return(tibble(slope=NA, slope_se=NA, p_value=NA, r2=NA))
      
      # derivative at median
      d_fun <- get_poly_derivative(fit)
      slope_val <- d_fun(x_m)
      
      # local regression for SE and R2
      d <- pred_df %>% filter(species==.y$species)
      k <- which.min(abs(d$bin_median - x_m))
      w0 <- max(1, k-5); w1 <- min(nrow(d), k+5)
      window <- d[w0:w1,,drop=FALSE]
      sm <- summary(lm(predicted ~ bin_median, data=window))
      
      tibble(slope=slope_val,
             slope_se=sm$coefficients["bin_median","Std. Error"],
             p_value=sm$coefficients["bin_median","Pr(>|t|)"],
             r2=sm$r.squared)
    }) %>%
    ungroup() %>%
    mutate(species=factor(x_at_mean$species, levels=species_levels),
           slope_abs=abs(slope),
           label_text=ifelse(p_value<0.05, sprintf("%.2f*", r2), sprintf("%.2f", r2)))
  
  p_bar <- ggplot(local_slope_data, aes(x=species, y=slope_abs, fill=species)) +
    geom_bar(stat="identity", width=0.7) +
    geom_errorbar(aes(ymin=pmax(0,slope_abs-slope_se), ymax=slope_abs+slope_se), width=0.2) +
    geom_text(aes(y=slope_abs/2, label=label_text), size=5, color="black") +
    scale_fill_manual(values=cb_palette, guide="none") +
    labs(x="", y="absolute slope") +
    ggtitle("(c)") +
    theme_minimal()
  
  # -----------------------------
  # Combine panels
  # -----------------------------
  final_plot <- (plot_mixed + (p_bar_x / p_bar)) +
    plot_layout(widths=c(2,1), guides="collect") &
    theme(legend.position="bottom", legend.title=element_blank())
  
  ggsave(figure_output, final_plot, width=10, height=8, dpi=300)
  
  # -----------------------------
  # Coefficient plot
  # -----------------------------
  coef_rows <- lapply(fits_all, function(ff) {
    sp <- ff$species; fit <- ff$best_fit; deg <- ff$best_deg
    if (is.null(fit) || is.na(deg)) return(NULL)
    sm <- summary(fit)$coefficients
    out <- as.data.frame(sm); out$term_raw <- rownames(out)
    out$species <- sp; out$degree <- deg; out
  }) %>% bind_rows()
  
  if (!is.null(coef_rows) && nrow(coef_rows) > 0) {
    coef_data <- coef_rows %>%
      mutate(term = dplyr::case_when(
        term_raw=="(Intercept)" ~ "a",
        grepl("poly\\(bin_median,\\s*1\\)1", term_raw) ~ "b",
        grepl("poly\\(bin_median,\\s*2\\)1", term_raw) ~ "b",
        grepl("poly\\(bin_median,\\s*2\\)2", term_raw) ~ "c",
        grepl("poly\\(bin_median,\\s*3\\)1", term_raw) ~ "b",
        grepl("poly\\(bin_median,\\s*3\\)2", term_raw) ~ "c",
        grepl("poly\\(bin_median,\\s*3\\)3", term_raw) ~ "d",
        TRUE ~ term_raw),
        species = factor(species, levels=species_levels),
        value = Estimate,
        p_value = `Pr(>|t|)`,
        label_text = ifelse(p_value<0.05, "*", sprintf("%.2f", p_value))
      ) %>%
      filter((term %in% c("a","b") & degree >= 1) |
               (term %in% c("c") & degree >= 2) |
               (term %in% c("d") & degree >= 3))
    
    plot_coeff <- ggplot(coef_data, aes(x=species, y=value, fill=species)) +
      geom_bar(stat="identity", position=position_dodge(0.8), width=0.7) +
      geom_text(aes(label=label_text, y=value/2), color="black", size=5, position=position_dodge(0.8)) +
      facet_wrap(~ term, scales="free_y") +
      scale_fill_manual(values=cb_palette) +
      labs(title="Transpiration deficit ~ Soil water potential", x="", y="coefficient") +
      theme_minimal()
    
    ggsave(coef_output, plot_coeff, width=10, height=8, dpi=300)
  }
  
  return(list(combined_plot=final_plot, coeff_plot=plot_coeff, aic_plot=p_aic))
}

# -----------------------------
# Run plotting
# -----------------------------
plot_TDiff_PSIbin_slope(
  data = data,
  coef_output = "results_rootzone/Figures_till2022/supplementary/TDiff_PSIbin_coeff_till2022.png",
  figure_output = "results_rootzone/Figures_till2022/supplementary/TDiff_PSIbin_slope_till2022.png",
  aic_barplot_fig = "results_rootzone/Figures_till2022/supplementary/TDiff_PSIbin_AIC_till2022.png"
)
