###### MAIN FIGURE 2

library(tidyverse)
library(cowplot)
library(fitdistrplus)
library(zoo)
library(ggprism)
library(data.table)
library(scales)
  
source("src/util.R")

##read processed data
rtk <- read_csv(file = "data/hk/hk-Rt-kt.csv")
epicurve <- read_csv(file = "data/hk/hk-epicurves.csv")
rt_npi_sse <- read_csv(file = "data/hk/hk-Rt-npi-sse.csv")
fits <- read_csv(file = "data/hk/hk-waves-R-k-bootstrap.csv")

COVID <- read_rds(file = "data/hk/rds/covid_dated_offspring.rds")
SARS <- read_rds(file = "data/hk/rds/sars_dated_offspring.rds")

offspring <- list(COVID[[1]], COVID[[2]], SARS[[1]], SARS[[2]])
offspring_names <- c(names(COVID), names(SARS))
names(offspring) <- offspring_names


pal <- RColorBrewer::brewer.pal(12, "Paired") # palette for plots
pal_epi <- RColorBrewer::brewer.pal(9, "BuGn") # palette for plots


## Fig.2a - 2c
pA <- fits |> (function(x) {
  x |> 
    mutate(period = case_when(period == "COVID Wave 1" ~ "C1",
                              period == "COVID Wave 2" ~ "C2",
                              period == "COVID Wave 3" ~ "C3",
                              TRUE ~ "S")) |> 
    group_by(sensitivity, period) |>
    reframe(mu_median = quantile(mu, 0.5),
            mu_low = quantile(mu, 0.025),
            mu_high = quantile(mu, 0.975),
            size_median = quantile(size, 0.5),
            size_low = quantile(size, 0.05),
            size_high = quantile(size, 0.95)) |> 
    mutate(sensitivity = case_when(str_detect(sensitivity, pattern = "primary") ~ "Sensitivity 1",
                                   TRUE ~ "Sensitivity 2")) |> 
    ggplot() +
    geom_hline(yintercept = 1, linetype = 2, alpha = 0.2) +
    geom_point(aes(x = period, y = mu_median, colour = sensitivity, fill = sensitivity), shape = 21, position = position_dodge(0.25)) +
    geom_linerange(aes(x = period, ymax = mu_high, ymin = mu_low, colour = sensitivity), position = position_dodge(0.25)) +
    # scale_y_continuous(breaks = c(0,1,2)) +
    theme_classic() +
    theme(
      plot.background = element_blank(),
      strip.background = element_blank(),
      strip.text = element_blank(),
      legend.position = 'bottom',
      legend.text.align = 0) +
    scale_fill_manual(labels = c("Primary", "Sensitivity"), values = c(pal[6], pal[5])) +
    scale_color_manual(labels = c("Primary", "Sensitivity"), values = c(pal[6], pal[5])) +
    labs(x = "",
         y = expression(paste(italic(R))),
         fill = "",
         colour = "") +
    coord_cartesian(ylim = c(0, 2)) +
    guides(y = "prism_offset")
  
}) ()
pB <- fits |> (function(x) {
  x |> 
    mutate(period = case_when(period == "COVID Wave 1" ~ "C1",
                              period == "COVID Wave 2" ~ "C2",
                              period == "COVID Wave 3" ~ "C3",
                              TRUE ~ "S")) |> 
    group_by(sensitivity, period) |>
    reframe(mu_median = quantile(mu, 0.5),
            mu_low = quantile(mu, 0.025),
            mu_high = quantile(mu, 0.975),
            size_median = quantile(size, 0.5),
            size_low = quantile(size, 0.05),
            size_high = quantile(size, 0.95)) |> 
    mutate(sensitivity = case_when(str_detect(sensitivity, pattern = "primary") ~ "Sensitivity 1",
                                   TRUE ~ "Sensitivity 2")) |> 
    ggplot() +
    geom_hline(yintercept = 1, linetype = 2, alpha = 0.2) +
    geom_point(aes(x = period, y = size_median, colour = sensitivity, fill = sensitivity), shape = 21, position = position_dodge(0.25)) +
    geom_linerange(aes(x = period, ymax = size_high, ymin = size_low, colour = sensitivity), position = position_dodge(0.25)) +
    scale_y_log10(breaks = c(0.001, 0.01, 0.1, 1.0, 10), labels = c(0.001, 0.01, 0.1, 1.0, 10)) +
    # scale_x_discrete(labels = c("COVIDWave 1", "COVIDWave 2", "COVIDWave 3", "SARS\nWave")) +
    theme_classic() +
    theme(
      plot.background = element_blank(),
      strip.background = element_blank(),
      strip.text = element_blank(),
      legend.position = 'bottom',
      legend.text.align = 0) +
    scale_fill_manual(labels = c("Primary", "Sensitivity"), values = c(pal[2], pal[1])) +
    scale_color_manual(labels = c("Primary", "Sensitivity"), values = c(pal[2], pal[1])) +
    labs(x = "",
         y = expression(paste(italic(k))),
         fill = "",
         colour = "") +
    coord_cartesian(ylim = c(0.001, 10)) +
    guides(y = "prism_offset")
  
}) ()
pC <- rtk |> (function(x) {
  
 
  x |> 
    ggplot() +
    geom_hline(yintercept = 1, linetype = 2, alpha = 0.2) +
    geom_violin(aes(x = period, y = kt, colour = sensitivity, fill = sensitivity), scale = "width", width = 0.5, alpha = 0.25, trim = FALSE) +
    geom_boxplot(aes(x = period, y = kt, colour = sensitivity, fill = sensitivity), outlier.shape = NA, alpha = 0.7, width = 0.15, position = position_dodge(0.5)) +
    scale_y_log10(breaks = c(0.001, 0.01, 0.1, 1.0, 10), labels = c(0.001, 0.01, 0.1, 1.0, 10)) +
    theme_classic() +
    theme(
      plot.background = element_blank(),
      strip.background = element_blank(),
      strip.text = element_blank(),
      legend.position = 'bottom',
      legend.text.align = 0) +
    scale_fill_manual(labels = c("Primary", "Sensitivity"), values = c(pal[2], pal[1])) +
    scale_color_manual(labels = c("Primary", "Sensitivity"), values = c(pal[2], pal[1])) +
    labs(x = "",
         y = expression(paste(italic(k[t]))),
         fill = "Dataset",
         colour = "Dataset") +
    coord_cartesian(ylim = c(0.001, 10)) +
    guides(y = "prism_offset")
  
  
}) ()

## Fig.2d
pD <- rt_npi_sse |> (function(x){
  
  xl <- x |> 
    group_split(period, sensitivity)
  
  xy <- map(xl, function(y) {
    
    avg_control <- mean(y$ce)
    
    y |> 
      mutate(control = case_when(ce < avg_control ~ "bc",
                                 ce > avg_control ~ "ac"))
    
  }
  ) |> 
    rbindlist() |> 
    mutate(path = case_when(str_detect(string = period, pattern = "COVID") ~ "COVID",
                            TRUE ~ "SARS"))
  
  
  pNEW <- xy |> 
    group_by(period, sensitivity, control) |> 
    reframe(Kt = median(kt),
            Kt_lower = quantile(kt, probs = 0.05),  
            Kt_upper = quantile(kt, probs = 0.95),          
            Rt = mean(rt),
            Rt_lower = quantile(rt, probs = 0.025),          
            Rt_upper = quantile(rt, probs = 0.975),          
    ) |> 
    pivot_wider(names_from = "control",
                values_from = c("Kt", "Kt_lower", "Kt_upper", "Rt", "Rt_lower", "Rt_upper")) |> 
    ggplot() +
    geom_hline(yintercept = 1, linetype = 2, alpha = 0.1) +
    geom_vline(xintercept = 1, linetype = 2, alpha = 0.1) +
    geom_linerange(aes(x = Kt_bc, y = Rt_bc, ymax = Rt_upper_bc, ymin = Rt_lower_bc, colour = sensitivity),
                   linewidth = 0.5, 
                   alpha = 0.3) +
    geom_linerange(aes(x = Kt_bc, y = Rt_bc, xmax = Kt_upper_bc, xmin = Kt_lower_bc, colour = sensitivity),  
                   linewidth = 0.5, 
                   alpha = 0.3) +
    geom_linerange(aes(x = Kt_ac, y = Rt_ac, ymax = Rt_upper_ac, ymin = Rt_lower_ac, colour = sensitivity),
                   linewidth = 0.5, 
                   alpha = 0.3) +
    geom_linerange(aes(x = Kt_ac, y = Rt_ac, xmax = Kt_upper_ac, xmin = Kt_lower_ac, colour = sensitivity),
                   linewidth = 0.5, 
                   alpha = 0.3) +
    ggarchery::geom_arrowsegment(aes(
      x = Kt_bc, 
      y = Rt_bc, 
      xend = Kt_ac,
      yend = Rt_ac,
      colour = sensitivity, 
      fill = sensitivity),
      arrows = arrow(type = 'closed', length = unit(0.1, "inches"))) +
    geom_point(aes(x = Kt_bc, y = Rt_bc, colour = sensitivity, fill = sensitivity), size = 2, shape = 21) +
    theme_classic() +
    scale_x_log10() +
    scale_y_log10() +
    facet_wrap(~period, ncol = 4, scales = "free") +
    theme(
      plot.background = element_blank(),
      strip.background = element_blank(),
      legend.position = 'none',
      legend.text.align = 0) +
    scale_fill_manual(labels = c("Primary", "Sensitivity"), values = c(pal[10], pal[9])) +
    scale_color_manual(labels = c("Primary", "Sensitivity"), values = c(pal[10], pal[9])) +
    labs(x = expression(paste(italic(k[t]))),
         y = expression(paste(italic(R[t]))),
         fill = "",
         colour = "") +
    guides(colour = 'none')
  
  
  ### custom legend, lines for primary sensivitiy, circle or diamond more or less
  
  pLegend <- crossing(
    tibble(sensitivity = c("Primary", "Sensitivity")),
    tibble(control = c("More control", "Less control"))) |> 
    ggplot() +
    geom_line(aes(x = NA, y = NA, colour = sensitivity)) +
    geom_point(aes(x = NA, y = NA, shape = control)) +
    scale_color_manual(labels = c("Primary", "Sensitivity"), values = c(pal[10], pal[9])) +
    scale_shape_manual(values = c(16, 17)) +
    theme_void() +
    theme(
      plot.background = element_blank(),
      aspect.ratio = 0.1, 
      legend.position = 'bottom') +
    labs(colour = "",
         shape = "")
  
  pLegend <- get_legend(pLegend)
  
  plot_grid(pNEW, pLegend, ncol = 1, rel_heights = c(1,0.05))
  
  
}) ()

### Fig.2e 
pE <- epicurve |> (function(x) {
  
  x |> 
    group_by(t, sse_case_classification, period) |> 
    count() |> 
    filter(!str_detect(sse_case_classification, pattern = "Missing")) |> 
    pivot_wider(names_from = sse_case_classification,
                values_from = n) |> 
    janitor::clean_names() |> 
    mutate(
      across(everything(), ~replace_na(.x, 0))
    ) |> 
    mutate(total = primary_sse_exposure_cases + non_sse_cases + other_sse_associated_cases) |> 
    mutate(sse_prop = primary_sse_exposure_cases / total,
           o_sse_prop = other_sse_associated_cases / total,
           non_sse_prop = non_sse_cases / total) |> 
    dplyr::select(- c(total, primary_sse_exposure_cases, non_sse_cases, other_sse_associated_cases)) |> 
    ungroup() |> 
    pivot_longer(cols = c("sse_prop", "o_sse_prop", "non_sse_prop"),
                 names_to = "case_classification",
                 values_to = "proportion") |> 
    group_by(period, case_classification) |> 
    mutate(
      proportion = rollapply(proportion, 
                             width = 7, by = 1, FUN = mean, align = "left", fill = "extend")) |> 
    ggplot(aes(x = t, y = proportion, colour = case_classification)) +
    geom_line() +
    scale_x_date(
      date_breaks = "1 month",
      date_labels = "%b"
    ) +
    scale_colour_manual(
      labels = c("Non-SSE", "Other SSE", "Primary SSE"),
      values = c(pal_epi[[3]], pal_epi[[5]], pal_epi[[7]])) +
    theme_classic() +
    theme(
      plot.background = element_blank(),
      legend.position = "bottom",
      strip.background = element_blank(),
      strip.placement = "outside",
      legend.text = element_text(hjust = 0),
      panel.spacing = unit(1.15, "lines")
    ) +
    facet_wrap(~period, ncol = 1, scales = "free") +
    labs(
      y = expression(paste("Proportion of cases, ", italic(n))),
      x = expression(paste(italic(t))),
      colour = "") +
    coord_cartesian(ylim = c(0,1)) +
    guides(y = "prism_offset") 
    
}) () 

### Fig.2f
pF <- rt_npi_sse |>  (function(x) {
  
  gpal <- ghibli::ghibli_palettes$PonyoMedium
  
  x |>
    ggplot(aes(y = sse_prop, x = kt , colour = sensitivity, fill = sensitivity)) +
    geom_smooth(method = "lm", alpha = 0.1) + 
    stat_binscatter(bins = 20, geom = "pointrange", alpha = 0.3) +
    geom_hline(yintercept = 0.5, linetype = 2, alpha = 0.25) +
    geom_vline(xintercept = 1, linetype = 2, alpha = 0.25) +
    theme_classic() +
    theme(
      plot.background = element_blank(),
      strip.background = element_blank(),
      legend.position = 'bottom',
      legend.text.align = 0,
      plot.title = element_text(hjust = 0.5, size = 10, face = 'italic')
    ) +
    scale_fill_manual(labels = c("Primary", "Sensitivity"), values = c(gpal[2], gpal[3])) +
    scale_color_manual(labels = c("Primary", "Sensitivity"), values = c(gpal[2], gpal[3])) +
    labs(y = expression(paste("Primary superspreading proportion, ", italic(SSE[t]))),
         x = expression(paste("Dispersion, ", italic(k[t]))),
         fill = "",
         colour = "") +
    coord_cartesian(expand = TRUE, clip = "on", xlim = c(0.01, 10), ylim = c(0,1)) +
    scale_x_log10() +
    facet_wrap(~path, scales = "free", ncol = 1) +
    guides(y = "prism_offset",
           x = "prism_offset")
  
}) ()

### Fig.2g
pG <- rt_npi_sse |>  (function(x) {
  
  
  pal_pu <- RColorBrewer::brewer.pal(9, "PuRd") # palette for plots
  
  gpal <- ghibli::ghibli_palettes$PonyoMedium
  
  x |>
    ggplot(aes(y = ce, x = kt , colour = sensitivity, fill = sensitivity)) +
    geom_smooth(method = "lm", alpha = 0.1) + 
    stat_binscatter(bins = 20, geom = "pointrange", alpha = 0.3) +
    geom_hline(yintercept = 0.5, linetype = 2, alpha = 0.25) +
    geom_vline(xintercept = 1, linetype = 2, alpha = 0.25) +
    theme_classic() +
    theme(
      plot.background = element_blank(),
      strip.background = element_blank(),
      legend.position = 'bottom',
      legend.text.align = 0,
      plot.title = element_text(hjust = 0.5, size = 10, face = 'italic')
    ) +
    scale_fill_manual(labels = c("Primary", "Sensitivity"), values = c(gpal[4], gpal[5])) +
    scale_color_manual(labels = c("Primary", "Sensitivity"), values = c(gpal[4], gpal[5])) +
    labs(y = expression(paste("Control effort, ", italic(c[t]))),
         x = expression(paste("Dispersion, ", italic(k[t]))),
         fill = "",
         colour = "") +
    coord_cartesian(expand = TRUE, clip = "on", xlim = c(0.01, 10), ylim = c(0,1)) +
    scale_x_log10() +
    facet_wrap(~path, scales = "free", ncol = 1) +
    guides(y = "prism_offset",
           x = "prism_offset")
  
  
}) ()

### Fig.2h
pH <- rt_npi_sse |>  (function(x) {
  
  gpal <- ghibli::ghibli_palettes$PonyoMedium
  
  x |>
    ggplot(aes(y = growth_rate, x = kt , colour = sensitivity, fill = sensitivity)) +
    geom_smooth(method = "lm", alpha = 0.1) + 
    stat_binscatter(bins = 20, geom = "pointrange", alpha = 0.3) +
    geom_hline(yintercept = 0, linetype = 2, alpha = 0.25) +
    geom_vline(xintercept = 1, linetype = 2, alpha = 0.25) +
    theme_classic() +
    theme(
      plot.background = element_blank(),
      strip.background = element_blank(),
      legend.position = 'bottom',
      legend.text.align = 0,
      plot.title = element_text(hjust = 0.5, size = 10, face = 'italic')
    ) +
    scale_fill_manual(labels = c("Primary", "Sensitivity"), values = c(gpal[6], gpal[7])) +
    scale_color_manual(labels = c("Primary", "Sensitivity"), values = c(gpal[6], gpal[7])) +
    labs(y =expression(paste("Growth rate, ", italic(r[t]))),
         x = expression(paste("Dispersion, ", italic(k[t]))),
         fill = "",
         colour = "") +
    coord_cartesian(expand = TRUE, clip = "on", xlim = c(0.01, 10), ylim = c(-0.2,0.2)) +
    scale_x_log10() +
    facet_wrap(~path, scales = "free", ncol = 1) +
    guides(y = "prism_offset",
           x = "prism_offset")
  
}) ()

## Plot panels
pABC <- plot_grid(pA, pB, pC, ncol = 3, rel_widths = c(0.75, 0.75 ,2.5), labels = c("a", "b", "c"))

pABCD <- plot_grid(pABC, pD, ncol = 1, rel_heights = c(2.5, 3), labels = c("", "d"))

pFGH <- plot_grid(pF, pG, pH, ncol = 3, labels = c("f", "g", "h"))

pE_FGH <- plot_grid(pE, NULL, pFGH, ncol = 3, rel_widths = c(0.8, 0.1, 3), labels = c("e", ""))

p2 <- plot_grid(pABCD, NULL, pE_FGH, ncol = 1, rel_heights = c(2, 0.1, 2.5))

save_plot(plot = p2, filename = "plots/Fig.2.png", base_height = 11, base_width = 11, dpi = 600)

