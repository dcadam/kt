###### MAIN FIGURE 4

library(tidyverse)
library(data.table)
library(ggprism)
library(ggridges)
library(cowplot)

source("src/util.R")


#### read data
mixture_k <- read_csv(file = "data/simulations/sim-mixture-k1.csv")
mixture_mu <- read_csv(file = "data/simulations/sim-mixture-mu.csv")
variance <- read_csv(file = "data/simulations/sim-variance-k1.csv")
relative_offspring <- read_rds(file = "data/simulations/rds/sim-rel-offspring.rds")
zinb <- read_rds(file = "data/simulations/rds/sim-zinb-R-k.rds")


### Figure labels
k_labels <- c(expression(paste(italic(k), " = 0.1")),
              expression(paste(italic(k), " = 0.5")),
              expression(paste(italic(k), " = 1.0")),
              expression(paste(italic(k), " = 10")),
              expression(paste(italic(k), " = 100")))

k_hat_labels <- c(expression(paste(italic(hat(k)), " = 0.08")),
                  expression(paste(italic(hat(k)), " = 0.35")),
                  expression(paste(italic(hat(k)), " = 0.57")),
                  expression(paste(italic(hat(k)), " = 1.51")),
                  expression(paste(italic(hat(k)), " = 1.84")))

pal <- RColorBrewer::brewer.pal(12, "Paired") # palette for plots


#### Fig.4a
pA <- mixture_k |> (function(x) {
  
    ggplot(x) +
    geom_hline(yintercept = 1, linetype = 2, alpha = 0.2) +
    geom_vline(xintercept = 0, linetype = 2, alpha = 0.2) +
    geom_pointrange(aes(x = R_change, y = k, ymax = k_upper, ymin = k_lower, colour = method, fill = method), 
                    position = position_dodge(0.35),
                    shape = 21,
                    alpha = 1) +
    theme_classic() +
    theme(aspect.ratio = 1,
          legend.position = 'bottom') +
    scale_y_log10(breaks = c(0.3, 0.6, 1, 3)) +
    scale_x_continuous(breaks = seq(-2.5, 2.5, by = 1)) +
    annotation_logticks(sides = "l") +
    scale_fill_manual(labels = c("NegBin", "Mixture"), values = c(pal[2], pal[1])) +
    scale_color_manual(labels = c("NegBin", "Mixture"), values = c(pal[2], pal[1])) +
    coord_cartesian(ylim = c(0.3, 3), xlim = c(-2.75, 2.75)) +
    labs(x = expression(paste(italic(R[A] - R[B]))),
         y = expression(paste(italic(hat(k)))),
         colour = "",
         fill = "") 
  
})()

#### Fig.4b
pB <- mixture_mu |> (function(x) {
  
x |> 
      ggplot() +
      geom_abline(slope = 1, linetype = 2, alpha = 0.2) +
      geom_pointrange(aes(x = exp_R, y = mu, ymax = mu_upper, ymin = mu_lower, colour = method, fill = method), 
                      position = position_dodge(0),
                      shape = 21,
                      alpha = 1) +
      theme_classic() +
      theme(aspect.ratio = 1,
            legend.position = 'bottom',
            plot.background = element_blank(),
            strip.background = element_blank(),
            strip.text = element_blank()) +
      scale_y_continuous(breaks = seq(0, 4, by = 0.5)) +
      scale_x_continuous(breaks = seq(0, 4, by = 0.5)) +
      coord_cartesian(ylim = c(0, 3.5), xlim = c(0, 3.5)) +
      labs(x = expression(paste(italic((R[A] + R[B])) / 2)),
           y = expression(paste(italic(hat(R)))),
           colour = "",
           fill = "") +
      scale_fill_manual(labels = c("NegBin", "Mixture"), values = c(pal[6], pal[5])) +
      scale_color_manual(labels = c("NegBin", "Mixture"), values = c(pal[6], pal[5]))
    
  
})()

#### Fig.4c
pC <- variance |> (function(x) {
  
  x |> 
    ggplot() +
    geom_hline(yintercept = 0, linetype = 2, alpha = 0.2) +
    geom_vline(xintercept = 0, linetype = 2, alpha = 0.2) +
    geom_smooth(aes(x = R_change, y = var), se = FALSE, linewidth = 0.5, linetype = 1, colour = "black") +
    geom_smooth(aes(x = R_change, y = var_lower), se = FALSE, linewidth = 0.5, linetype = 3, colour = "black") +
    geom_smooth(aes(x = R_change, y = var_upper), se = FALSE, linewidth = 0.5, linetype = 3, colour = "black") +
    scale_x_continuous(expand = c(0,0), breaks = seq(-2.5, 2.5, by = 1)) +
    scale_y_continuous(breaks = seq(-2,6, by = 1)) +
    theme_classic() +
    theme(aspect.ratio = 1) +
    coord_cartesian(ylim = c(-2, 6)) +
    labs(x = expression(paste(italic(R[A] - R[B]))),
         y = expression(paste(italic(s ^2 - sigma ^2))))
  
}) ()


### Representative fits of underlying and observed densities when R1 = 2 and R2 = 0.5 and k = 1
k <- c(0.1, 0.5, 1.0, 10.0, 100)
single_points <- map(k, function(k) {
  
  R1_points <- data.table(x = 0:10,
                          y = dnbinom(x = 0:10, size = k, mu = 2),
                          source = "Uncontrolled")
  
  R2_points <- data.table(x = 0:10,
                          y = dnbinom(x = 0:10, size = k, mu = 0.5),
                          source = "Controlled")
  
  rbind(R1_points, R2_points) |> 
    mutate(k = k)
  
}) |> 
  rbindlist() |> 
  mutate(k_f = factor(k, labels = k_labels),
         k_h = factor(k, labels = k_hat_labels))

k_hat <- c(0.08, 0.35, 0.57, 1.51, 1.84)
joint_points <- map(k_hat, function(k) {
  
  k_character <- as.character(1.84)
  
  single_points_peak <- single_points[str_detect(single_points$k_h, k_character)] |> 
    group_by(x) |> 
    reframe(y = sum(y))
  
  
  joint_points <- data.table(x = 0:10,
                             y = dnbinom(x = 0:10, size = k, mu = 1)) |> 
    mutate(y_norm = custom_min_max_normalize(y, a = single_points_peak$y[[11]], b = single_points_peak$y[[1]]),
           k = k)
}) |> 
  rbindlist() |> 
  mutate(k_f = factor(k, labels = k_labels),
         k_h = factor(k, labels = k_hat_labels))


#### Fig.4d
pD <- relative_offspring |> (function(x) {
  
  x |> 
    filter(str_detect(k_h, "0.57")) |> 
    mutate(Z_lim = case_when(
      Z >= 10 ~ 10,
      TRUE ~ Z
    )) |>
    ggplot() +
    geom_bar(aes(x = Z_lim, colour = source, fill = source, y = after_stat(prop)),
             size = 0.75, 
             width = 0.75,
             alpha = 0.5
    ) +
    geom_line(stat = 'smooth',
              data = joint_points |> 
                filter(str_detect(k_h, "0.57")),
              mapping = aes(x = x, y = y_norm),
              
              method = "lm", 
              formula = y ~ poly(x, 10),
              se = FALSE,
              linetype = 2, 
              linewidth = 0.75, 
              colour = "#1E6CA8",
              alpha = 0.9
    ) + 
    geom_point(
      data = joint_points |> 
        filter(str_detect(k_h, "0.57")),
      mapping = aes(
        x = x,
        y = y_norm
      ),
      fill = "#1E6CA8",
      shape = 21, 
      alpha = 0.8,
      size = 2.5
    ) +
    facet_wrap(~k_h, ncol = 1, scales = "free", labeller = label_parsed) +
    scale_x_continuous(limits = c(-0.5, 10.5), breaks = c(0:10), labels = c(as.character(0:9), "10+")) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_classic() +
    theme(
      aspect.ratio = 1,
      legend.position = "bottom",
      strip.background = element_blank(),
      strip.placement = "outside",
      legend.text = element_text(hjust = 0),
      panel.spacing = unit(1.15, "lines")
    ) +
    scale_fill_manual(values = c("#1E6CA8", "#1E6CA8"),
                      breaks = c("Uncontrolled", "Controlled")) +
    scale_color_manual(values = c("#1E6CA8", "#1E6CA8"),
                       breaks = c("Uncontrolled", "Controlled")) +
  labs(
    fill = "",
    colour = "",
    shape = "",
    x = expression(paste("Number of secondary cases (", italic(Z), ")")),
    y = "Observed density"
  ) +
    guides(y = "prism_offset") 
  
})()

#### Fig.4e
pE <- relative_offspring |> (function(x) {
  
  pal <- 7
  
  
  x |> 
    filter(str_detect(k_h, "0.57")) |> 
    mutate(Z_lim = case_when(
      Z >= 10 ~ 10,
      TRUE ~ Z
    )) |>
    ggplot() +
    geom_bar(aes(x = Z_lim, colour = source, fill = source, y = after_stat(prop)),
             position = position_dodge2(preserve = "single",
                                        padding = 0.2,
                                        width = 1),
             size = 0.75, 
             alpha = 0.75
    ) +
    geom_line(stat = 'smooth',
              data = single_points |> 
                filter(str_detect(k_h, "0.57")),
              mapping = aes(x = x, y = y, colour = source),
              
              method = "lm", 
              formula = y ~ poly(x, 10),
              se = FALSE,
              linetype = 2, 
              linewidth = 0.75, 
              alpha = 0.9, 
    ) +  
    geom_point(
      data = single_points |> 
        filter(str_detect(k_h, "0.57")),
      mapping = aes(
        x = x,
        y = y,
        fill = source
      ),
      shape = 21, 
      alpha = 0.8,
      size = 2.5
    ) +
    facet_wrap(~k_f, ncol = 1, scales = "free", labeller = label_parsed) +
    scale_x_continuous(limits = c(-0.5, 10.5), breaks = c(0:10), labels = c(as.character(0:9), "10+")) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_classic() +
    theme(
      aspect.ratio = 1,
      legend.position = "bottom",
      strip.background = element_blank(),
      strip.placement = "outside",
      legend.text = element_text(hjust = 0),
      panel.spacing = unit(1.15, "lines")
    ) +
    scale_fill_brewer(
      breaks = c("Uncontrolled", "Controlled"),
      type = "qual",
      palette = pal,
      direction = 1
    ) +
    scale_colour_brewer(
      breaks = c("Uncontrolled", "Controlled"),
      type = "qual",
      palette = pal,
      direction = 1
    ) +
    scale_shape_manual(values = c(21, 22)) +
    labs(
      fill = "",
      colour = "",
      shape = "",
      x = expression(paste("Number of secondary cases (", italic(Z), ")")),
      y = "Underlying density"
    ) +
    guides(y = "prism_offset")
  
  
})()


#### Fig.4f
pF <- bind_rows(
  zinb$zinb |> mutate(method = "2 param estimation"),
  zinb$zinb_fixed |> mutate(method = "1 param estimation")
) |> 
  (function(x) {
    
    x |> 
      mutate(
        pstr0_fct = factor(pstr0, labels = c("0.00", "0.05", "0.10", "0.15", "0.20", "0.25", "0.30", "0.35", "0.40", "0.45", "0.50")),
        method_f = as.factor(method),
        k_input_f = factor(k, labels = k_labels)
      ) |>
      filter(k == 1) |> 
      ggplot() +
      geom_vline(aes(xintercept = k), alpha = 0.2, linetype = 2) +
      geom_density_ridges(
        aes(
          y = pstr0_fct,
          x = size,
          fill = mse_size_normal,
          linetype = method_f,
          group = paste(pstr0_fct, method_f)
        ),
        scale = 0.98,
        rel_min_height = 0.01,
        alpha = 0.8
      ) +
      scale_x_log10(breaks = c(0.3, 0.6, 1, 3)) +
      annotation_logticks(sides = "b") +
      scale_y_discrete(guide = "prism_offset") +
      theme_classic() +
      theme(
        aspect.ratio = 1,
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.text.align = 0,
        panel.spacing = unit(1.15, "lines")
      ) +
      scale_fill_distiller(values = c(0, 0.2, 1), palette = "RdYlGn") +
      scale_linetype_manual(
        values = c(2, 1),
        labels = c(
          expression(paste(italic(k), " only")),
          expression(paste(italic(R - k)))
        )
      ) +
      labs(
        x = expression(paste(italic(hat(k)))),
        y = expression(paste("Probability of structural zeros (", italic(phi), ")")),
        fill = expression(paste("Error"))) +
      guides(linetype = 'none') +
      facet_wrap(~k_input_f, ncol = 1, scales = "free", labeller = label_parsed)
    
    
  })()


p4 <- plot_grid(pA, pB, pC, pD, pE, pF, ncol = 3, labels = c("a", "b", "c", "d", "e", "f"), align = "hv", axis = "tblr")

save_plot(plot = p4, filename = "plots/Fig.4.png", base_height = 8, base_width = 10, dpi = 600)
