###### MAIN FIGURE 3

library(tidyverse)
library(cowplot)
library(fitdistrplus)
library(zoo)
library(ggprism)
library(data.table)
library(scales)
library(RColorBrewer)

# Read data
t_data_ss <- read_csv(file = "data/simulations/sim-epicurves.csv")
rt_kt_est <- read_csv(file = "data/simulations/sim-Rt-kt.csv")

### Figure labels
t_labels <- c(expression(paste(italic(t), " = 1")),
              expression(paste(italic(t), " = 7")),
              expression(paste(italic(t), " = 14")))

k_labels <- c(expression(paste(italic(k), " = 0.1")),
              expression(paste(italic(k), " = 0.5")),
              expression(paste(italic(k), " = 1.0")),
              expression(paste(italic(k), " = 10")),
              expression(paste(italic(k), " = 100")))

### Fig.3a for k = 1
pA <- t_data_ss |> (function(x) {
  
  R1_title <- expression(paste(italic(R), " = 2"))
  R2_title <- expression(paste(italic(R), " = 0.5"))
  
  
  x |>
    filter(k_input == 1) |> 
    group_by(sse) |>
    mutate(nn = rollapply(n,
                          width = 7, by = 1, FUN = mean, align = "center", fill = "extend"
    )) |>
    mutate(nn = case_when(
      sse == "ssn" ~ NA_real_,
      TRUE ~ nn
    )) |>
    mutate(nn_plot = coalesce(nn, n)) |>
    ggplot() +
    geom_vline(xintercept = 60, linetype = 2, alpha = 0.3) +
    geom_histogram(aes(
      x = t,
      y = nn_plot,
      fill = sse
    ),
    stat = "identity",
    position = "dodge",
    colour = "grey",
    linewidth = 0.1
    ) +
    scale_fill_brewer(
      type = "seq",
      palette = 2,
      labels = c("Other cases", "Infected by primary SSE cases", "Primary SSE cases (>6 secondary infections)")
    ) +
    theme_classic() +
    scale_x_continuous(limits = c(0, 120), expand = c(0, 0), breaks = seq(0, 120, by = 10)) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, 400, by = 50)) +
    labs(
      fill = "Case classification",
      y = expression(paste("Mean case count, ", italic(n))),
      x = expression(paste(italic(t)))
    ) +
    theme(
      legend.position = 'bottom',
      strip.background = element_blank(),
      strip.placement = "outside",
      plot.title = element_text(hjust = 0.5, size = 10, face = "italic")
    ) +
    coord_cartesian(ylim = c(0,400), xlim = c(0,125), clip = "off") +
    guides(y = "prism_offset")
}) ()

### Fig.3b for k = 1
pB <- rt_kt_est |> (function(x) {
  
  pal <- RColorBrewer::brewer.pal(n = 9, name = "Reds")
  
  R_input <- tibble(t = 0:120,
                    R = c(rep(2, 60), rep(0.5, 61)))
  
  x |> 
    mutate(k_input_f = factor(k_input, labels = k_labels)) |> 
    filter(k_input == 1) |> 
    mutate(window_f = factor(window_t, labels = t_labels)) |> 
    
    ggplot() +
    geom_vline(xintercept = 60, linetype = 2, alpha = 0.2) +
    geom_line(data = R_input,
              mapping = aes(x = t, y = R)) +
    geom_ribbon(aes(x = t, ymax = Rc_upper, ymin = Rc_lower, fill = factor(window_t)), alpha = 0.1) +
    
    geom_line(aes(x = t, y = Rc, colour = factor(window_t))) + 
    facet_wrap(~k_input_f, ncol = 1, scales = "free", labeller = label_parsed) +
    scale_x_continuous(limits = c(0, 120), expand = c(0, 0), breaks = seq(0, 120, by = 20)) +
    theme_classic() +
    scale_color_manual(values = c(pal[4], pal[6], pal[7])) +
    scale_fill_manual(values = c(pal[4], pal[6], pal[7])) +
    labs(fill = "Window size",
         colour = "Window size",
         x = expression(paste(italic(t))),
         y = expression(paste(italic(hat(R[t]))))) +
    theme(panel.spacing = unit(1.5, "lines"),
          legend.position = "bottom",
          strip.background = element_blank(), 
          strip.placement = "outside") +
    coord_cartesian(expand = TRUE, xlim = c(0, 125), ylim = c(0, 5)) +
    facet_wrap(~window_f, ncol = 3, scales = "free", labeller = label_parsed)
  
})() 

### Fig.3c for k = 1
pC <- rt_kt_est |> (function(x) {
  
  x |> 
    mutate(k_input_f = factor(k_input, labels = k_labels)) |> 
    mutate(window_f = factor(window_t, labels = t_labels)) |> 
    filter(k_input == 1) |> 
    ggplot() +
    geom_vline(xintercept = 60, linetype = 2, alpha = 0.2) +
    geom_hline(yintercept = 1) +
    geom_ribbon(aes(x = t, ymax = kc_upper, ymin = kc_lower, fill = factor(window_t)), alpha = 0.1) +
    geom_line(aes(x = t, y = kc, colour = factor(window_t))) + 
    facet_wrap(~k_input_f, ncol = 1, scales = "free", labeller = label_parsed) +
    scale_x_continuous(limits = c(0, 120), expand = c(0, 0), breaks = seq(0, 120, by = 20)) +
    scale_y_log10() +
    theme_classic() +
    scale_color_manual(values = colorRampPalette(brewer.pal(9, "YlGnBu"))(12)[6:12]) +
    scale_fill_manual(values = colorRampPalette(brewer.pal(9, "YlGnBu"))(12)[6:12]) +
    labs(fill = "Window size",
         colour = "Window size",
         x = expression(paste(italic(t))),
         y = expression(paste(italic(hat(k[t]))))) +
    theme(panel.spacing = unit(1.5, "lines"),
          legend.position = "bottom",
          strip.background = element_blank(), 
          strip.placement = "outside") +
    coord_cartesian(expand = TRUE, xlim = c(0, 125), ylim = c(0.1, 100)) +
    facet_wrap(~window_f, ncol = 3, scales = "free", labeller = label_parsed)
  
})()

p3 <- plot_grid(pA, pB, pC, ncol = 1, labels = c("a", "b", "c"), align = 'hv', axis = 'tblr', rel_heights = c(2.5, 2, 2))

save_plot(plot = p3, filename = "plots/Fig.3.png", base_width = 10, base_height = 9, dpi = 600)
save_plot(plot = p3, filename = "plots/Fig.3.pdf", base_width = 10, base_height = 9)



