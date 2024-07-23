###### MAIN FIGURE 3

library(tidyverse)
library(cowplot)
library(fitdistrplus)
library(zoo)
library(ggprism)
library(data.table)
library(scales)
library(RColorBrewer)

t_data <- read_rds(file = "code/simulations/data/raw_sim_data.rds")  ## [[i = 1:5]][[j = 1:300 (to split) ]]


### Get average epidemic size and characterise by link to SSEs by R & k
ss_cases <- map(t_data, function(x) {
  
  x |> 
    mutate(ss_linked = case_when(obs_n >= 6 ~ 1,
                                 TRUE ~ 0),
           ss_case_id = case_id) |> 
    filter(ss_linked == 1) |> 
    dplyr::select(epidemic_id, ss_case_id, ss_linked)
  
})
t_data_ss <- map2(t_data, ss_cases, function(x, y) {
  
  x |> 
    mutate(ss = case_when(obs_n >= 6 ~ 1,
                          TRUE ~ 0)) |> 
    mutate(ss_case_id = as.numeric(source)) |> 
    left_join(y, by = c("ss_case_id", "epidemic_id")) |>
    mutate(ss_linked = replace_na(ss_linked, 0)) |> 
    mutate(ss_linked = case_when(ss == 1 & ss_linked == 1 ~ 0,
                                 TRUE ~ ss_linked)) |> 
    group_by(epidemic_id, t_inf) |>
    reframe(ssn = sum(ss),
            ss_link = sum(ss_linked),
            nn = n() - ssn - ss_link) |> 
    group_by(t_inf) |> 
    # reframe(ssn = quantile(x = ssn, 0.5),
    #         ss_link = quantile(x = ss_link, 0.5),
    #         nn = quantile(x = nn, 0.5)) |> 
    reframe(ssn = mean(x = ssn),
            ss_link = mean(x = ss_link),
            nn = mean(x = nn)) |>
    pivot_longer(cols = c("ssn", "ss_link", "nn"),
                 names_to = "sse",
                 values_to = "n") |> 
    rename(t = t_inf) |> 
    mutate(k_input = unique(x$k))
})




### Figure labels
t_labels <- c(expression(paste(italic(t), " = 1")),
              expression(paste(italic(t), " = 7")),
              expression(paste(italic(t), " = 14")))


pA <- t_data_ss[[3]] |> (function(x) {
  
  R1_title <- expression(paste(italic(R), " = 2"))
  R2_title <- expression(paste(italic(R), " = 0.5"))
  
  
  x |>
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
    # annotate(geom = "rect", xmin = 60, xmax = 120, ymin = 0, ymax = 400, fill = "grey", alpha = 0.1) +
    # annotate(geom = "text", x = 30, y = 410, label = R1_title) +
    # annotate(geom = "text", x = 90, y = 410, label = R2_title) +
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


pA

## READ data
rc_kc_inf <- read_rds(file = "code/simulations/data/out_rc_kc.rds")

nlim <- 150
rc_kc_est <- map(rc_kc_inf, function(x) {
  
  x |> 
    group_by(t_start, t_end, window, k_input) |> 
    mutate(n = n()) |> 
    filter(n >= 150) |> 
    reframe(
      Rc =  mean(r, na.rm = TRUE),
      Rc_lower = quantile(r, 0.025, na.rm = TRUE),
      Rc_upper = quantile(r, 0.975, na.rm = TRUE),
      kc =  quantile(k, 0.5, na.rm = TRUE),
      kc_lower = quantile(k, 0.05, na.rm = TRUE),
      kc_upper = quantile(k, 0.95, na.rm = TRUE)
    ) |> 
    group_by(window, k_input) |>
    mutate(
      kc = rollapply(kc,
                     width = 3, by = 1, FUN = mean, align = "center", fill = "extend"),
      kc_upper = rollapply(kc_upper,
                           width = 3, by = 1, FUN = mean, align = "center", fill = "extend"),
      kc_lower = rollapply(kc_lower,
                           width = 3, by = 1, FUN = mean, align = "center", fill = "extend"),
      Rc = rollapply(Rc,
                     width = 3, by = 1, FUN = mean, align = "center", fill = "extend"),
      Rc_upper = rollapply(Rc_upper,
                           width = 3, by = 1, FUN = mean, align = "center", fill = "extend"),
      Rc_lower = rollapply(Rc_lower,
                           width = 3, by = 1, FUN = mean, align = "center", fill = "extend")
    ) |>
    ungroup() |>
    mutate(window_t = case_when(str_detect(window, "daily") ~ 1,
                                str_detect(window, "weekly") ~ 7,
                                TRUE ~ 14)) |> 
    mutate(t = case_when(window_t == 1 ~ t_start,
                         TRUE ~ t_start + (window_t/2)))
  
  
})

pB <- rbindlist(rc_kc_est) |> (function(x) {
  
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
    # guides(y = "prism_offset") +
    facet_wrap(~window_f, ncol = 3, scales = "free", labeller = label_parsed)
  
})() ## Rt

pB

pC <- rbindlist(rc_kc_est) |> (function(x) {
  
  x |> 
    mutate(k_input_f = factor(k_input, labels = k_labels)) |> 
    mutate(window_f = factor(window_t, labels = t_labels)) |> 
    filter(k_input == 1) |> 
    ggplot() +
    geom_vline(xintercept = 60, linetype = 2, alpha = 0.2) +
    geom_hline(yintercept = 1) +
    geom_ribbon(aes(x = t, ymax = kc_upper, ymin = kc_lower, fill = factor(window_t)), alpha = 0.1) +
    # geom_smooth(aes(x = t, y = kc, colour = factor(window_t))) +
    geom_line(aes(x = t, y = kc, colour = factor(window_t))) + 
    facet_wrap(~k_input_f, ncol = 1, scales = "free", labeller = label_parsed) +
    scale_x_continuous(limits = c(0, 120), expand = c(0, 0), breaks = seq(0, 120, by = 20)) +
    # annotation_logticks(sides = "l") +
    # scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
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
    # guides(y = "prism_offset") +
    facet_wrap(~window_f, ncol = 3, scales = "free", labeller = label_parsed)
  
})() ## Kt

ABC <- plot_grid(pA, pB, pC, ncol = 1, labels = c("a", "b", "c"), align = 'hv', axis = 'tblr', rel_heights = c(2.5, 2, 2))

save_plot(plot = ABC, filename = "manuscript/Fig.3.pdf", base_width = 10, base_height = 9)



########### SUPPLEMENTARY



sA1 <- rbindlist(t_data_ss) |>
  (function(x) {
  
 
  x |>
    mutate(k_input_f = factor(k_input, labels = k_labels)) |> 
    group_by(sse) |>
    mutate(nn = rollapply(n,
                          width = 7, by = 1, FUN = mean, align = "center", fill = "extend"
    )) |>
    mutate(nn = case_when(
      sse == "ssn" ~ NA_real_,
      TRUE ~ nn
    )) |>
    mutate(nn_plot = coalesce(nn, n)) |>
      filter(k_input != 1) |> 
    ggplot() +
    # annotate(geom = "rect", xmin = 60, xmax = 120, ymin = 0, ymax = 400, fill = "grey", alpha = 0.1) +
    # annotate(geom = "text", x = 30, y = 410, label = R1_title) +
    # annotate(geom = "text", x = 90, y = 410, label = R2_title) +
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
    scale_y_continuous(expand = c(0, 0)) +
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
    guides(y = "prism_offset", 
           fill = guide_legend(ncol = 1)) +
  facet_wrap(~k_input_f, ncol = 1, scales = "free", labeller = label_parsed)
  
}) ()

sA2 <- rbindlist(t_data_ss) |> 
   (function(x) {
  
  
  x |>
    mutate(k_input_f = factor(k_input, labels = k_labels)) |> 
    group_by(sse) |>
    mutate(nn = rollapply(n,
                          width = 7, by = 1, FUN = mean, align = "center", fill = "extend"
    )) |>
    mutate(nn = case_when(
      sse == "ssn" ~ NA_real_,
      TRUE ~ nn
    )) |>
    mutate(nn_plot = coalesce(nn, n)) |>
       filter(k_input != 1) |> 
    ggplot() +
    # annotate(geom = "rect", xmin = 60, xmax = 120, ymin = 0, ymax = 400, fill = "grey", alpha = 0.1) +
    # annotate(geom = "text", x = 30, y = 410, label = R1_title) +
    # annotate(geom = "text", x = 90, y = 410, label = R2_title) +
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
    scale_y_continuous(expand = c(0, 0)) +
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
    guides(y = "prism_offset", 
           fill = guide_legend(ncol = 1)) +
    facet_wrap(~k_input_f, ncol = 1, scales = "free_x", labeller = label_parsed) +
    coord_cartesian(xlim = c(0,120))
  
}) ()

sA12 <- plot_grid(sA1, sA2)

save_plot(plot = sA12, filename = "manuscript/supplementary/sim_epicurves.pdf", base_width = 10, base_height = 8)


sB <- rbindlist(rc_kc_est) |> (function(x) {
  
  pal <- RColorBrewer::brewer.pal(n = 9, name = "Reds")
  
  R_input <- tibble(t = 0:120,
                    R = c(rep(2, 60), rep(0.5, 61)))
  
  x |> 
    mutate(k_input_f = factor(k_input, labels = k_labels)) |> 
    mutate(window_f = factor(window_t, labels = t_labels)) |> 
    filter(k_input != 1) |> 
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
    coord_cartesian(expand = TRUE, xlim = c(0, 125)) +
    # guides(y = "prism_offset") +
    facet_wrap(~k_input_f*window_f, ncol = 3, scales = "free", labeller = label_parsed)
  
})() ## Rt

save_plot(plot = sB, filename = "manuscript/supplementary/sim_Rt_window.pdf", base_width = 10, base_height = 10)


sC1 <- rbindlist(rc_kc_est) |> (function(x) {
  
  x |> 
    mutate(k_input_f = factor(k_input, labels = k_labels)) |> 
    mutate(window_f = factor(window_t, labels = t_labels)) |> 
    filter(k_input != 1) |> 
    ggplot() +
    geom_vline(xintercept = 60, linetype = 2, alpha = 0.2) +
    geom_hline(aes(yintercept = k_input)) +
    geom_ribbon(aes(x = t, ymax = kc_upper, ymin = kc_lower, fill = factor(window_t)), alpha = 0.1) +
    # geom_smooth(aes(x = t, y = kc, colour = factor(window_t))) +
    geom_line(aes(x = t, y = kc, colour = factor(window_t))) + 
    scale_x_continuous(limits = c(0, 120), expand = c(0, 0), breaks = seq(0, 120, by = 20)) +
    # annotation_logticks(sides = "l") +
    # scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
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
    coord_cartesian(expand = TRUE, xlim = c(0, 125)) +
    facet_wrap(~k_input_f*window_f, ncol = 3, scales = "free", labeller = label_parsed)
  
})() ## Kt

save_plot(plot = sC1, filename = "manuscript/supplementary/sim_kt_window_ribbon.pdf", base_width = 10, base_height = 10)


sC2 <- rbindlist(rc_kc_est) |> (function(x) {
  
  x |> 
    mutate(k_input_f = factor(k_input, labels = k_labels)) |> 
    mutate(window_f = factor(window_t, labels = t_labels)) |> 
    ggplot() +
    geom_vline(xintercept = 60, linetype = 2, alpha = 0.2) +
    geom_hline(aes(yintercept = k_input)) +
    # geom_ribbon(aes(x = t, ymax = kc_upper, ymin = kc_lower, fill = factor(window_t)), alpha = 0.1) +
    # geom_smooth(aes(x = t, y = kc, colour = factor(window_t))) +
    geom_line(aes(x = t, y = kc, colour = factor(window_t))) + 
    scale_x_continuous(limits = c(0, 120), expand = c(0, 0), breaks = seq(0, 120, by = 20)) +
    # annotation_logticks(sides = "l") +
    # scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
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
    coord_cartesian(expand = TRUE, xlim = c(0, 125)) +
    facet_wrap(~k_input_f*window_f, ncol = 3, scales = "free", labeller = label_parsed)
  
})() ## Kt

save_plot(plot = sC2, filename = "manuscript/supplementary/sim_kt_window.pdf", base_width = 12, base_height = 12)

ABC <- plot_grid(pA, pB, pC, ncol = 1, labels = c("a", "b", "c"), align = 'hv', axis = 'tblr', rel_heights = c(2.5, 2, 2))

















### import figure data
data_1a <- read_rds(file = "code/simulations/finished scripts/data/nb_sims_error.rds")
data_1b <- read_rds(file = "code/simulations/data/data_1b.rds")
data_1c <- read_rds(file = "code/simulations/data/data_1c.rds")
data_1d <- read_rds(file = "code/simulations/data/data_1d.rds")

t_data <- read_rds(file = "code/simulations/data/raw_sim_data.rds")  ## [[i = 1:5]][[j = 1:300 (to split) ]]
sim_curves <- read_rds(file = "code/simulations/finished scripts/data/epicurve_sims_t.rds")


data_1a <- bind_cols(data_1a, 
                     bind_rows(
                       tibble(
                         k_h = factor(k, labels = k_hat_labels)),
                       tibble(
                         k_h = factor(k, labels = k_hat_labels))))

k <- c(0.1, 0.5, 1, 10, 100)


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




### 3a
pA <- data_1a |> (function(x) {
  
  x |> 
  filter(str_detect(k_h, "0.57")) |> 
    mutate(control = factor(control, levels = c("Controlled", "Uncontrolled"), labels = c("Observed", "Expected"))) |> 
    unnest(dist) |> 
    mutate(dist_change = case_when(dist >= 10 ~ 10,
                                   TRUE ~ dist)) |> 
    ggplot(aes(x = dist_change, colour  = control, fill = control, y = after_stat(prop))) +
    geom_bar(position = position_dodge2(padding = 0.2), width = 0.9, alpha = 0.75) +
    facet_wrap(~k_f, ncol = 1, scales = "free", labeller = label_parsed) +
    scale_x_continuous(limits = c(-0.5,10.5), breaks = c(0:10), labels = c(as.character(0:9), "10+")) +
    scale_y_continuous(limits = c(0,0.6)) +
    theme_classic() +
    theme(aspect.ratio = 1,
          legend.position = "bottom",
          strip.background = element_blank(),
          strip.placement = "outside",
          legend.text = element_text(hjust = 0),
          panel.spacing = unit(1.15, "lines")
    ) +
    scale_fill_brewer(breaks = c("Expected", "Observed"),
                      labels = c("Expected", expression(paste("Observed (", italic(k), " = 0.57)"))),
                      type = 'qual', 
                      palette = 3,
                      direction = -1) +
    scale_colour_brewer(breaks = c("Expected", "Observed"),
                        labels = c("Expected", expression(paste("Observed (", italic(k), " = 0.57)"))),
                        type = 'qual', 
                        palette = 3,
                        direction = -1) +
    labs(fill = "",
         colour = "",
         x = expression(paste("Number of secondary cases (", italic(Z), ")")), 
         y = "Proportion of secondary cases") +
    guides(y = "prism_offset")
  
}) ()
 
### 3b
pB <- list_rbind(data_1b) |> (function(x) {
    
    x |> 
    filter(str_detect(k_h, "0.57")) |> 
    mutate(control = factor(control, levels = c("Controlled", "Uncontrolled"), labels = c("Observed", "Expected"))) |>
    ggplot(aes(x = x, y = gy, colour = factor(control))) +
    geom_vline(xintercept = 5.4, linetype = 2, alpha = 0.2) +
    
    geom_line() +
    scale_x_continuous(breaks = seq(0,6, by = 1)) +
    scale_y_continuous(breaks = seq(0,2, by = 0.5)) +
    theme_classic() +
    facet_wrap(~k_f, ncol = 1, scales = "free", labeller = label_parsed) +
    theme(aspect.ratio = 1,
          legend.position = "bottom",
          strip.background = element_blank(),
          strip.placement = "outside",
          legend.text = element_text(hjust = 0),
          panel.spacing = unit(1.15, "lines")
    ) +
    scale_fill_brewer(breaks = c("Expected", "Observed"),
                      labels = c("Expected", expression(paste("Observed (", italic(k), " = 0.57)"))),
                      type = 'qual',
                      palette = 3,
                      direction = -1) +
    scale_colour_brewer(breaks = c("Expected", "Observed"),
                        labels = c("Expected", expression(paste("Observed (", italic(k), " = 0.57)"))),
                        type = 'qual',
                        palette = 3,
                        direction = -1) +
    labs(fill = "",
         colour = "",
         x = expression(paste("Individual reproductive number (", italic(v), ")")),
         y = "Gamma probability distribution function") +
    guides(y = "prism_offset",
           x = "prism_offset") +
    coord_cartesian(expand = TRUE,
                    clip = "on",
                    xlim = c(0,6),
                    ylim = c(0,2))
    
  }) ()

### 3c
pC <- data_1c |> (function(x) {
  
  x |> 
    filter(k == 1 | k == 0.57) |>
    mutate(control = factor(control, levels = c("Controlled", "Uncontrolled"), labels = c("Observed", "Expected"))) |>
    ggplot(aes(x = proportion_of_infectious, y = exp_transmission, colour = factor(control))) +
    geom_abline(slope = 1, linetype = 2, alpha = 0.2) +
    geom_vline(xintercept = 0.2, linetype = 2, alpha = 0.2) +
    geom_line() +
    theme_classic() +
    facet_wrap(~k_f, ncol = 1, scales = "free", labeller = label_parsed) +
    theme(aspect.ratio = 1,
          legend.position = "bottom",
          strip.background = element_blank(),
          strip.placement = "outside",
          legend.text = element_text(hjust = 0),
          panel.spacing = unit(1.15, "lines")
    ) +
    scale_x_continuous(breaks = seq(0,1, by = 0.2)) +
    scale_y_continuous(breaks = seq(0,1, by = 0.2)) +
    scale_fill_brewer(breaks = c("Expected", "Observed"),
                      labels = c("Expected", expression(paste("Observed (", italic(k), " = 0.57)"))),
                      type = 'qual',
                      palette = 3,
                      direction = -1) +
    scale_colour_brewer(breaks = c("Expected", "Observed"),
                        labels = c("Expected", expression(paste("Observed (", italic(k), " = 0.57)"))),
                        type = 'qual',
                        palette = 3,
                        direction = -1) +
    labs(fill = "",
         colour = "",
         x = "Proportion of infectious cases (ranked)",
         y = "Proportion of transmission") +
    guides(y = "prism_offset",
           x = "prism_offset") +
    coord_cartesian(expand = TRUE, xlim = c(0,1), ylim = c(0,1)) +
    annotate(geom = "text",
             angle = 45,
             size = 2.5,
             x = 0.5,
             y = 0.5,
             vjust = 1.5,
             label = "Homogeneous population")
  
}) ()
  
### 3d
pD <- data_1d |> (function(x) { 
  x |> 
    filter(str_detect(k_h, "0.57")) |> 
    mutate(control = factor(control, levels = c("Controlled", "Uncontrolled"), labels = c("Observed", "Expected"))) |> 
    ggplot(aes(fill = control, x = generation)) +
    geom_line(aes(y = q, colour = control)) +
    geom_point(aes(y = gs, colour = control), alpha = 0.7, shape = 21) +
    theme_classic() +
    facet_wrap(~k_f, ncol = 1, scales = "free", labeller = label_parsed) +
    scale_x_continuous(breaks = seq(1,10, by = 1)) +
    scale_y_continuous(breaks = seq(0,0.8, by = 0.2)) +
    theme(aspect.ratio = 1,
          legend.position = "bottom",
          strip.background = element_blank(),
          strip.placement = "outside",
          legend.text = element_text(hjust = 0),
          panel.spacing = unit(1.15, "lines")
    ) +
    scale_fill_brewer(breaks = c("Expected", "Observed"),
                      labels = c("Expected", expression(paste("Observed (", italic(k), " = 0.57)"))),
                      type = 'qual', 
                      palette = 3,
                      direction = -1) +
    scale_colour_brewer(breaks = c("Expected", "Observed"),
                        labels = c("Expected", expression(paste("Observed (", italic(k), " = 0.57)"))),
                        type = 'qual', 
                        palette = 3,
                        direction = -1) +
    labs(fill = "",
         colour = "",
         x = expression(paste("Generation, ", 
                              italic(n))),
         y = expression(paste("Probability stochastic extinction, ", 
                              italic(q[n])))
    ) +
    guides(y = "prism_offset",
           x = "prism_offset") +
    coord_cartesian(expand = TRUE, 
                    clip = "on", 
                    xlim = c(1,10), 
                    ylim = c(0,0.8))
    }) ()
  

ABCD <- plot_grid(pA, pB, pC, pD, ncol = 4, labels = "auto", align = 'hv', axis = 'tblr')






### RT COME BACK TO THIS, FIND A WAY TO SMOOTH ESTIMATES IN A CLEAN WAY
## for now plot blanks



## READ data
rc_kc_inf <- read_rds(file = "code/simulations/data/out_rc_kc.rds")

nlim <- 150
rc_kc_est <- map(rc_kc_inf, function(x) {
  
  x |> 
    group_by(t_start, t_end, window, k_input) |> 
    mutate(n = n()) |> 
    filter(n >= nlim) |> 
    reframe(
      Rc =  mean(r, na.rm = TRUE),
      Rc_lower = quantile(r, 0.025, na.rm = TRUE),
      Rc_upper = quantile(r, 0.975, na.rm = TRUE),
      kc =  quantile(k, 0.5, na.rm = TRUE),
      kc_lower = quantile(k, 0.05, na.rm = TRUE),
      kc_upper = quantile(k, 0.95, na.rm = TRUE)
    ) |> 
  group_by(window, k_input) |>
  mutate(
    kc = rollapply(kc,
                   width = 3, by = 1, FUN = mean, align = "center", fill = "extend"),
    kc_upper = rollapply(kc_upper,
                         width = 3, by = 1, FUN = mean, align = "center", fill = "extend"),
    kc_lower = rollapply(kc_lower,
                         width = 3, by = 1, FUN = mean, align = "center", fill = "extend"),
    Rc = rollapply(Rc,
                   width = 3, by = 1, FUN = mean, align = "center", fill = "extend"),
    Rc_upper = rollapply(Rc_upper,
                         width = 3, by = 1, FUN = mean, align = "center", fill = "extend"),
    Rc_lower = rollapply(Rc_lower,
                         width = 3, by = 1, FUN = mean, align = "center", fill = "extend")
  ) |>
  ungroup() |>
  mutate(window_t = case_when(str_detect(window, "daily") ~ 1,
                              str_detect(window, "weekly") ~ 7,
                              TRUE ~ 14)) |> 
    mutate(t = case_when(window_t == 1 ~ t_start,
                         TRUE ~ t_start + (window_t/2)))
  
  
})

pF <- rbindlist(rc_kc_est) |> (function(x) {
  
  pal <- RColorBrewer::brewer.pal(n = 9, name = "Reds")
  
  R_input <- tibble(t = 0:120,
                    R = c(rep(2, 60), rep(0.5, 61)))
  
  x |> 
    mutate(k_input_f = factor(k_input, labels = k_labels)) |> 
    filter(k_input == 1) |> 
    ggplot() +
    geom_vline(xintercept = 60, linetype = 2, alpha = 0.2) +
    geom_line(data = R_input,
              mapping = aes(x = t, y = R)) +
    geom_ribbon(aes(x = t, ymax = Rc_upper, ymin = Rc_lower, fill = factor(window_t)), alpha = 0.1) +

    geom_line(aes(x = t, y = Rc, colour = factor(window_t))) + 
    facet_wrap(~k_input_f, ncol = 1, scales = "free", labeller = label_parsed) +
    scale_x_continuous(limits = c(0, 120), expand = c(0, 0), breaks = seq(0, 120, by = 10)) +
    # scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
    theme_classic() +
    scale_color_manual(values = c(pal[4], pal[6], pal[7])) +
    scale_fill_manual(values = c(pal[4], pal[6], pal[7])) +
    labs(fill = "Window (t)",
         colour = "Window (t)",
         x = expression(paste(italic(t))),
         y = expression(paste("Reproductive number, ", italic(R), " (", italic(hat(R[t])), ")"))) +
    theme(panel.spacing = unit(1.5, "lines"),
          legend.position = "bottom",
          strip.text = element_blank(),
          strip.background = element_blank(), 
          strip.placement = "outside") +
    coord_cartesian(expand = TRUE, xlim = c(0, 120), ylim = c(0, 5)) +
  # guides(y = "prism_offset") +
  facet_wrap(~window_t, ncol = 3, scales = "free")
  
})() ## Rt

pG <- rbindlist(rc_kc_est) |> (function(x) {
  
  x |> 
    mutate(k_input_f = factor(k_input, labels = k_labels)) |> 
    filter(k_input == 1) |> 
    ggplot() +
    geom_vline(xintercept = 60, linetype = 2, alpha = 0.2) +
    geom_hline(yintercept = 1) +
    geom_ribbon(aes(x = t, ymax = kc_upper, ymin = kc_lower, fill = factor(window_t)), alpha = 0.1) +
    # geom_smooth(aes(x = t, y = kc, colour = factor(window_t))) +
    geom_line(aes(x = t, y = kc, colour = factor(window_t))) + 
    facet_wrap(~k_input_f, ncol = 1, scales = "free", labeller = label_parsed) +
    scale_x_continuous(limits = c(0, 120), expand = c(0, 0), breaks = seq(0, 120, by = 10)) +
    # scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10() +
    theme_classic() +
    scale_color_manual(values = colorRampPalette(brewer.pal(9, "YlGnBu"))(12)[6:12]) +
    scale_fill_manual(values = colorRampPalette(brewer.pal(9, "YlGnBu"))(12)[6:12]) +
    labs(fill = "Window (t)",
         colour = "Window (t)",
         x = expression(paste(italic(t))),
         y = expression(paste("Estimated time-varying ", italic(k), " (", italic(hat(k[t])), ")"))) +
    theme(panel.spacing = unit(1.5, "lines"),
          legend.position = "bottom",
          strip.text = element_blank(),
          strip.background = element_blank(), 
          strip.placement = "outside") +
    coord_cartesian(expand = TRUE, xlim = c(0, 120), ylim = c(0.1, 100)) +
  # guides(y = "prism_offset") +
  facet_wrap(~window_t, ncol = 3, scales = "free")
  
})() ## Kt

### fit a smooth estimates





EFG <- plot_grid(pE, pF, pG, ncol = 1, labels = c("e", "f", "g"), align = 'hv', axis = 'tblr')

EFG

### Now to show, window period


# bind_rows(plot_data_sliding, plot_data_daily)
pal <- RColorBrewer::brewer.pal(12, name = "Paired")
pal_k <- pal[9:10]  
pal_R <- pal[5:6]  
pal_w <- pal[1:2]


pH <- plot_data_sliding |> (function(x) {
  
  x |> 
    filter(cases == "Controlled", x >= 39 & x <= 81) |>
    ggplot() +
    geom_vline(xintercept = 60, linetype = 2, alpha = 0.3) +
    geom_ribbon(aes(y = y, 
                    xmin = x_start,
                    xmax = x_end, 
                    fill = factor(window_t), 
                    colour = factor(window_t)),
                alpha = 0.1,
                colour = NA,
                # show.legend = FALSE
    ) +
    geom_line(aes(x = x, y = y, group =  factor(window_t), colour = factor(window_t)),
              size = 0.75,
    ) +
    # geom_line(aes(x = x_start, y = y, colour = factor(window)), linetype = 3, alpha = 0.7) +
    # geom_line(aes(x = x_end, y = y, colour =factor(window)), linetype = 3, alpha = 0.7) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      panel.spacing = unit(1.5, "lines"),
      aspect.ratio = 1,
      strip.placement = "outside"
    ) +
    scale_x_continuous(breaks = seq(39, 81, by = 7)) +
    coord_cartesian(xlim = c(40, 80)) +
    guides(
      y = "prism_offset"
      # x = "prism_offset"
    ) +
    scale_color_manual(values = pal_w) +
    scale_fill_manual(values = pal_w) +
    # scale_color_manual(values = colorRampPalette(brewer.pal(9, "YlGnBu"))(12)[5:12]) +
    labs(
      y = expression(paste(italic(P[CONTROL]))), 
      x = expression(paste("Time (", italic(t), ")")),
      colour = "Window length (days)",
      fill = "Window length (days)"
    )
}) () ## from  fig1.R files
 
pI <- bind_rows(cases_slinding14_m_prop_group, cases_slinding7_m_prop_group) |> (function(x) {
  
  x |>
    filter(!is.na(k_error)) |>
    ggplot(aes(
      y = r,
      x = factor(prop_cut, labels = c("0", "0 < 1", "1")),
      # colour = factor(window_t),
      fill = factor(window_t)
    )) +
    geom_hline(yintercept = c(2, 0.5), linetype = 2, alpha = 0.3) +
    # geom_boxplot(size = 0.5, alpha = 0.9, outlier.shape = NA) +
    geom_point(
      data = jitter_sample,
      aes(
        y = r,
        x = factor(prop_cut, labels = c("0", "0 < 1", "1")),
        group = factor(window_t),
        colour = factor(window_t),
        fill = factor(window_t)
      ),
      shape = 21,
      alpha = 0.5,
      position = position_jitterdodge(jitter.width = 0.15),
      show.legend = F
    ) +
    geom_boxplot(
      alpha = 0.9, 
      outlier.shape = NA
    ) +
    theme_classic() +
    scale_y_continuous(breaks = seq(0, 3.5, by = 0.5)) +
    theme(
      legend.position = "bottom",
      panel.spacing = unit(1.5, "lines"),
      aspect.ratio = 1.75,
      strip.placement = "outside"
    ) +
    guides(y = "prism_offset") +
    scale_colour_manual(values = pal_R) +
    scale_fill_manual(values = pal_R) +
    labs(
      y = expression(paste(italic(hat(R[t])))),
      x = expression(paste(italic(P[CONTROL]))),
      fill = "Window length (days)"
    ) +
    annotate("text", size = 3, 
             x = 3.8,
             y = 2.01,
             label = expression(paste(italic(t), " < 60"))) +
    annotate("text", size = 3, 
             x = 3.8,
             y = 0.51,
             label = expression(paste(italic(t), " > 60"))) +
    coord_cartesian(xlim = c(1, 3), ylim = c(0, 3.5), clip = "off")
  
}) ()
pJ <- bind_rows(cases_slinding14_m_prop_group, cases_slinding7_m_prop_group) |> (function(x) {
  
  x |> 
    filter(!is.na(k_error)) |>
    ggplot(aes(
      y = k,
      x = factor(prop_cut, labels = c("0", "0 < 1", "1")),
      # colour = factor(window_t),
      fill = factor(window_t)
    )) +
    geom_hline(yintercept = 1, linetype = 2, alpha = 0.3) +
    # geom_boxplot(size = 0.5, alpha = 0.9, outlier.shape = NA) +
    geom_point(
      data = jitter_sample,
      aes(
        y = k,
        x = factor(prop_cut, labels = c("0", "0 < 1", "1")),
        group = factor(window_t),
        colour = factor(window_t),
        fill = factor(window_t)
      ),
      shape = 21,
      alpha = 0.5,
      position = position_jitterdodge(jitter.width = 0.15),
      show.legend = F
    ) +
    geom_boxplot(
      alpha = 0.9, 
      outlier.shape = NA
    ) +
    scale_y_log10(breaks = c(0.1, 0.3, 1, 3, 10)) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      panel.spacing = unit(1.5, "lines"),
      aspect.ratio = 1.75,
      strip.placement = "outside"
    ) +
    guides(y = "prism_offset") +
    scale_colour_manual(values = pal_k) +
    scale_fill_manual(values = pal_k) +
    # scale_colour_brewer(palette = "Paired") +
    # scale_fill_brewer(palette = "Paired") +
    # scale_color_manual(values = colorRampPalette(brewer.pal(9, "YlGnBu"))(12)[5:12]) +
    # scale_fill_manual(values = colorRampPalette(brewer.pal(9, "YlGnBu"))(12)[5:12]) +
    labs(
      y = expression(paste(italic(hat(k[t])))),
      x = expression(paste(italic(P[CONTROL]))),
      fill = "Window length (days)"
    ) +
    annotate("text", size = 3, 
             x = 3.8,
             y = 1.01,
             label = expression(paste(italic(t), " > 0"))) +
    coord_cartesian(xlim = c(1, 3), ylim = c(0.1, 10), clip = "off")
  
}) ()
 

HIJ <- plot_grid(pH, pI, pJ, ncol = 1, labels = c("h", "i", "j"), align = 'hv', axis = 'tblr')


EFG_HIJ <- plot_grid(EFG, HIJ, ncol = 2, rel_widths = c(3,1))


p3 <- plot_grid(ABCD, EFG_HIJ, ncol = 1, rel_heights = c(1,3))

save_plot(plot = p3, filename = "plots/Fig.3.pdf", base_height = 13, base_width = 13)



