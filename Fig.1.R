###### MAIN FIGURE 1

library(tidyverse)
library(scales)
library(cowplot)
library(data.table)

##read data 
epicurve <- read_csv(file = "data/hk/hk-epicurves.csv")
rtk <- read_csv(file = "data/hk/hk-Rt-kt.csv")
epinow <- read_csv(file = "data/hk/hk-epinow2.csv")

COVID <- read_rds(file = "data/hk/rds/covid_dated_offspring.rds")
SARS <- read_rds(file = "data/hk/rds/sars_dated_offspring.rds")

offspring <- list(COVID[[1]], COVID[[2]], SARS[[1]], SARS[[2]])
offspring_names <- c(names(COVID), names(SARS))
names(offspring) <- offspring_names


pal <- RColorBrewer::brewer.pal(12, "Paired") # palette for plots
pal_epi <- RColorBrewer::brewer.pal(9, "BuGn") # palette for plots

figure_1 <- list() #empty list for grip

## Fig.1a
pA <- epicurve |> (function(x) {
  ggplot(x) +
    geom_histogram(aes(x = t, fill = sse_case_classification),
                   colour = "white",
                   linewidth = 0.3,
                   binwidth = 3
    ) +
    geom_hline(yintercept = 1, linetype = 2, alpha = 0) +
    scale_x_date("Symptom onset date",
                 date_breaks = "1 month",
                 date_labels = "%b %y"
    ) +
    scale_y_continuous(
      expand = c(0, 0)
    ) +
    scale_fill_manual(
  
      values = c(pal_epi[[1]], pal_epi[[3]], pal_epi[[5]], pal_epi[[7]])) +
    theme_classic() +
    theme(aspect.ratio = 0.75,
          legend.position = "bottom",
          strip.background = element_blank(),
          strip.placement = "outside",
          legend.text = element_text(hjust = 0),
          panel.spacing = unit(1.15, "lines")
    ) +
    facet_wrap(~period, ncol = 4, scales = "free") +
    labs(
      y = expression(paste("Confirmed cases, ", italic(n))),
      fill = "Case classification")
  
})()

## Fig.1b
pB <- map2(offspring, offspring_names, function(x, y) {
  
  x |> 
    mutate(period = case_when(
      t_inf >= "2020-01-20" & t_inf <= "2020-05-01" ~ "COVID Wave 1",
      t_inf >= "2020-06-20" & t_inf <= "2020-10-24" ~ "COVID Wave 2",
      t_inf >= "2020-10-25" ~ "COVID Wave 3",
      t_inf < "2020-01-01" ~ "SARS Wave",
      TRUE ~ NA_character_
    )) |> 
    filter(!is.na(period)) |> 
    mutate(sensitivity = y)}
  ) |> rbindlist() |>
  (function(x) {
      
      x |> 
      mutate(sensitivity = case_when(str_detect(sensitivity, pattern = "primary") ~ "Sensitivity 1",
                                     TRUE ~ "Sensitivity 2")) |> 
        mutate(dist_change = case_when(offspring_count >= 10 ~ 10,
                                       TRUE ~ offspring_count)) |> 
        ggplot(aes(x = dist_change, colour  = sensitivity, fill = sensitivity, y = after_stat(prop))) +
        geom_bar(position = position_dodge2(padding = 0.2), width = 0.9, alpha = 0.75) +
        facet_wrap(~k_f, ncol = 1, scales = "free", labeller = label_parsed) +
        scale_x_continuous(limits = c(-0.5,10.5), breaks = c(0:10), labels = c(as.character(0:9), "10+")) +
        scale_y_continuous(limits = c(0,1)) +
        theme_classic() +
        theme(aspect.ratio = 0.75,
              legend.position = "bottom",
              strip.background = element_blank(),
              strip.placement = "outside",
              legend.text = element_text(hjust = 0),
              panel.spacing = unit(1.15, "lines")
        ) +
        scale_fill_manual(labels = c("Primary", "Sensitivity"), values = c(pal[10], pal[9])) +
        scale_color_manual(labels = c("Primary", "Sensitivity"), values = c(pal[10], pal[9])) +
        labs(fill = "Dataset",
             colour = "Dataset",
             x = expression(paste("Number of secondary cases, ", italic(Z))), 
             y = "Proportion of secondary cases") +
      coord_cartesian(expand = TRUE) +
        facet_wrap(~period, ncol = 4, scales = "free")
      
    }) ()

## Fig.1c
pC <- rtk |> (function(x) {
  align_df <- bind_rows(
    epicurve |>
      group_by(period) |>
      reframe(t = max(t)) |>
      mutate(
        k = c(100, 100, 100, 100),
        r = c(4, 4, 4, 35)
      ),
    epicurve |>
      group_by(period) |>
      reframe(t = min(t)) |>
      mutate(
        k = c(0.005, 0.005, 0.005, 0.005),
        r = c(0, 0, 0, 0)
      )
  )
  
  x |>
    filter(sensitivity == "Sensitivity 1") |> 
    mutate(method = "Z") |>
    bind_rows(epinow) |>
    ggplot() +
    geom_point(data = align_df, aes(x = t, y = r), colour = "white") +
    geom_hline(yintercept = 1, alpha = 0.5, linetype = 2) +
    geom_ribbon(aes(x = t, ymax = rt_upper, ymin = rt_lower, fill = method), alpha = 0.1) +
    geom_line(aes(x = t, y = rt, colour = method, linetype = method)) +
    scale_x_date(
      name = expression(paste(italic(t))),
      date_breaks = "1 month",
      date_labels = "%b %y",
      minor_breaks = NULL
    ) +
    theme_classic() +
    theme(aspect.ratio = 0.75,
          legend.position = "bottom",
          strip.background = element_blank(),
          strip.placement = "outside",
          legend.text = element_text(hjust = 0),
          panel.spacing = unit(1.15, "lines")
    ) +
    facet_wrap(~period, ncol = 4, scales = "free") +
    labs(
      colour = "Source",
      fill = "Source",
      linetype = "Source",
      y = expression(paste("Reproductive number, ", italic(R[t])))
    ) +
    scale_color_manual(
      labels = c("EpiNow2", "Offspring distribution (Primary)"),
      values = c(pal[5], pal[6])
    ) +
    scale_fill_manual(
      labels = c("EpiNow2", "Offspring distribution (Primary)"),
      values = c(pal[5], pal[6])
    ) +
    scale_linetype_manual(
      labels = c("EpiNow2", "Offspring distribution (Primary)"),
      values = c("dashed", "solid")
    )
})()

## Fig.1d
pD <- rtk |> (function(x) {
  align_df <- bind_rows(
    epicurve |>
      group_by(period) |>
      reframe(t = max(t)) |>
      mutate(
        k = c(100, 100, 100, 100),
        r = c(3, 3, 3, 35)
      ),
    epicurve |>
      group_by(period) |>
      reframe(t = min(t)) |>
      mutate(
        k = c(0.005, 0.005, 0.005, 0.005),
        r = c(0, 0, 0, 0)
      )
  )
  
  
  x |>
    ggplot() +
    geom_point(data = align_df, aes(x = t, y = k), colour = "white") +
    geom_hline(yintercept = 1, alpha = 0.5, linetype = 2) +
    geom_ribbon(aes(x = t, ymax = kt_upper, ymin = kt_lower, fill = sensitivity), alpha = 0.1) +
    geom_line(aes(x = t, y = kt, colour = sensitivity)) +
    scale_x_date(
      name = expression(paste(italic(t))),
      date_breaks = "1 month",
      date_labels = "%b %y",
      minor_breaks = NULL
    ) +
    scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
    theme_classic() +
    theme(aspect.ratio = 0.75,
          legend.position = "bottom",
          strip.background = element_blank(),
          strip.placement = "outside",
          legend.text = element_text(hjust = 0),
          panel.spacing = unit(1.15, "lines")
    ) +
    facet_wrap(~period, ncol = 4, scales = "free") +
    labs(
      colour = "Dataset",
      fill = "Dataset",
      shape = "",
      linetype = "",
      y = expression(paste("Dispersion, ", italic(k[t])))
    ) +
    scale_color_manual(values = c(pal[2], pal[1]), labels = c("Primary", "Sensitivity")) +
    scale_fill_manual(values = c(pal[2], pal[1]), labels = c("Primary", "Sensitivity"))
})()

# save Fig.1
p1 <- plot_grid(pA, pB, pC, pD, labels = "auto", ncol = 1, align = "hv", axis = "b")

save_plot(plot = p1, filename = "plots/Fig.1.png", base_height = 13, base_width = 13, dpi = 600)
save_plot(plot = p1, filename = "plots/Fig.1.pdf", base_height = 13, base_width = 13)

