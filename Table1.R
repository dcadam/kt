###### MAIN TABLE 1

library(tidyverse)
library(zoo)
library(data.table)
library(meta)
library(metasens)
library(fitdistrplus)
library(insight)

source('src/p80.R')
source('src/t20.R')

covid_sars_bc_ac <- read_rds(file = "data/hk/rds/hk-offspring-control.rds")

set.seed(123)
covid_tbl <- covid_sars_bc_ac$covid |> (function(x) {
  
  xl <- x |> 
    group_split(dataset, control, period)
  
  xl_fit <- map(xl, function(y) {
    
    fit <- fitdist(data = y$offspring_count, distr = 'nbinom')
    
    rbind(
      data.table(
        var = "R",
        mean = unname(fit$estimate[2]),
        se = unname(fit$sd[2]),
        n = length(y$offspring_count),
        period = unique(y$period),
        control = unique(y$control),
        dataset = unique(y$dataset)
      ),
      data.table(
        var = "k",
        mean = unname(fit$estimate[1]),
        se = unname(fit$sd[1]),
        n = length(y$offspring_count),
        period = unique(y$period),
        control = unique(y$control),
        dataset = unique(y$dataset)
      )
    ) 
  }) |> 
    rbindlist() |> 
    mutate(sd = se * sqrt(n),
           n = as.double(n)) |> 
    group_split(var, control, dataset)
  
  ### inverse variance analysis
  
  out <- map(xl_fit, function(y) {
    
    mean.m <- metamean(data = y,
                       n = n,
                       mean = mean,
                       sd = sd,
                       sm = "MLN", 
                       method.random.ci = "HK", 
                       common = FALSE)
    
    data.table(
      var = unique(y$var),
      value = exp(mean.m$TE.random),
      CI_low = exp(mean.m$lower.random),
      CI_high = exp(mean.m$upper.random),
      dataset = unique(y$dataset),
      control = unique(y$control)
    )
    
  }) |>
    rbindlist() |> 
    mutate(control = case_when(control == "bc" ~ "less",
                               control == "ac" ~ "more"))
 
  
  return(out)
  
  
}) ()


sars_tbl <- covid_sars_bc_ac$sars |> (function(x) {
  
  xl <- x |> 
    group_split(dataset, control, period)
  
  out <- map(xl, function(y) {
    
    fit <- fitdist(data = y$offspring_count, distr = 'nbinom')
    boot <- bootdist(fit, niter = 1000, bootmethod = 'nonparam', parallel = 'multicore', ncpus = 4)
    
    size <- quantile(boot$estim$size, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
    mu <- quantile(boot$estim$mu, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
    
    rbind(
      data.table(
        var = "R",
        CI_low = mu[1],
        value = mu[2],
        CI_high = mu[3],
        dataset = unique(y$dataset),
        control = unique(y$control)
      ),
      data.table(
        var = "k",
        CI_low = size[1],
        value = size[2],
        CI_high = size[3],
        dataset = unique(y$dataset),
        control = unique(y$control)
      )
    ) 
  }, .progress = TRUE) |> 
    rbindlist() |> 
    mutate(control = case_when(control == "bc" ~ "less",
                               control == "ac" ~ "more"))
   
  
  return(out)
  
  
}) ()

## write table
out <- bind_rows(covid_tbl, sars_tbl) |> 
format_table(ci_brackets = c("(",")")) |> 
  unite("res", value:CI, sep = " ") |> 
  pivot_wider(names_from = c("var"), values_from = res) |> 
  arrange(dataset, control)


write_csv(out,
          file = "output/table1.csv")



#### calculate p80


p80 <- bind_rows(covid_tbl, sars_tbl) |> 
  pivot_wider(names_from = var, values_from = c(value, CI_low, CI_high)) |> 
  mutate(value = propresponsible(value_R, value_k, 0.8),
         CI_low = propresponsible(CI_low_R, CI_low_k, 0.8),
         CI_high = propresponsible(CI_high_R, CI_high_k, 0.8)) |> 
  dplyr::select(dataset, control, value:CI_high) |> 
  mutate(across(value:CI_high, ~ round(x = .x*100, digits = 1))) |> 
  mutate(value = paste0(value, "%")) |> 
  format_table(ci_brackets = c("(",")"), digits = 1, ci_digits = 1) |> 
  unite("p80", value:CI, sep = " ") |> 
  arrange(dataset, control)




#### t20


t20 <- bind_rows(covid_tbl, sars_tbl) |> 
  pivot_wider(names_from = "var", values_from = value:CI_high) |> 
  arrange(dataset, control)


  data.table(
    dataset = t20$dataset, 
    control = t20$control,
    t20 = map2_dbl(.x = t20$value_R, .y = t20$value_k, function(x, y) {
      proportion_transmission2(R = x, k = y, prop = 0.2)
    }),
    CI_high = map2_dbl(.x = t20$CI_low_R, .y = t20$CI_low_k, function(x, y) {
      proportion_transmission2(R = x, k = y, prop = 0.2)
    }),
    CI_low = map2_dbl(.x = t20$CI_high_R, .y = t20$CI_high_k, function(x, y) {
      proportion_transmission2(R = x, k = y, prop = 0.2)
    })
  ) |> 
      left_join(p80, by = c("dataset", "control")) |> 
    mutate(across(t20:CI_low, ~ round(x = .x*100, digits = 1))) |> 
    mutate(t20 = paste0(t20, "%")) |> 
    format_table(ci_brackets = c("(",")"), digits = 1, ci_digits = 1) |> 
    unite("t20", t20:CI, sep = " ") |> 
    as_tibble() |> 
    knitr::kable()
  