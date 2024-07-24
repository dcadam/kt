###### MAIN TABLE 1

library(tidyverse)
library(zoo)
library(data.table)
library(meta)
library(metasens)
library(fitdistrplus)
library(insight)

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


format_table(ci_brackets = c("(",")")) |> 
  unite("res", value:CI, sep = " ") |> 
  pivot_wider(names_from = c("var"), values_from = res) |> 
  arrange(dataset, control)


format_table(ci_brackets = c("(",")")) |> 
  unite("res", value:CI, sep = " ") |> 
  pivot_wider(names_from = c("var"), values_from = res) |> 
  arrange(dataset, control)

write_csv(
  x = bind_rows(covid_tbl, sars_tbl),
  file = "output/table1.csv")

