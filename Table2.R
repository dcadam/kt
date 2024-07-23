###### MAIN TABLE 2

library(fitdistrplus)
library(tidyverse)
library(data.table)
library(epicontacts)
library(VGAM)
library(insight)
library(meta)

set.seed(123)

source("code/simulations/util.R")


COVID_dated <- read_rds(file = "code/empirical/finished scripts/data/covid_dated_offspring.rds")
SARS_dated <- read_rds(file = "code/empirical/finished scripts/data/sars_dated_offspring.rds")


c_names <- names(COVID_dated)
s_names <- names(SARS_dated)

for (i in 1:2) {
  
  COVID_dated[[i]]$dataset <- c_names[[i]]
  SARS_dated[[i]]$dataset <- s_names[[i]]

}

covid_offspring <- rbindlist(COVID_dated) |> 
  mutate(period = case_when(
    t_inf >= "2020-01-20" & t_inf <= "2020-05-01" ~ "COVID Wave 1",
    t_inf >= "2020-06-20" & t_inf <= "2020-10-24" ~ "COVID Wave 2",
    t_inf >= "2020-10-25" ~ "COVID Wave 3",
    t_inf < "2020-01-01" ~ "SARS Wave",
    TRUE ~ NA_character_)) |> 
  filter(!is.na(t_inf)) |> 
  group_split(dataset, period)

sars_offspring <- rbindlist(SARS_dated) |> 
  mutate(period = case_when(
    t_inf >= "2020-01-20" & t_inf <= "2020-05-01" ~ "COVID Wave 1",
    t_inf >= "2020-06-20" & t_inf <= "2020-10-24" ~ "COVID Wave 2",
    t_inf >= "2020-10-25" ~ "COVID Wave 3",
    t_inf < "2020-01-01" ~ "SARS Wave",
    TRUE ~ NA_character_)) |> 
  filter(!is.na(t_inf)) |> 
  group_split(dataset, period)


#### ZINB by wave
pstr0 <- c(0, 0.1, 0.2, 0.3)

set.seed(123)
covid_tbl <- map(covid_offspring, function(x) {
  
    fit_list <- vector(mode = "list", length = length(pstr0))
    for (i in 1:length(pstr0)) {
      

        fit <- fitdist(x$offspring_count, 
                       distr = "zinegbin", 
                       start = list(munb = 1, size = 1), 
                       fix.arg = list(pstr0 = pstr0[[i]]))
        
        
        fit_list[[i]] <- rbind(
          data.table(
            pstr0 = pstr0[[i]], 
            var = "R",
            mean = unname(fit$estimate[1]),
            se = unname(fit$sd[1]),
            n = length(x$offspring_count),
            period = unique(x$period),
            dataset = unique(x$dataset)
          ),
          data.table(
            pstr0 = pstr0[[i]], 
            var = "k",
            mean = unname(fit$estimate[2]),
            se = unname(fit$sd[2]),
            n = length(x$offspring_count),
            period = unique(x$period),
            dataset = unique(x$dataset)
          )
        ) |> 
          mutate(sd = se * sqrt(n),
                 n = as.double(n))
        
    }
    
    rbindlist(fit_list)
    
}) |> 
  rbindlist() |> 
  group_split(var, dataset, pstr0) |> 
  (function(x) map(x, function(y) {

      mean.m <- metamean(data = y,
                         n = n,
                         mean = mean,
                         sd = sd,
                         sm = "MLN", 
                         common = FALSE)
      
      data.table(
        pstr0 = unique(y$pstr0),
        var = unique(y$var),
        value = exp(mean.m$TE.random),
        CI_low = exp(mean.m$lower.random),
        CI_high = exp(mean.m$upper.random),
        dataset = unique(y$dataset)
        )
      
  })
  ) () |> 
  rbindlist()


sars_tbl <- map(sars_offspring, function(x) {

    fit_list <- vector(mode = "list", length = length(pstr0))
  for (i in 1:length(pstr0)) {
    
    fit <- fitdist(x$offspring_count, 
                   distr = "zinegbin", 
                   start = list(munb = 1, size = 1), 
                   fix.arg = list(pstr0 = pstr0[[i]]))
    
    boot <- bootdist(fit, niter = 1000, bootmethod = 'nonparam', parallel = 'multicore', ncpus = 4)
    
    size <- quantile(boot$estim$size, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
    mu <- quantile(boot$estim$munb, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
    
    fit_list[[i]] <- rbind(
      data.table(
        pstr0 = pstr0[[i]],
        var = "R",
        value = mu[2],
        CI_low = mu[1],
        CI_high = mu[3],
        dataset = unique(x$dataset)
        ),
      data.table(
        pstr0 = pstr0[[i]],
        var = "k",
        value = size[2],
        CI_low = size[1],
        CI_high = size[3],
        dataset = unique(x$dataset)
      )
    )
  
  }
  
  rbindlist(fit_list)
  
}) |>  
  rbindlist()



out <- bind_rows(covid_tbl, sars_tbl) |> 
insight::format_table(ci_brackets = c("(",")")) |> 
  unite("res", value:CI, sep = " ") |> 
  pivot_wider(names_from = c("var"), values_from = res) |> 
  arrange(dataset, pstr0)


write_csv(out, file = "code/empirical/finished scripts/data/zinb_fits.csv")


knitr::kable(
  out,
  booktabs = TRUE)


write_csv(out_p80, file = "code/empirical/finished scripts/data/zinb_p80.csv")

bind_rows(
  sars_tbl |> 
  pivot_wider(names_from = var, values_from = c(value, CI_low, CI_high)) |> 
  mutate(joint_median = propresponsible(value_R, value_k, 0.8),
         joint_lower = propresponsible(CI_low_R, CI_low_k, 0.8),
         joint_upper = propresponsible(CI_high_R, CI_high_k, 0.8),
         fixed_median = propresponsible(3, value_k, 0.8),
         fixed_lower = propresponsible(3, CI_low_k, 0.8),
         fixed_upper = propresponsible(3, CI_high_k, 0.8)) |> 
  dplyr::select(pstr0, dataset, joint_median:fixed_upper) |> 
  mutate(across(joint_median:fixed_upper, ~ round(x = .x*100, digits = 0))) |> 
  pivot_longer(cols = joint_median:fixed_upper, names_to = "fit_method") |> 
  separate(fit_method, into = c("fit_method", "estimate"), sep = "_") |> 
  pivot_wider(names_from = estimate, values_from = value) |> 
  rename(value = median, 
         CI_low = lower,
         CI_high = upper) |> 
  mutate(value = paste0(value, "%")) |> 
  format_table(ci_brackets = c("(",")"), digits = 0, ci_digits = 0) |> 
  unite("res", value:CI, sep = " ") |> 
  pivot_wider(names_from = fit_method, values_from = res),
  covid_tbl |> 
  pivot_wider(names_from = var, values_from = c(value, CI_low, CI_high)) |> 
  mutate(joint_median = propresponsible(value_R, value_k, 0.8),
       joint_lower = propresponsible(CI_low_R, CI_low_k, 0.8),
       joint_upper = propresponsible(CI_high_R, CI_high_k, 0.8),
       fixed_median = propresponsible(2, value_k, 0.8),
       fixed_lower = propresponsible(2, CI_low_k, 0.8),
       fixed_upper = propresponsible(2, CI_high_k, 0.8)) |> 
  dplyr::select(pstr0, dataset, joint_median:fixed_upper) |> 
  mutate(across(joint_median:fixed_upper, ~ round(x = .x*100, digits = 0))) |> 
  pivot_longer(cols = joint_median:fixed_upper, names_to = "fit_method") |> 
  separate(fit_method, into = c("fit_method", "estimate"), sep = "_") |> 
  pivot_wider(names_from = estimate, values_from = value) |> 
  rename(value = median, 
         CI_low = lower,
         CI_high = upper) |> 
  mutate(value = paste0(value, "%")) |> 
  format_table(ci_brackets = c("(",")"), digits = 0, ci_digits = 0) |> 
  unite("res", value:CI, sep = " ") |> 
  pivot_wider(names_from = fit_method, values_from = res)
) |> knitr::kable(booktabs = TRUE)

