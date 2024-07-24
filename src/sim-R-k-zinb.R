

library(tidyverse)
library(VGAM)
library(fitdistrplus)
library(data.table)
library(ggridges)
library(scales)
library(insight)

k_steps <- c(0.1, 0.5, 1, 10, 100)

### 300 epidemics, sample 500 observations. 
t_data <- read_rds(file = "data/simulations/rds/sim-output-raw.rds")

t_data_l <- map(t_data, function(x) {
  
  x |> 
    group_split(epidemic_id)
  
}) 

k_steps <- c(0.1, 0.5, 1, 10, 100)
out_nb <- vector(mode = "list", length = 5)
for (i in 1:5) {
  
  print(paste("ANALYZING DATASET", i, "..."))
  
  out_nb[[i]] <- map(t_data_l[[i]], .progress = TRUE, function(x) {
    
    obs <- x$obs_n
    
    obs <- sample(obs, size = 500)
    
    epi_id <- unique(x$epidemic_id)
    
    fit <- possfitdist(obs, distr = "nbinom", start = list(mu = 1, size = 1))
    
    est_fit <- fit$estimate
    
    est_names <- names(est_fit)
    
    fit_list <- tibble(est_names, est_fit) |>
      pivot_wider(names_from = "est_names",
                  values_from = "est_fit") |> 
      mutate(pstr0 = "nb", k = k_steps[[i]], R = 2, epidemic_id = epi_id, AIC = fit$aic) |> 
      as.data.table()
    
    return(fit_list) 
    
    
  }) |> 
    keep(function(x) length(x) > 1) |> 
    list_rbind() ## end of map
  
} 

set.seed(123)
out_nb_fixed <- vector(mode = "list", length = 5)
for (i in 1:5) {
  
  print(paste("ANALYZING DATASET", i, "..."))
  
  out_nb_fixed[[i]] <- map(t_data_l[[i]], .progress = TRUE, function(x) {
    
    obs <- x$obs_n
    
    obs <- sample(obs, size = 500)
    
    epi_id <- unique(x$epidemic_id)
    
    fit <- possfitdist(obs, distr = "nbinom", start = list(size = 1), fix.arg = list(mu = 2))
    
    est_fit <- c(mu = fit$fix.arg$mu, fit$estimate)
    
    est_names <- names(est_fit)
    
    fit_list <- tibble(est_names, est_fit) |>
      pivot_wider(names_from = "est_names",
                  values_from = "est_fit") |> 
      mutate(pstr0 = "nb", k = k_steps[[i]], R = 2, epidemic_id = epi_id, AIC = fit$aic) |> 
      as.data.table()
    
    return(fit_list) 
    
    
  }) |> 
    keep(function(x) length(x) > 1) |> 
    list_rbind() ## end of map
  
} ### DONE

set.seed(123)
out_zinb <- vector(mode = 'list', length = length(t_data))
for (i in 1:length(t_data)) {
  
  print(paste("ANALYZING DATASET", i, "..."))
  
  out_zinb[[i]] <- map(t_data_l[[i]], function(x) {
    
    obs <- x$obs_n
    
    obs <- sample(obs, size = 500, replace = FALSE) ## a sample of 500 cases what is the k distribution
    
    epi_id <- unique(x$epidemic_id)
    
    pstr0 <- seq(0, 0.5, by = 0.05) ##pstr0 is equivalent to nbinom
    
    #### fixed range of PSTR0 up to max 0.5 when R = 2
    fit_list <- vector(mode = "list", length = length(pstr0))
    for (k in 1:length(pstr0)) {
      
      fit <- possfitdist(obs, distr = "zinegbin", start = list(munb = 1, size = 1), fix.arg = list(pstr0 = pstr0[[k]]))
      
      if (fit[[1]][[1]] == "Error") {
        
        fit_list[[k]] <- fit
        
      } else {
        
        est_fit <- fit$estimate
        
        est_names <- names(est_fit)
        
        fit_list[[k]] <- tibble(est_names, est_fit) |>
          pivot_wider(names_from = "est_names",
                      values_from = "est_fit") |> 
          mutate(pstr0 = pstr0[[k]], k = k_steps[[i]], R = 2, epidemic_id = epi_id, AIC = fit$aic) |> 
          as.data.table()
        
        
      }
      
    }
    
    fit_list |> 
      keep(function(x) length(x) > 1) |> 
      list_rbind()
    
    
    
    
  }, .progress = TRUE) |> 
    keep(function(x) length(x) > 1) |> 
    list_rbind()
  
  
}

set.seed(123)
out_zinb_fixed <- vector(mode = 'list', length = length(t_data))
for (i in 1:length(t_data)) {
  
  print(paste("ANALYZING DATASET", i, "..."))
  
  out_zinb_fixed[[i]] <- map(t_data_l[[i]], function(x) {
    
    obs <- x$obs_n
    
    obs <- sample(obs, size = 500, replace = FALSE) ## a sample of 500 cases what is the k distribution
    
    epi_id <- unique(x$epidemic_id)
    
    pstr0 <- seq(0, 0.5, by = 0.05) ##pstr0 is equivalent to nbinom
    
    #### fixed range of PSTR0 up to max 0.5 when R = 2
    fit_list <- vector(mode = "list", length = length(pstr0))
    for (k in 1:length(pstr0)) {
      
      fit <- possfitdist(obs, distr = "zinegbin", start = list(size = 1), fix.arg = list(munb = 2, pstr0 = pstr0[[k]]))
      
      if (fit[[1]][[1]] == "Error") {
        
        fit_list[[k]] <- fit
        
      } else {
        
        est_fit <-  c(munb = fit$fix.arg$munb, fit$estimate)
        
        est_names <- names(est_fit)
        
        fit_list[[k]] <- tibble(est_names, est_fit) |>
          pivot_wider(names_from = "est_names",
                      values_from = "est_fit") |> 
          mutate(pstr0 = pstr0[[k]], k = k_steps[[i]], R = 2, epidemic_id = epi_id, AIC = fit$aic) |> 
          as.data.table()
        
        
      }
      
    }
    
    fit_list |> 
      keep(function(x) length(x) > 1) |> 
      list_rbind()
    
    
    
    
  }, .progress = TRUE) |> 
    keep(function(x) length(x) > 1) |> 
    list_rbind()
  
  
}

res <- list(nb = rbindlist(out_nb),
            nb_fixed = rbindlist(out_nb_fixed),
            zinb = rbindlist(out_zinb),
            zinb_fixed = rbindlist(out_zinb_fixed))

### calculate errors
res_error <- map(res, function(x) {
  
  x |> 
    mutate(k_input_f = factor(k, labels = c("k = 0.1", "k = 0.5", "k = 1", "k = 10", "k = 100")),
           pstr0_f = factor(pstr0),
           error_size = (log10(k) - log10(size))^2) |> 
    group_by(pstr0, k) |> 
    mutate(mse_size = mean(error_size)) |> 
    group_by(k) |> 
    mutate(mse_size_normal = normalize(mse_size)) |> 
    as.data.table()
  
  
})

### TABLE S2 (Negative Binomial)

res_error$zinb |> 
  filter(pstr0 == 0) |> 
  group_by(k) |> 
  reframe(mu_median = quantile(munb, 0.5),
          mu_low = quantile(munb, 0.025),
          mu_high = quantile(munb, 0.975),
          size_median = quantile(size, 0.5),
          size_low = quantile(size, 0.05),
          size_high = quantile(size, 0.95)) |> 
  as.data.table() |> 
  melt(id.vars = "k", 
       measure.vars = patterns(value = "median" ,
                               CI_low = "low",
                               CI_high = "high"),
       variable.name = "var") |> 
  mutate(var = case_when(var == 1 ~ "R", 
                         var == 2 ~ "size")) |> 
  format_table(ci_brackets = c("(",")")) |> 
  unite("res", value:CI, sep = " ") |> 
  pivot_wider(names_from = var, values_from = res) |> 
  knitr::kable()



## TABLE S2 (Zero-Inflated Negative Binomial)

res_error$zinb |> 
  group_by(k, pstr0) |> 
  reframe(mu_median = quantile(munb, 0.5),
          mu_low = quantile(munb, 0.025),
          mu_high = quantile(munb, 0.975),
          size_median = quantile(size, 0.5),
          size_low = quantile(size, 0.05),
          size_high = quantile(size, 0.95),
          mse = mean(mse_size_normal)) |> 
  group_by(k) |> 
  arrange(mse) |> 
  slice(1) |> 
  as.data.table() |> 
  melt(id.vars = c("k", "pstr0"), 
       measure.vars = patterns(value = "median" ,
                               CI_low = "low",
                               CI_high = "high"),
       variable.name = "var") |> 
  mutate(var = case_when(var == 1 ~ "R", 
                         var == 2 ~ "size")) |> 
  format_table(ci_brackets = c("(",")")) |> 
unite("res", value:CI, sep = " ") |> 
  pivot_wider(names_from = var, values_from = res) |> 
  knitr::kable()





