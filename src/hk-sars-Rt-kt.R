library(epicontacts)
library(zoo)
library(tidyverse)
library(fitdistrplus)
library(lubridate)


options(show.error.messages = TRUE) ## suppress expected MLE failures from console

input <- read_rds(file = "data/hk/rds/sars-dated-offspring.rds")

### DAILY t = 1

output_daily <- map(input, function(x) {
  
  vec_t <- x |> 
    arrange(t_inf) |> 
    pull(t_inf) |> 
    unique()
  
  tmp_dt <- as.data.table(x)
  
  estimate <- tibble()
  for (t in 1:length(vec_t)) {
    
    tryCatch( #error catch
      {
        
        set.seed(as.numeric(vec_t[t]))
        
        tmp <- tmp_dt[t_inf == vec_t[t]]
        
        fit_joint <- fitdist(tmp$offspring_count, 'nbinom')
        
        fit_joint_boot <- fit_joint %>%
          bootdist(parallel = "multicore", ncpus = ncpu, bootmethod = "nonparam")
        
        est_joint <- tibble(method = "joint", 
                            window = "daily",
                            k = unname(fit_joint_boot$estim$size %>%
                                         quantile(0.50, na.rm = TRUE)),
                            k_lower = unname(fit_joint_boot$estim$size %>%
                                               quantile(0.05, na.rm = TRUE)),
                            k_upper = unname(fit_joint_boot$estim$size %>%
                                               quantile(0.95, na.rm = TRUE)),
                            k_sd = fit_joint$sd[1],
                            k_boot = fit_joint_boot$estim$size %>% list(),
                            start = vec_t[t],
                            end = vec_t[t],
                            rt = unname(fit_joint_boot$estim$mu %>%
                                          quantile(0.50, na.rm = TRUE)),
                            rt_lower = unname(fit_joint_boot$estim$mu %>%
                                                quantile(0.05, na.rm = TRUE)),
                            rt_upper = unname(fit_joint_boot$estim$mu %>%
                                                quantile(0.95, na.rm = TRUE)),
                            r_sd = fit_joint$sd[2],
                            rt_boot = fit_joint_boot$estim$mu %>% list(),
                            n = nrow(tmp))
        
        
        
        est <- est_joint
        
        estimate <- bind_rows(estimate, est)
        
        print(est)
        
      },
      error = function(e) {}
    )
    
  }
  
  return(estimate)
  
  
  
})

output_daily <- imap(output_daily, ~ .x |> 
                       mutate(sensitivity = .y)) |> 
  list_rbind()


#### SLIDING ESTIMATION t = 7

output_sliding7 <- map(input, function(x) {
  
  vec_t <- x |> 
    arrange(t_inf) |> 
    transmute(start = t_inf, end = t_inf + 7) |> 
    distinct()
  
  tmp_dt <- as.data.table(x)
  
  
  estimate <- tibble()
  for (t in 1:nrow(vec_t)) {
    
    tryCatch( #error catch
      {
        
        set.seed(as.numeric(vec_t$start[t]))
        
        tmp <- tmp_dt[t_inf >= vec_t$start[t] & t_inf <= vec_t$end[t]]
        
        fit_joint <- fitdist(tmp$offspring_count, 'nbinom')
        
        fit_joint_boot <- fit_joint %>%
          bootdist(parallel = "multicore", ncpus = ncpu, bootmethod = "nonparam")
        
        est_joint <- tibble(method = "joint", 
                            window = "sliding7",
                            k = unname(fit_joint_boot$estim$size %>%
                                         quantile(0.50, na.rm = TRUE)),
                            k_lower = unname(fit_joint_boot$estim$size %>%
                                               quantile(0.05, na.rm = TRUE)),
                            k_upper = unname(fit_joint_boot$estim$size %>%
                                               quantile(0.95, na.rm = TRUE)),
                            k_sd = fit_joint$sd[1],
                            k_boot = fit_joint_boot$estim$size %>% list(),
                            start = vec_t$start[t],
                            end = vec_t$end[t],
                            rt = unname(fit_joint_boot$estim$mu %>%
                                          quantile(0.50, na.rm = TRUE)),
                            rt_lower = unname(fit_joint_boot$estim$mu %>%
                                                quantile(0.05, na.rm = TRUE)),
                            rt_upper = unname(fit_joint_boot$estim$mu %>%
                                                quantile(0.95, na.rm = TRUE)),
                            r_sd = fit_joint$sd[2],
                            rt_boot = fit_joint_boot$estim$mu %>% list(),
                            n = nrow(tmp))
        
        
        est <- est_joint
        
        estimate <- bind_rows(estimate, est)
        
        print(est)
        
      },
      error = function(e) {}
    )
    
  }
  
  return(estimate)
  
  
  
})

output_sliding7 <- imap(output_sliding7, ~ .x |> 
                          mutate(sensitivity = .y)) |> 
  list_rbind()


#### SLIDING ESTIMATION t = 14

output_sliding14 <- map(input, function(x) {
  
  vec_t <- x |> 
    arrange(t_inf) |> 
    transmute(start = t_inf, end = t_inf + 14) |> 
    distinct()
  
  tmp_dt <- as.data.table(x)
  
  estimate <- tibble()
  for (t in 1:nrow(vec_t)) {
    
    tryCatch( #error catch
      {
        
        set.seed(as.numeric(vec_t$start[t]))
        
        tmp <- tmp_dt[t_inf >= vec_t$start[t] & t_inf <= vec_t$end[t]]
        
        fit_joint <- fitdist(tmp$offspring_count, 'nbinom')
        
        fit_joint_boot <- fit_joint %>%
          bootdist(parallel = "multicore", ncpus = ncpu, bootmethod = "nonparam")
        
        est_joint <- tibble(method = "joint", 
                            window = "sliding14",
                            k = unname(fit_joint_boot$estim$size %>%
                                         quantile(0.50, na.rm = TRUE)),
                            k_lower = unname(fit_joint_boot$estim$size %>%
                                               quantile(0.05, na.rm = TRUE)),
                            k_upper = unname(fit_joint_boot$estim$size %>%
                                               quantile(0.95, na.rm = TRUE)),
                            k_sd = fit_joint$sd[1],
                            k_boot = fit_joint_boot$estim$size %>% list(),
                            start = vec_t$start[t],
                            end = vec_t$end[t],
                            rt = unname(fit_joint_boot$estim$mu %>%
                                          quantile(0.50, na.rm = TRUE)),
                            rt_lower = unname(fit_joint_boot$estim$mu %>%
                                                quantile(0.05, na.rm = TRUE)),
                            rt_upper = unname(fit_joint_boot$estim$mu %>%
                                                quantile(0.95, na.rm = TRUE)),
                            r_sd = fit_joint$sd[2],
                            rt_boot = fit_joint_boot$estim$mu %>% list(),
                            n = nrow(tmp))
        
        
        est <- est_joint
        
        estimate <- bind_rows(estimate, est)
        
        print(est)
        
      },
      error = function(e) {}
    )
    
  }
  
  return(estimate)
  
  
  
})

output_sliding14 <- imap(output_sliding14, ~ .x |> 
                           mutate(sensitivity = .y)) |> 
  list_rbind()


### SAVE DATA

bind_rows(output_daily, output_sliding7, output_sliding14) |> 
  saveRDS(file = "data/hk/rds/hk-covid-Rt-kt-raw.rds")



######## MIXTURE ESIMTATION

### DAILY t = 1

####### FOREACH MIXTURE

library(foreach)
library(doParallel)


ncpus <- parallel::detectCores()

cluster <- makeCluster(ncpus) ## start parralel cluster within each epidemic j to control memory
registerDoParallel(cluster)

output_daily <- map(input, function(x) {
  
  vec_t <- x |> 
    arrange(t_inf) |> 
    pull(t_inf) |> 
    unique()
  
  x <- as.data.table(x)
  
  est <- foreach(t = 1:length(vec_t), .packages=c('tidyr', 'flexmix', 'countreg', 'purrr', 'readr', 'dplyr', 'fitdistrplus', 'data.table', 'VGAM')) %dopar% {
    
    tryCatch( #error catch
      {
        
        set.seed(as.numeric(vec_t[t]))
        
        tmp <- x[t_inf == vec_t[t]]
        
        
        if (sum(tmp$offspring_count) == 0) {
          
          next
          
        }
        
        ## NB fit with random cluster assignment cannot determine clustering. 
        flex_fit_t <- stepFlexmix(x ~ 1, 
                                  data = data.table(x = tmp$offspring_count),
                                  k = 1:2, 
                                  model = FLXMRnegbin(), nrep = 5, verbose = FALSE)
        
        if (class(flex_fit_t)[1] == "stepFlexmix") {
          
          flex_fit_t <- getModel(flex_fit_t)
          
        }
        
        params_t <- parameters(flex_fit_t)
        
        if (flex_fit_t@k0 == 1) {
          
          mixture <- tibble(method = "mixture", 
                            window = "daily",
                            r1 = exp(params_t[1,1]),
                            r2 = NA,
                            k1 = params_t[2,1],
                            k2 = NA,
                            n1 = sum(flex_fit_t@cluster),
                            n2 = 0,
                            start = vec_t[t],
                            end = vec_t[t])
          
          
        } else {
          
          
          mixture <- tibble(method = "mixture", 
                            window = "daily",
                            r1 = exp(params_t[1,1]),
                            r2 = exp(params_t[1,2]),
                            k1 = params_t[2,1],
                            k2 = params_t[2,2],
                            n1 = flex_fit_t@size[1],
                            n2 =flex_fit_t@size[2],
                            start = vec_t[t],
                            end = vec_t[t])
          
        }
        
        
        return(mixture)
        
        
      },
      error = function(e) {}
    )
    
    
  }
  
  
  estimate <- est |> 
    keep(function(x) length(x) > 1) |> 
    list_rbind() 
  
  
  return(estimate)
  
  
})


output_daily <- imap(output_daily, ~ .x |> 
                       mutate(sensitivity = .y)) |> 
  list_rbind()


### SLIDING t = 7

output_sliding7 <- map(input, function(x) {
  
  vec_t <- x |> 
    arrange(t_inf) |> 
    transmute(start = t_inf, end = t_inf + 7) |> 
    distinct()
  
  x <- as.data.table(x)
  
  est <- foreach(t = 1:nrow(vec_t), .packages=c('tidyr', 'flexmix', 'countreg', 'purrr', 'readr', 'dplyr', 'fitdistrplus', 'data.table', 'VGAM')) %dopar% {
    
    tryCatch( #error catch
      {
        
        set.seed(as.numeric(vec_t$start[t]))
        
        tmp <- x[t_inf >= vec_t$start[t] & t_inf <= vec_t$end[t]]
        
        
        if (sum(tmp$offspring_count) == 0) {
          
          next
          
        }
        
        ## NB fit with random cluster assignment cannot determine clustering. 
        flex_fit_t <- stepFlexmix(x ~ 1, 
                                  data = data.table(x = tmp$offspring_count),
                                  k = 1:2, 
                                  model = FLXMRnegbin(), nrep = 5, verbose = FALSE)
        
        if (class(flex_fit_t)[1] == "stepFlexmix") {
          
          flex_fit_t <- getModel(flex_fit_t)
          
        }
        
        
        params_t <- parameters(flex_fit_t)
        
        if (flex_fit_t@k0 == 1) {
          
          mixture <- tibble(method = "mixture", 
                            window = "sliding7",
                            r1 = exp(params_t[1,1]),
                            r2 = NA,
                            k1 = params_t[2,1],
                            k2 = NA,
                            n1 = sum(flex_fit_t@cluster),
                            n2 = 0,
                            start = vec_t$start[t],
                            end = vec_t$end[t])
          
          
        } else {
          
          
          mixture <- tibble(method = "mixture", 
                            window = "sliding7",
                            r1 = exp(params_t[1,1]),
                            r2 = exp(params_t[1,2]),
                            k1 = params_t[2,1],
                            k2 = params_t[2,2],
                            n1 = flex_fit_t@size[1],
                            n2 =flex_fit_t@size[2],
                            start = vec_t$start[t],
                            end = vec_t$end[t])
          
        }
        
        return(mixture)
        
        
      },
      error = function(e) {}
    )
    
  }
  
  estimate <- est |> 
    keep(function(x) length(x) > 1) |> 
    list_rbind() 
  
  
  return(estimate)
  
  
  
})

output_sliding7 <- imap(output_sliding7, ~ .x |> 
                          mutate(sensitivity = .y)) |> 
  list_rbind()


### SLIDING t = 14

output_sliding14 <- map(input, function(x) {
  
  vec_t <- x |> 
    arrange(t_inf) |> 
    transmute(start = t_inf, end = t_inf + 14) |> 
    distinct()
  
  x <- as.data.table(x)
  
  est <- foreach(t = 1:nrow(vec_t), .packages=c('tidyr', 'flexmix', 'countreg', 'purrr', 'readr', 'dplyr', 'fitdistrplus', 'data.table', 'VGAM')) %dopar% {
    
    tryCatch( #error catch
      {
        
        set.seed(as.numeric(vec_t$start[t]))
        
        tmp <- x[t_inf >= vec_t$start[t] & t_inf <= vec_t$end[t]]
        
        
        if (sum(tmp$offspring_count) == 0) {
          
          next
          
        }
        
        ## NB fit with random cluster assignment cannot determine clustering. 
        flex_fit_t <- stepFlexmix(x ~ 1, 
                                  data = data.table(x = tmp$offspring_count),
                                  k = 1:2, 
                                  model = FLXMRnegbin(), nrep = 5, verbose = FALSE)
        
        if (class(flex_fit_t)[1] == "stepFlexmix") {
          
          flex_fit_t <- getModel(flex_fit_t)
          
        }
        
        
        params_t <- parameters(flex_fit_t)
        
        if (flex_fit_t@k0 == 1) {
          
          mixture <- tibble(method = "mixture", 
                            window = "sliding14",
                            r1 = exp(params_t[1,1]),
                            r2 = NA,
                            k1 = params_t[2,1],
                            k2 = NA,
                            n1 = sum(flex_fit_t@cluster),
                            n2 = 0,
                            start = vec_t$start[t],
                            end = vec_t$end[t])
          
          
        } else {
          
          
          mixture <- tibble(method = "mixture", 
                            window = "sliding14",
                            r1 = exp(params_t[1,1]),
                            r2 = exp(params_t[1,2]),
                            k1 = params_t[2,1],
                            k2 = params_t[2,2],
                            n1 = flex_fit_t@size[1],
                            n2 =flex_fit_t@size[2],
                            start = vec_t$start[t],
                            end = vec_t$end[t])
          
        }
        
        return(mixture)
        
        
      },
      error = function(e) {}
    )
    
  }
  
  estimate <- est |> 
    keep(function(x) length(x) > 1) |> 
    list_rbind() 
  
  
  return(estimate)
  
  
  
})

output_sliding14 <- imap(output_sliding14, ~ .x |> 
                           mutate(sensitivity = .y)) |> 
  list_rbind()


bind_rows(output_daily, output_sliding7, output_sliding14) |> 
  saveRDS(file = "data/hk/rds/hk-covid-Rt-kt-mixture.rds")
