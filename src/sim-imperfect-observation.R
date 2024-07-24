### SIMULATING a PSUEDO BRANCHING PROCESS EPIDEMIC with constant k, a fixed epidemic size, and incomplete observation

library(tidyverse)
library(fitdistrplus)
library(data.table)
library(doSNOW)

reps <- 300
R_change <- seq(0.5, 3, by = 0.5)
k_change <- c(0.1, 0.5, 1, 10, 100)
nn <- 1000

obs_count <- c(30, 100, 300, 500, 1000) ## primary cases under observation
ct_prop <- seq(0.25, 1, by = 0.25) ## secondary cases under observation

C1 <- crossing(R1 = R_change, k1 = k_change)
C2 <- crossing(R2 = R_change, k2 = k_change)

rk_steps <- crossing(C1, C2, ct_prop, obs_count) |> 
  as.data.table() |> 
  filter(k1 == k2)

ncpus <- parallel::detectCores()

cluster <- makeCluster(ncpus)
registerDoParallel(cluster)

## SIMULATE

set.seed(123)
out <- vector(mode = 'list', length = reps)
for (i in 1:length(out)) {
  
  pb <- txtProgressBar(min = 1, max = length(out_Rk), style = 3)
  
  tmp_Rk <- foreach(j = 1:nrow(rk_steps), .packages=c('dplyr', 'fitdistrplus', 'data.table')) %dopar% {
    
    repeat {
      
      nb1 <- rnbinom(n = nn/2, mu = rk_steps$R1[[j]], size = rk_steps$k1[[j]])
      nb1 <- map_dbl(nb1, function(z) rbinom(n = 1, size = z, prob = rk_steps$ct_prop[[j]]))
      
      nb2 <- rnbinom(n = nn/2, mu = rk_steps$R2[[j]], size = rk_steps$k2[[j]])
      nb2 <- map_dbl(nb2, function(z) rbinom(n = 1, size = z, prob = rk_steps$ct_prop[[j]]))
      
      nb_mix <- sample(x = c(nb1,nb2), replace = FALSE, size = rk_steps$obs_count[[j]])
      
      nb_mix <- data.table(x = nb_mix)
      
      tryCatch( #error catch 
        {
          
          mix_fit <- fitdist(data = nb_mix$x, distr = 'nbinom', keepdata = FALSE)
          clear1 <- TRUE
          
          vx_mix <- var(nb_mix$x)
          
          
        },
        error = function(e) {clear1 <<- FALSE}
      )
      
      if (clear1 == TRUE) {
        
        tmp_Rk_j <- rk_steps[j] |> 
          mutate(
            exp_R = ((R1+R2)/2)*ct_prop,
            exp_k = 1,
            exp_var = nb_var(exp_R, exp_k),
            fit_var = vx_mix,
            fit_mu = mix_fit$estimate[2],
            fit_k = mix_fit$estimate[1],
            fit_aic = mix_fit$aic)
        
        rm(nb12, nb1, nb2, mix_fit)
        
        
        break ## breaks the repeat only when estimation is complete
      }
      
    }
    
    return(tmp_Rk_j)
    
  }
  
  tmp_Rk <- rbindlist(tmp_Rk) |> 
    mutate(reps = i)
  
  setTxtProgressBar(pb, i)
  
  out[[i]] <- tmp_Rk
  
  
}
