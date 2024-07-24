library(data.table)
library(fitdistrplus)
library(foreach)
library(doParallel)

## nb variance function 
nb_var <- function(mean, dispersion) {
  
  mean + mean^2 / dispersion #which equals variance
  
  
}


reps <- 100
set.seed(123)
R_change <- seq(0.5, 3, by = 0.5)
k_change <- 1
nn <- 1000

C1 <- crossing(R1 = R_change, k1 = k_change)
C2 <- crossing(R2 = R_change, k2 = k_change)

rk_steps <- crossing(C1, C2) |> 
  as.data.table()

cluster <- makeCluster(15) ## start parralel cluster within each epidemic j to control memory
registerDoParallel(cluster) 

out <- vector(mode = 'list', length = reps)
for (i in 1:length(out)) {
  
  pb <- txtProgressBar(min = 1, max = length(out), style = 3)
  
  tmp_Rk <- foreach(j = 1:nrow(rk_steps), .packages=c('dplyr', 'fitdistrplus', 'data.table')) %dopar% {
    
    
    repeat {
      
      nb1 <- rnbinom(n = nn/2, mu = rk_steps$R1[[j]], size = rk_steps$k1[[j]])
      nb2 <- rnbinom(n = nn/2, mu = rk_steps$R2[[j]], size = rk_steps$k2[[j]])
      
      nb12 <- sample(x = c(nb1, nb2), replace = FALSE, size = nn)
      
      offspring_x <- data.table(x = nb12)
      
      
      tryCatch( #error catch 1
        {
          
          single_fit <- fitdist(data = offspring_x$x, distr = 'nbinom', keepdata = FALSE)
          clear1 <- TRUE
          
          vx <- var(offspring_x$x)
          
          
        },
        error = function(e) {clear1 <<- FALSE}
      )
      
      if (clear1 == TRUE) {
        
        tmp_Rk_j <-rk_steps[j] |> 
          mutate(
            exp_R = (R1+R2)/2,
            exp_k = 1,
            exp_var = nb_var(exp_R, exp_k),
            x1_var = var(nb1),
            x2_var = var(nb2),
            fit_mu = single_fit$estimate[2],
            fit_k = single_fit$estimate[1],
            fit_var = vx,
            fit_aic = single_fit$aic)
        
        rm(single_fit)
        
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

saveRDS(object = out, file = "code/simulations/data/underlying_variance.rds")

