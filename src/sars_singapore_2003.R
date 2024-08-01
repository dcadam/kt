library(readr)
library(data.table)
library(fitdistrplus)
library(insight)

set.seed(1234)

# Percentile bias correction function
bias_correction <- function(vec, value_mle) {
  
  tmp <- data.table(vec = vec) |> 
    mutate(bca = case_when(vec > value_mle ~ 1,
                           TRUE ~ 0))
  
  # 90% CI
  quantile(tmp$vec,
           c(pnorm((2*qnorm(mean(tmp$bca))) - 1.65),
             pnorm((2*qnorm(mean(tmp$bca))) + 1.65)))
  
}

# Singapore (Lloyd-Smith et al JLS) 
# Original tree from Leo, Y. S. et al. Severe acute respiratory syndrome - Singapore, 2003. Morbid. Mortal. Wkly. Rep. 52, 405-411 (2003).

sars_singapore_2003 <- read_rds(file = "data/sars_singapore_2003.rds")

## Generations 1 to 3 only
# Estimate R and k
fit <- fitdist(data = sars_singapore_2003$z, distr = 'nbinom') # k = 0.16

# non-parametric bootsrap (higher number of replicates 50,000 vs 10,000 in JLS to allow for differences in seed but aiming to match JLS result
fit_boot <- bootdist(fit, bootmethod = 'nonparam', niter = 50000, parallel = 'multicore', ncpus = 4)

# Uncorrected 90% CI for k = 0.10, 0.36 matches result from JLS
quantile(fit_boot$estim$size, probs = c(0.05, 0.95)) |> 
  round(digits = 2)

# Bias-corrected 90% CI for k = 0.11, 0.64 matches result from JLS
ci_k <- bias_correction(vec = fit_boot$estim$size, value_mle = fit$estimate[1])
ci_R <- bias_correction(vec = fit_boot$estim$mu, value_mle = fit$estimate[2])


sars_before_control <- data.table(
  Z = "Before control",
  R = fit$estimate[2],
  k = fit$estimate[1]) |> 
  pivot_longer(cols = c(R, k)) |> 
  mutate(CI_low = c(ci_R[1],
                    ci_k[1]),
         CI_high = c(ci_R[2],
                      ci_k[2])
  )


## All generations 

## Generations 1 to 3 only
# Estimate R and k
fit_all <- fitdist(data = sars_singapore_2003$z_all, distr = 'nbinom') # k = 0.16

# non-parametric bootsrap (higher number of replicates 50,000 vs 10,000 in JLS to allow for differences in seed but aiming to match JLS result
fit_boot_all <- bootdist(fit_all, bootmethod = 'nonparam', niter = 50000, parallel = 'multicore', ncpus = 4)

# Uncorrected 90% CI for k = 0.10, 0.36 matches result from JLS
quantile(fit_boot_all$estim$size, probs = c(0.05, 0.95)) |> 
  round(digits = 2)

# Bias-corrected 90% CI for k = 0.11, 0.64 matches result from JLS
ci_all_k <- bias_correction(vec = fit_boot_all$estim$size, value_mle = fit_all$estimate[1])
ci_all_R <- bias_correction(vec = fit_boot_all$estim$mu, value_mle = fit_all$estimate[2])


sars_all <- data.table(
  Z = "Complete offspring",
  R = fit_all$estimate[2],
  k = fit_all$estimate[1]) |> 
  pivot_longer(cols = c(R, k)) |> 
  mutate(CI_low = c(ci_all_R[1],
                    ci_all_k[1]),
         CI_high = c(ci_all_R[2],
                     ci_all_k[2])
  )


bind_rows(sars_all, 
          sars_before_control) |> 
  format_table(ci_brackets = c("(",")"), digits = 2, ci_digits = 2) |> 
  unite("res", value:CI, sep = " ") |> 
  pivot_wider(names_from = name, values_from = res) |> 
  arrange(Z) |> 
  write_csv(file = "output/sars_singapore_2003.csv")



