library(purrr)
library(data.table)

### Generate representative data for Fig.4d and 4e (replicated pseudo branching process)
set.seed(12345)

epidemic_size <- 100000
R1 <- 2
R2 <- 0.5
k <- c(0.1, 0.5, 1, 10, 100)

relative_offspring <- map(k, function(x) {
  
  
  R1_size <- epidemic_size*0.333
  R2_size <- epidemic_size*0.666
  
  uncontrolled_cases <- rnbinom(n = R1_size, size = x, mu = R1)
  controlled_cases <- rnbinom(n = R2_size, size = x, mu = R2)
  
  data.table(k = x, 
             Z = c(uncontrolled_cases, controlled_cases), 
             source = c(rep("Uncontrolled", times = R1_size),
                        rep("Controlled", times = R2_size)))
  
  
})