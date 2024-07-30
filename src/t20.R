#### FUNCTION to calculate the proportion of transmission expected from the upper 20% most infectious cases (t_20)

### Defines the gamma density function in terms of the mean reproductive number (R) and the dispersion paramter (k)
dgammaRk <- function(x, R, k) {
  out <- dgamma(x, shape = k, scale = R / k)
  return(out)
}

### Defines the gamma distribution function in terms of the mean reproductive number (R) and the dispersion paramter (k)
pgammaRk <- function(x, R, k) {
  out <- pgamma(x, shape = k, scale = R / k)
  return(out)
}

### Generates a log scaled sequence of real numbers (useful for plotting later)
lseq <- function(from, to, length.out) {
  exp(seq(log(from), log(to), length.out = length.out))
}

### Gamma probability densitiy function (pdf) describing the individual reproductive number ν given R0 and k 
fvx <- function(x, R, k) {
  dgammaRk(x = x, R = R, k = k)
}

### Unitroot solver for the solution to u, or the value of x from the gamma CDF given the desired proportion of transmission.
solve_for_u <- function(prop, R, k) {
  
  f <- function(x, prop) {
    res <- 1 - pgammaRk(x, R, k)
    res - prop
  }
  
  # Initial interval for u
  lower <- 0
  upper <- 1
  
  # Find the root using uniroot
  root <- uniroot(f, prop = prop, interval = c(lower, upper), extendInt="yes")
  
  root$root
  
}

### Function: Calculates the expected proportion of transmission from the upper 'prop'% of most infectious cases
proportion_transmission2 <- function(R, k, prop) {
  
  u <- solve_for_u(prop = prop, R = R, k = k)
  
  integral_result <- integrate(function(u) u * fvx(u, R, k), lower = 0, upper = u)
  
  res <- integral_result$value / R
  
  1 - res
  
}

