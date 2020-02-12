
# Setup -------------------------------------------------------------------

library(tidyverse)
library(pracma)


# Question 2: Test CAPM ---------------------------------------------------


# Question 3: Plot Arbitrage ----------------------------------------------


# Question 5: Calibrating -------------------------------------------------

mean_var <- function(d, mu, sig){
  # Gets mean and variance of returns d (equal probabilities)
  # Compares them to the desired mu and sigma
  # d should be 1x2
  
  # Mean and variance for this guess
  mu_g = sum(d) / 2
  sig_g = (1 / 2) * ((d[1] - mu_g) ^ 2) + 
    (1 / 2) * ((d[2] - mu_g) ^ 2)
  
  dif = max(abs(mu_g - mu), abs(sig_g - sig))
  return(dif)
}

out <- optim(
  par = c(0.9, 1.1), # initial guess
  fn = mean_var, 
  mu = 1.07, sig = 0.04, # desired values
  method = "L-BFGS-B", # allows for constraints
  
  lower = c(0.0, 1.0), # lower bounds on d
  upper = c(1.0, Inf) # upper bounds on d
)

d_opt <- out$par

d_mat <- cbind(rep(1.01, 2), d_opt) # matrix of returns
prices <- rep(1.0, 2)
a_mat <- prices %*% solve(d_mat) # A-D State prices
rf <- 1 / sum(a_mat) # 1.01, as desired
qs <- a_mat / sum(a_mat) # risk-neutral probabilities
ms <- a_mat / 0.5 # pricing kernel 

# Question 6: HJ Bounds ---------------------------------------------------


