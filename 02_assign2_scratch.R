
# Setup -------------------------------------------------------------------

library(tidyverse)
library(pracma)
library(wrapr)
library(stringr)
library(lubridate)
library(stargazer)
inds <- read.csv("Industry49_data.csv", header = T)
theme_set(theme_bw())

# Clean Data --------------------------------------------------------------

# Get expected returns
ers <- inds %>% 
  mutate(date = paste(year, month, day, sep = "-"),
         date = ymd(date)) %>% 
  select(date, starts_with("eR"))

exp_ers <- ers %>% 
  pivot_longer(-date, names_to = "industry", values_to = "eR") %>% 
  mutate(industry = str_remove_all(industry, "eR_"),
         industry = str_remove_all(industry, "I_")) %>% 
  group_by(industry) %>% 
  summarise(
    r_bar = mean(eR, na.rm = T),
    r_bar = r_bar * 12 + 0.01 # convert to ann, add RfR
  )

# Inputs to efficient frontier function
R_bar <- exp_ers[["r_bar"]] # get vector, which is what we're after
names(R_bar) <- exp_ers[["industry"]] # want to keep track

names(ers) <- str_remove_all(names(ers), "eR_")
names(ers) <- str_remove_all(names(ers), "I_")
v_mat <- ers %>% 
  select(-date) %>% 
  as.matrix() %>% 
  cov(use = "complete.obs")

v_mat <- v_mat * 12 # convert to annual


# Functions ---------------------------------------------------------------

# Efficient frontier
efficient_front <- function(M = 30, assets = c(), N = 0,
                            risk_free = F){
  
  # This takes either a vector of assets, or a number (N) 
  # of assets to choose at random
  
  # Get number of assets if supplied, asset names if not
  if (length(assets) == 0){
    industries = names(R_bar)[names(R_bar) != "Market"]
    assets = sample(industries, N)
  } else {
    N = length(assets)
  }
  
  # Get the necessary matrices/vectors
  mus = seq(-0.05, 0.20, length.out = M) # return values 
  V = v_mat[assets, assets] # variance-covar
  V_inv = solve(V)
  R = R_bar[assets] # expected returns on assets
  one = rep(1.0, N)
  
  # Calculate efficient frontier portfolio weights
  if (risk_free){
    Rf = 0.01 # default risk-free rate is 1%
    Rf_1 = one * Rf
    
    opt_w_rf <- function(mu){
      w_opt = as.numeric((mu - Rf) / (t(R - Rf_1) %*% V_inv %*% (R - Rf_1))) * 
        V_inv %*% (R - Rf_1)
    }
    
    eff_front = t(sapply(mus, opt_w_rf)) # get optimal w
    colnames(eff_front) = assets
    
  } else {
    # Get the optimal portfolios
    w_r = V_inv %*% R / as.numeric(one %*% V_inv %*% R)
    w_1 = V_inv %*% one / as.numeric(one %*% V_inv %*% one)
    
    # Calculate alpha, as function of mu
    get_alpha = function(mu){
      alpha = as.numeric(mu - t(w_1) %*% R) / 
        as.numeric(t(w_r) %*% R - t(w_1) %*% R)
    }
    
    get_alpha = Vectorize(get_alpha)
    
    # Calculate efficient frontier: optimal portfolios at each return
    alphas = get_alpha(mus)
    
    wts_1 = alphas %*% t(w_r)
    wts_2 = (1 - alphas) %*% t(w_1)
    eff_front = wts_1 + wts_2
  }
  
  port_vars = vector("numeric")
  for (i in 1:M){
    v = t(eff_front[i, ]) %*% V %*% eff_front[i, ]
    port_vars = append(port_vars, v)
  }
  
  eff_front = cbind.data.frame(mus, eff_front)
  eff_front[["p_sd"]] = sqrt(port_vars)
  
  indiv_assets = data.frame(
    indust = assets,
    eR = R, 
    var = sqrt(diag(V))
  )
  
  # Return args
  return(list(eff_front = eff_front, indiv_assets = indiv_assets))
  
}

# Optimal Sharpe Ratio
sharpe_ratio_opt <- function(assets, Rf = 0.01){
  N = length(assets)
  V = v_mat[assets, assets] # variance-covar
  V_inv = solve(V)
  R = R_bar[assets] # expected returns on assets
  one = rep(1.0, N)
  Rf_one = one * Rf
  
  w_star_sharpe = (V_inv %*% (R - Rf_one)) * 
    as.numeric(1 / t(one) %*% V_inv %*% (R - Rf_one))
}

# Question 2: Test CAPM ---------------------------------------------------

# Run multiple regression
y <- select(ers, -c(date, eR_Market)) %>% as.matrix()
rm <- ers["eR_Market"] %>% as.matrix()
models <- lm(formula = y ~ rm)

# Get coefficient values
coeffs <- data.frame(t(models$coefficients))
names(coeffs) <- qc(alpha, beta)
coeffs["industry"] <- rownames(coeffs)
rownames(coeffs) <- c()

# Get p-values
all_coefs <- coef(summary(models))
all_pvals <- data.frame(t(sapply(all_coefs, function(m) m[, 4])))
names(all_pvals) <- qc(p_alpha, p_beta)
all_pvals["industry"] <- str_remove_all(rownames(all_pvals), "Response ")
rownames(all_pvals) <- c()

# merge results
capm_results <- merge(coeffs, all_pvals, by = "industry")
capm_results <- mutate(capm_results, 
                       industry = str_remove_all(industry, "eR_I_"))

# export results
capm_results_out <- capm_results %>% 
  select(industry, beta, alpha, p_alpha) %>% 
  stargazer(., summary = F)

write(capm_results_out, file = "assignment_writeups/02_assign/table_1_capm.txt")

# One issue: do we have autocorrelation here? may need to try something 
# different (MLE, Sandwich, etc.) to get the right standard errors.
# The estimates of alpha should be consistent, though.

# B 
capm_results %>% 
  merge(., r_bar, by = "industry") %>% 
  ggplot(aes(beta, r_bar * 100)) + 
  geom_point() + 
  labs(
    title = "Fama-Macbeth Test",
    x = expression(beta),
    y = expression(r[i]*"(%)")
  ) + 
  geom_smooth(method = "lm", se = T, linetype = "dashed", color = "black",
              size = 0.5)

# Question 3: Plot Arbitrage ----------------------------------------------

P <- c(1.575, 1.35, 3.425)
D <- matrix(c(2, 1, 4, 1, 3, 3), nrow = 2, byrow = T)

# Save this for last, per his instructions

# Question 5: Calibrating -------------------------------------------------

pricing_kernel = function(mu, sig, prices = rep(1.0, 2)){
  # Function to get the following, given mu and sigma:
  #   - d_mat: matrix of returns
  #   - a_mat: A-D state prices
  #   - qs: risk-neutral probabilities
  #   - ms: pricing kernel
  # assumes price vector of ones
  # mu is total return (e.g. 1 + mu + rf)
  
  mean_var = function(d, mu, sig){
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
  
  out = optim(
    par = c(0.9, 1.1), # initial guess
    fn = mean_var, 
    mu = mu, sig = sig, # desired values
    method = "L-BFGS-B", # allows for constraints
    
    lower = c(0.0, 1.0), # lower bounds on d
    upper = c(1.0, Inf) # upper bounds on d
  )
  
  d_opt = out$par
  
  d_mat = cbind(rep(1.01, 2), d_opt) # matrix of returns
  a_mat = prices %*% solve(d_mat) # A-D State prices
  rf = 1 / sum(a_mat) # 1.01, as desired
  qs = a_mat / sum(a_mat) # risk-neutral probabilities
  ms = a_mat / 0.5 # pricing kernel 
  
  return(list(
    d_mat = d_mat, 
    a_mat = a_mat, 
    rf_check = rf, 
    qs = qs, ms = ms
  ))
}

baseline_result <- pricing_kernel(1.07, 0.04)

# f) Plotting
mus <- seq(1.01, 1.11, length.out = 10)
m_grid <- t(sapply(mus, function(mu) pricing_kernel(mu, sig = 0.04)$ms))
m_ratio <- data.frame(
  mu = mus,
  m_rat = m_grid[, 1] / m_grid[, 2]
)
ggplot(m_ratio, aes(mu, m_rat)) + 
  geom_point() + 
  geom_line()

sig_m <- data.frame(
  mu = mus,
  sigm = 0.5 * (m_grid[, 1] - mus) ^ 2 + 
    0.5 * (m_grid[, 2] - mus) ^ 2
)

ggplot(sig_m, aes(mu, sigm)) + 
  geom_point() + 
  geom_line()

sigs <- seq(0.01, 0.10, length.out = 10)
m_grid_sig <- t(sapply(sigs, function(sig) pricing_kernel(mu = 1.07, sig)$ms))
m_ratio_sig <- data.frame(
  sig = sigs,
  m_rat = m_grid_sig[, 1] / m_grid_sig[, 2]
)

ggplot(m_ratio_sig, aes(sig, m_rat)) + 
  geom_point() + 
  geom_line()

sig_m_sig <- data.frame(
  mu = mus,
  sigm = 0.5 * (m_grid_sig[, 1] - mus) ^ 2 + 
    0.5 * (m_grid_sig[, 2] - mus) ^ 2
)

ggplot(sig_m_sig, aes(mu, sigm)) + 
  geom_point() + 
  geom_line()

# Question 6: HJ Bounds ---------------------------------------------------

# a)
assets <- c("BusSv", "Agric") # looks fine
out <- efficient_front(assets = assets, risk_free = F)
eff_front <- out$eff_front
indiv_assets <- out$indiv_assets

ggplot(eff_front, aes(p_sd, mus)) +
  geom_path() +
  geom_point(data = indiv_assets, aes(var, eR, color = indust)) +
  labs(
    title = "Efficient frontier",
    subtitle = paste0(length(assets), " industries, no risk-free asset"),
    color = "Individual industries",
    x = expression(sigma[p]), y = expression(mu[p])
  ) +
  theme(legend.position = "bottom")

# b) 
rfs <- seq(1.01, 1.05, by = 0.01) - 1

opt_sharpes <- t(sapply(rfs, function(rf) sharpe_ratio_opt(assets, Rf = rf)))

get_sharpe_ratio <- function(assets, w, Rf = 0.01){
  # given assets and portfolio weights in w, calculate the Sharpe ratio
  R = R_bar[assets]
  V = v_mat[assets, assets]
  
  eer = w %*% R
  sdr = sqrt(t(w) %*% V %*% w)
  sharpe = as.numeric((eer - Rf) / sdr)
  return(sharpe)
}

sharpes <- vector("numeric", nrow(opt_sharpes))
for (i in 1:nrow(opt_sharpes)){
  sharpes[i] <- get_sharpe_ratio(assets, opt_sharpes[i, ], Rf = rfs[i])
}

hjb <- data.frame(
  Rf = rfs,
  sharpe = sharpes
)

hjb <- mutate(hjb,
              bound = sharpe / (1 + Rf),
              Em = 1 / (1 + Rf))

ggplot(hjb, aes(Em, bound)) + 
  geom_point() + 
  geom_line()

