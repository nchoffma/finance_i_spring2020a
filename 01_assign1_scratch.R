
# Data and packages -------------------------------------------------------

library(tidyverse)
library(wrapr)
library(stringr)
library(lubridate)
inds <- read.csv("Industry49_data.csv", header = T)

# Question 5 --------------------------------------------------------------

# Clean the excess returns data
ers <- inds %>% 
  mutate(date = paste(year, month, day, sep = "-"),
         date = ymd(date)) %>% 
  select(date, starts_with("eR")) %>% 
  pivot_longer(-date, names_to = "industry", values_to = "eR") %>% 
  mutate(industry = str_remove_all(industry, "eR_"),
         industry = str_remove_all(industry, "I_"))

# a) Expected return by industry
exp_ers <- ers %>% 
  group_by(industry) %>% 
  summarise(
    r_bar = mean(eR, na.rm = T),
    r_bar = r_bar * 12 + 0.01 # convert to ann, add RfR
  )

R_bar <- exp_ers[["r_bar"]] # get vector, which is what we're after
names(R_bar) <- exp_ers[["industry"]] # want to keep track

# b) Variance-covariance matrix
v_mat <- ers %>% 
  pivot_wider(names_from = industry, values_from = eR) %>% 
  select(-date) %>% 
  as.matrix() %>% 
  cov(use = "complete.obs")

v_mat <- v_mat * 12 # convert to annual

# c) Sharpe ratio
sds <- sqrt(diag(v_mat)) # standard deviations
sharpe <- R_bar / sds

# Function to calculate efficient frontier --------------------------------

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

# Question 6 --------------------------------------------------------------

# a)
assets_a <- qc(Util, Txtls)
out1 <- efficient_front(assets = assets_a, risk_free = F)
eff_front <- out1$eff_front
indiv_assets <- out1$indiv_assets

ggplot(eff_front, aes(p_sd, mus)) +
  geom_path() +
  geom_point(data = indiv_assets, aes(var, eR, color = indust)) +
  labs(
    title = "Efficient frontier",
    subtitle = paste0(length(assets_a), " industries, no risk-free asset"),
    color = "Individual industries",
    x = expression(sigma[p]), y = expression(mu[p])
  ) +
  theme(legend.position = "bottom")

# b)
assets_b <- qc("Gold", "Clths", "Rubbr", "Drugs", "Util", "FabPr")
out1 <- efficient_front(assets = assets_b, risk_free = F)
eff_front <- out1$eff_front
indiv_assets <- out1$indiv_assets

ggplot(eff_front, aes(p_sd, mus)) +
  geom_path() +
  geom_point(data = indiv_assets, aes(var, eR, color = indust)) +
  labs(
    title = "Efficient frontier",
    subtitle = paste0(length(assets_b), " industries, no risk-free asset"),
    color = "Individual industries",
    x = expression(sigma[p]), y = expression(mu[p])
  ) +
  theme(legend.position = "bottom")

# c)
out1 <- efficient_front(assets = assets_a, risk_free = T)
eff_front <- out1$eff_front
indiv_assets <- out1$indiv_assets

ggplot(eff_front, aes(p_sd, mus)) +
  geom_path() +
  geom_point(data = indiv_assets, aes(var, eR, color = indust)) +
  labs(
    title = "Efficient frontier",
    subtitle = paste0(length(assets_a), " industries, with risk-free asset"),
    color = "Individual industries",
    x = expression(sigma[p]), y = expression(mu[p])
  ) +
  theme(legend.position = "bottom")

out1 <- efficient_front(assets = assets_b, risk_free = T)
eff_front <- out1$eff_front
indiv_assets <- out1$indiv_assets

ggplot(eff_front, aes(p_sd, mus)) +
  geom_path() +
  geom_point(data = indiv_assets, aes(var, eR, color = indust)) +
  labs(
    title = "Efficient frontier",
    subtitle = paste0(length(assets_b), " industries, with risk-free asset"),
    color = "Individual industries",
    x = expression(sigma[p]), y = expression(mu[p])
  ) +
  theme(legend.position = "bottom")


# Question 7 --------------------------------------------------------------

# Repeat question 6, with randomly selected industries
n_assets_a <- 2
n_assets_b <- 6

# a)
out1 <- efficient_front(N = n_assets_a, risk_free = F)
eff_front <- out1$eff_front
indiv_assets <- out1$indiv_assets

ggplot(eff_front, aes(p_sd, mus)) +
  geom_path() +
  geom_point(data = indiv_assets, aes(var, eR, color = indust)) +
  labs(
    title = "Efficient frontier",
    subtitle = paste0(n_assets_a, " industries, with risk-free asset"),
    color = "Individual industries",
    x = expression(sigma[p]), y = expression(mu[p])
  ) +
  theme(legend.position = "bottom")

# b)
out1 <- efficient_front(N = n_assets_b, risk_free = F)
eff_front <- out1$eff_front
indiv_assets <- out1$indiv_assets

ggplot(eff_front, aes(p_sd, mus)) +
  geom_path() +
  geom_point(data = indiv_assets, aes(var, eR, color = indust)) +
  labs(
    title = "Efficient frontier",
    subtitle = paste0(n_assets_b, " industries, with risk-free asset"),
    color = "Individual industries",
    x = expression(sigma[p]), y = expression(mu[p])
  ) +
  theme(legend.position = "bottom")  

# c)
out1 <- efficient_front(N = n_assets_a, risk_free = T)
eff_front <- out1$eff_front
indiv_assets <- out1$indiv_assets

ggplot(eff_front, aes(p_sd, mus)) +
  geom_path() +
  geom_point(data = indiv_assets, aes(var, eR, color = indust)) +
  labs(
    title = "Efficient frontier",
    subtitle = paste0(n_assets_a, " industries, with risk-free asset"),
    color = "Individual industries",
    x = expression(sigma[p]), y = expression(mu[p])
  ) +
  theme(legend.position = "bottom")

# d)
out1 <- efficient_front(N = n_assets_b, risk_free = T)
eff_front <- out1$eff_front
indiv_assets <- out1$indiv_assets

ggplot(eff_front, aes(p_sd, mus)) +
  geom_path() +
  geom_point(data = indiv_assets, aes(var, eR, color = indust)) +
  labs(
    title = "Efficient frontier",
    subtitle = paste0(n_assets_b, " industries, with risk-free asset"),
    color = "Individual industries",
    x = expression(sigma[p]), y = expression(mu[p])
  ) +
  theme(legend.position = "bottom")

# Question 8 --------------------------------------------------------------

industries <- names(R_bar)[names(R_bar) != "Market"]
industry_pairs <- expand.grid(industries, industries)
names(industry_pairs) <- qc(ind1, ind2)
industry_pairs <- filter(industry_pairs, ind1 != ind2) %>% 
  arrange(ind1)

sharpe_ratio_opt <- function(assets){
  V = v_mat[assets, assets] # variance-covar
  V_inv = solve(V)
  R = R_bar[assets] # expected returns on assets
  one = rep(1.0, N)
  
  w_star_sharpe = 
}
