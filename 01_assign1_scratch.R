
# Data and packages -------------------------------------------------------

library(tidyverse)
library(wrapr)
library(stringr)
library(lubridate)
library(stargazer)
inds <- read.csv("Industry49_data.csv", header = T)
to_show <- scan("assignment_writeups/01_assign/inds_to_show.txt", 
                what = "character")

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
sharpe <-( R_bar - 0.01) / sds # should 

returns_data <- cbind.data.frame(
  Industry = names(R_bar), 
  excess_r = R_bar * 100, 
  var = diag(v_mat) * 100,
  stand_dev = sds * 100, 
  sharpe = sharpe * 100
)

to_show <- to_show
return_table_tex <- stargazer(filter(returns_data, Industry %in% to_show), 
                              summary = F, rownames = F)
write(return_table_tex, file = "assignment_writeups/01_assign/table_1_returns.txt")

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

figpath <- "assignment_writeups/01_assign/"

# a)
assets_a <- qc(Util, Agric)
out1 <- efficient_front(assets = assets_a, risk_free = F)
eff_front <- out1$eff_front
indiv_assets <- out1$indiv_assets

theme_set(theme_bw())
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
ggsave(paste0(figpath, "plot_6a.eps"))

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
ggsave(paste0(figpath, "plot_6b.eps"))

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
ggsave(paste0(figpath, "plot_6c.eps"))


out1 <- efficient_front(assets = assets_b, risk_free = T)
eff_front <- out1$eff_front
indiv_assets <- out1$indiv_assets

# example_weights <- stargazer(eff_front, summary = F, rownames = F)
# write(example_weights, file = "assignment_writeups/01_assign/table_3_weights.txt")

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
ggsave(paste0(figpath, "plot_6d.eps"))

# Question 8 --------------------------------------------------------------

# Moved before q7, to avoid any confusion with what Rbar is

# Get all combinations of industries
industries <- names(R_bar)[names(R_bar) != "Market"]
industry_pairs <- t(combn(industries, m = 2))
colnames(industry_pairs) <- qc(ind1, ind2)

# Function to calculate portfolio weights for maximal Sharpe Ratio
sharpe_ratio_opt <- function(assets){
  N = length(assets)
  V = v_mat[assets, assets] # variance-covar
  V_inv = solve(V)
  R = R_bar[assets] # expected returns on assets
  one = rep(1.0, N)
  Rf = 0.01
  Rf_one = one * Rf
  
  w_star_sharpe = (V_inv %*% (R - Rf_one)) * 
    as.numeric(1 / t(one) %*% V_inv %*% (R - Rf_one))
}

# How many involve no short selling?
npair <- nrow(industry_pairs)
opt_sharps <- 0L
for (i in seq_len(npair)){
  w_star_sharpe <- sharpe_ratio_opt(industry_pairs[i, ])
  opt_sharps <- opt_sharps + all(w_star_sharpe >= 0)
}

sharpe_prob <- opt_sharps / npair # ~1/2

# Question 7 --------------------------------------------------------------

# Repeat question 6, with different subsets of data
# Breaking rules here with copy/paste, but so be it
# Filtering post-2000

# a) Expected return by industry
exp_ers <- ers %>% 
  filter(year(date) >= 2000) %>% 
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
sharpe <-( R_bar - 0.01) / sds # should 

returns_data <- cbind.data.frame(
  Industry = names(R_bar), 
  excess_r = R_bar * 100, 
  var = diag(v_mat) * 100,
  stand_dev = sds * 100, 
  sharpe = sharpe * 100
)

to_show <- to_show
return_table_tex <- stargazer(filter(returns_data, Industry %in% to_show), 
                              summary = F, rownames = F)
write(return_table_tex, file = "assignment_writeups/01_assign/table_2_returns.txt")

# a)
assets_a <- qc(Util, Agric)
out1 <- efficient_front(assets = assets_a, risk_free = F)
eff_front <- out1$eff_front
indiv_assets <- out1$indiv_assets

ggplot(eff_front, aes(p_sd, mus)) +
  geom_path() +
  geom_point(data = indiv_assets, aes(var, eR, color = indust)) +
  labs(
    title = "Efficient frontier, post-2000",
    subtitle = paste0(length(assets_a), " industries, no risk-free asset"),
    color = "Individual industries",
    x = expression(sigma[p]), y = expression(mu[p])
  ) +
  theme(legend.position = "bottom")
ggsave(paste0(figpath, "plot_7a.eps"))

# b)
assets_b <- qc("Gold", "Clths", "Rubbr", "Drugs", "Util", "FabPr")
out1 <- efficient_front(assets = assets_b, risk_free = F)
eff_front <- out1$eff_front
indiv_assets <- out1$indiv_assets

ggplot(eff_front, aes(p_sd, mus)) +
  geom_path() +
  geom_point(data = indiv_assets, aes(var, eR, color = indust)) +
  labs(
    title = "Efficient frontier, post-2000",
    subtitle = paste0(length(assets_b), " industries, no risk-free asset"),
    color = "Individual industries",
    x = expression(sigma[p]), y = expression(mu[p])
  ) +
  theme(legend.position = "bottom")
ggsave(paste0(figpath, "plot_7b.eps"))

# c)
out1 <- efficient_front(assets = assets_a, risk_free = T)
eff_front <- out1$eff_front
indiv_assets <- out1$indiv_assets

ggplot(eff_front, aes(p_sd, mus)) +
  geom_path() +
  geom_point(data = indiv_assets, aes(var, eR, color = indust)) +
  labs(
    title = "Efficient frontier, post-2000",
    subtitle = paste0(length(assets_a), " industries, with risk-free asset"),
    color = "Individual industries",
    x = expression(sigma[p]), y = expression(mu[p])
  ) +
  theme(legend.position = "bottom")
ggsave(paste0(figpath, "plot_7c.eps"))


out1 <- efficient_front(assets = assets_b, risk_free = T)
eff_front <- out1$eff_front
indiv_assets <- out1$indiv_assets

ggplot(eff_front, aes(p_sd, mus)) +
  geom_path() +
  geom_point(data = indiv_assets, aes(var, eR, color = indust)) +
  labs(
    title = "Efficient frontier, post-2000",
    subtitle = paste0(length(assets_b), " industries, with risk-free asset"),
    color = "Individual industries",
    x = expression(sigma[p]), y = expression(mu[p])
  ) +
  theme(legend.position = "bottom")
ggsave(paste0(figpath, "plot_7d.eps"))

