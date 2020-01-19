
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

# Question 6 --------------------------------------------------------------

# Initialize
N <- 6 # number of industries
M <- 30 # number of mu values to show on efficient frontier
industries <- names(R_bar)[names(R_bar) != "Market"]
assets <- sample(industries, N)

# Get the necessary matrices/vectors
V <- v_mat[assets, assets] # variance-covar
V_inv <- solve(V)
R <- R_bar[assets] # expected returns 
one <- rep(1.0, N)

# Get the optimal portfolios
w_r <- V_inv %*% R / as.numeric(one %*% V_inv %*% R)
w_1 <- V_inv %*% one / as.numeric(one %*% V_inv %*% one)

# Calculate alpha, as function of mu
get_alpha <- function(mu){
  alpha = as.numeric(mu - t(w_1) %*% R) / 
    as.numeric(t(w_r) %*% R - t(w_1) %*% R)
}

get_alpha <- Vectorize(get_alpha)

# Calculate efficient frontier: optimal portfolios at each return
mus <- seq(0, 0.20, length.out = M) # return values 
alphas <- get_alpha(mus)

wts_1 <- alphas %*% t(w_r)
wts_2 <- (1 - alphas) %*% t(w_1)
wts_opt <- wts_1 + wts_2

eff_front <- cbind.data.frame(mus, wts_opt)
port_vars <- vector("numeric")
for (i in 1:M){
  v = t(wts_opt[i, ]) %*% V %*% wts_opt[i, ]
  port_vars <- append(port_vars, v)
}

eff_front[["p_sd"]] <- sqrt(port_vars)
indiv_assets <- data.frame(
  indust = assets,
  eR = R,
  var = sqrt(diag(V))
)

theme_set(theme_bw())
ggplot(eff_front, aes(p_sd, mus)) + 
  geom_point() + 
  geom_path() + 
  geom_point(data = indiv_assets, aes(var, eR, color = indust)) + 
  labs(
    title = "Efficient frontier",
    subtitle = paste0(N, " industries"),
    color = "Individual industries",
    x = expression(sigma), y = expression(mu)
  ) + 
  theme(legend.position = "bottom")

# TODO: wrap the above in a function, so that it can quickly generate plots
# TODO: add capability to allow for a risk-free asset

# Question 7 --------------------------------------------------------------

# Just replicating 5 and 6 for various subsets

# Question 8 --------------------------------------------------------------

industry_pairs <- expand.grid(industries, industries)
names(industry_pairs) <- qc(ind1, ind2)
industry_pairs <- filter(industry_pairs, ind1 != ind2) %>% 
  arrange(ind1)
