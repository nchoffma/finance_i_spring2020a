# Archive -----------------------------------------------------------------


# Initialize
N <- 7 # number of industries
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
    x = expression(sigma[p]), y = expression(mu[p])
  ) + 
  theme(legend.position = "bottom")


# add risk-free asset
# Initialize
N <- 4 # number of industries
M <- 30 # number of mu values to show on efficient frontier
industries <- names(R_bar)[names(R_bar) != "Market"]
assets <- sample(industries, N)

# Get the necessary matrices/vectors
V <- v_mat[assets, assets] # variance-covar
V_inv <- solve(V)
R <- R_bar[assets] # expected returns 
one <- rep(1.0, N)
Rf <- 0.01
Rf_1 <- one * Rf

opt_w_rf <- function(mu){
  w_opt = as.numeric((mu - Rf) / (t(R - Rf_1) %*% V_inv %*% (R - Rf_1))) * 
    V_inv %*% (R - Rf_1)
}

mus <- seq(-0.05, 0.20, length.out = M) # return values 
eff_front <- t(sapply(mus, opt_w_rf)) # get optimal w
colnames(eff_front) <- assets
port_vars <- vector("numeric")
for (i in 1:M){
  v = t(eff_front[i, ]) %*% V %*% eff_front[i, ]
  port_vars <- append(port_vars, v)
}

eff_front <- cbind.data.frame(mus, eff_front)
eff_front[["p_sd"]] <- sqrt(port_vars)
indiv_assets <- data.frame(
  indust = assets,
  eR = R,
  var = sqrt(diag(V))
)

ggplot(eff_front, aes(p_sd, mus)) + 
  geom_point() + 
  geom_path() + 
  geom_point(data = indiv_assets, aes(var, eR, color = indust)) + 
  labs(
    title = "Efficient frontier",
    subtitle = paste0(N, " industries"),
    color = "Individual industries",
    x = expression(sigma[p]), y = expression(mu[p])
  ) + 
  theme(legend.position = "bottom")