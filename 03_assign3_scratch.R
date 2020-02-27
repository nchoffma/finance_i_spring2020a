

# Question 1' -------------------------------------------------------------

obj <- function(a){
  ((10000 ^ a) / a) - (1 / 2) * ((9900 ^ a) / a) - 
    (1 / 2) * ((10110 ^ a) / a)
}

curve(obj, from = -10, to = 0)

avals <- seq(-10, 0, length.out = 101)
guess <- data.frame(
  a = avals,
  v = obj(avals)
)

r <- uniroot(obj, c(-9,-8))
alpha_max <- r$root

obj2 <- function(G){
  a = alpha_max - 0.01 # to ensure refusing previous lottery
  10000 ^ a - (1 / 2) * (9000 ^ a) - (1 / 2) * ((10000 + G) ^ a)
}

curve(obj2, -1000, 1000)
abline(h = 0) # seems to asypmtotically approach 0

# Question 2 --------------------------------------------------------------

# Payoffs and prices
P_up <- matrix(c(1, 1.9444), ncol = 1)
D_up <- matrix(c(1, 2.5, 1, 1.5), nrow = 2, byrow = T)
P_dn <- matrix(c(0.96, 1.2742), ncol = 1)
D_dn <- matrix(c(1, 2, 1, 1), nrow = 2, byrow = T)
P_0 <- matrix(c(0.9487, 1.5285), ncol = 1)
D_0 <- rbind(t(P_up), t(P_dn))

# Cond'l A-D state prices
a_zuu <- t(P_up) %*% solve(D_up) %*% c(1, 0)
a_zud <- t(P_up) %*% solve(D_up) %*% c(0, 1)
a_zdu <- t(P_dn) %*% solve(D_dn) %*% c(1, 0)
a_zdd <- t(P_dn) %*% solve(D_dn) %*% c(0, 1)

# Date-0 prices
a0_zuu <- t(P_0) %*% solve(D_0) %*% c(a_zuu, 0)
a0_zud <- t(P_0) %*% solve(D_0) %*% c(0, a_zud)
a0_zdu <- t(P_0) %*% solve(D_0) %*% c(a_zdu, 0)
a0_zdd <- t(P_0) %*% solve(D_0) %*% c(0, a_zdd)

# Prices for Z_1 = {z_u = {uu, ud}, z_d = {du, dd}} 
a_z1 <- t(P_0) %*% diag(2) %*% solve(D_0)

# SDFs
a0 <- rbind(a0_zuu, a0_zud, a0_zdu, a0_zdd) # unconditional AD prices
a1 <- rbind(a_zuu, a_zud, a_zdu, a_zdd) # conditional (t = 1) AD prices
M_2 <- a1 * 2
M_1 <- t(a_z1) * 2

# Pricing kernels
m_2 <- M_2 / rep(M_1, each = 2)
m_1 <- M_1

# Risk-free rate
Rf_0 <- 1 / sum(a_z1) # ~ 3.1%
Rf_1_zu <- 1 / (a_zuu + a_zud) # 0%--asset 1 has same price in all 3 nodes
Rf_1_zd <- 1 / (a_zdu + a_zdd) # 4%

# Expected and excess returns on asset 2
Er_2_0 <- mean(c(1.9444, 1.2742)) / 1.5285
Eer_2_0 <- Er_2_0 - Rf_0 + 1 # excess return, t = 0

Er_2_up <- mean(c(2.5, 1.5)) / 1.9444
Eer_2_up <- Er_2_up - Rf_1_zu + 1 # t = 1, up

Er_2_dn <- 1.5 / 1.2742
Eer_2_dn <- Er_2_dn - Rf_1_zd + 1 # t = 1, down

# Optimal consumption
c2_s <- (1 / a0) / 4
num_a0 <- c2_s * a0
