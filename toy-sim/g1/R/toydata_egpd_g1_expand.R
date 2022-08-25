library(readr)
library(Matrix)
library(splines)
library(tidyverse)
library(spdep)
library(spatialreg)

g1_cdf_inv <- function(u, sigma = sigma, xi = xi, kappa = kappa) {
  (sigma/xi) * ((1-u^(1/kappa))^-xi - 1)
}
g1_random <- function(n = n, sigma = 1, xi = 0.5, kappa = 5) {
  u <- runif(n)
  return(g1_cdf_inv(u, sigma, xi, kappa))
}

# now adding in time component
# n <- 2000
t <- 1000
r <- 15 # regions
p <- 7 # parameters

# create correlation matrix from 3 levels of relationships using real ecoregions
load(file = "~/Desktop/csu/research/region_key.RData")
mod_reg_key <- as_tibble(region_key) %>% 
  mutate(region = sprintf("reg%d", 1:84),
         NA_L2CODE = as.factor(NA_L2CODE),
         NA_L1CODE = as.factor(NA_L1CODE)) %>% 
  select(3:5)
level3 <- matrix(0, 84, 84)
level2 <- matrix(0, 84, 84)
level1 <- matrix(0, 84, 84)
# rownames(corr_mat) <- region_key$NA_L3NAME
for(i in 1:84) {
  for(j in 1:84) {
    if (region_key[j, "NA_L3CODE"] == region_key[i, "NA_L3CODE"]) {
      level3[i, j] = 1 # i = j, diagonal of 1s
    } else if (region_key[j, "NA_L2CODE"] == region_key[i, "NA_L2CODE"]) {
      level1[i, j] = 1 # indicator for correlation at level 1
      level2[i, j] = 1 # indicator for correlation at level 2
    } else if (region_key[j, "NA_L1CODE"] == region_key[i, "NA_L1CODE"]) {
      level1[i, j] = 1 # indicator for correlation at level 1
    } else {
      level3[i, j] = 0
      level1[i, j] = 0
      level2[i, j] = 0
    }
  }
}
l3 <- level3[1:r, 1:r]
l2 <- level2[1:r, 1:r]
l1 <- level1[1:r, 1:r]
rho1 <- 0.54; rho2 <- 0.45
corr <- l3 + rho2 * l2 + rho1 * l1
chol_corr <- chol(corr)

# add in AR prior to betas -------
bp_kappa <- runif(1)/2
# bp_nu <- runif(1)/2
# bp_xi <- runif(1)/2
# create indicator matrices to create AR(1) covariance matrix for use in data generation and in stan model
zeroes <- matrix(0, p, p)
equal <- diag(p)
bp_lin <- zeroes
bp_square <- zeroes
bp_cube <- zeroes
bp_quart <- zeroes

for(i in 3:p) { # indices 1 and 2 are for the intercept column and the linear column
  for(j in 3:p) {
    if (i + 1 == j | i - 1 == j) {
      bp_lin[i, j] = 1
    } else if (i + 2 == j | i - 2 == j) {
      bp_square[i, j] = 1
    } else if (i + 3 == j | i - 3 == j) {
      bp_cube[i, j] = 1
    } else if (i + 4 == j | i - 4 == j) {
      bp_quart[i, j] = 1
    } else {
      bp_lin[i, j] = 0
      bp_square[i, j] = 0
      bp_cube[i, j] = 0
      bp_quart[i, j] = 0
    }
  }
}
cov_ar1_kappa <- equal + bp_kappa * bp_lin + bp_kappa^2 * bp_square + bp_kappa^3 * bp_cube + bp_kappa^4 * bp_quart
# cov_ar1_nu <- equal + bp_nu * bp_lin + bp_nu^2 * bp_square + bp_nu^3 * bp_cube + bp_nu^4 * bp_quart
# cov_ar1_xi <- equal + bp_xi * bp_lin + bp_xi^2 * bp_square + bp_xi^3 * bp_cube + bp_xi^4 * bp_quart
chol_ar1_kappa <- chol(cov_ar1_kappa)
# chol_ar1_nu <- chol(cov_ar1_nu)
# chol_ar1_xi <- chol(cov_ar1_xi)


# -----
# normal process --------
Z_kappa <- matrix(rnorm(r * p), r, p)
# Z_nu <- matrix(rnorm(r * p), r, p)
# Z_xi <- matrix(rnorm(r * p), r, p)

# create betas from std normal with AR(1) covariance matrix and 4 region correlation matrix
betas_kappa <- t(t(chol_corr) %*% Z_kappa %*% chol_ar1_kappa)
# betas_nu <- t(t(chol_corr) %*% Z_nu %*% chol_ar1_nu)
# betas_xi <- t(t(chol_corr) %*% Z_xi %*% chol_ar1_xi)

# ---------
# no time component ----
genSpline <- function(x, n, df = 5, degree) {
  basis <- bs(x = x, df = df, degree = degree, 
              Boundary.knots = range(x), intercept = FALSE)
  full_design <- cbind(rep(1, n), x, basis)
  return(full_design)
}

X <- matrix(rnorm(n*1), n, 1)
X_full <- genSpline(X, n, df = 5, degree = 3)
df_kappa <- X_full %*% betas_kappa

kappa_effects <- t(df_kappa) %>% as_tibble() %>% mutate(time = c(1:100))
  rename_with(., ~ gsub("V", "reg", .x, fixed = TRUE)) %>%
  mutate(design = as.vector(X)) %>%
  # change cols to appropriate number of regions; if 15 regions, make 1:15, etc
  pivot_longer(cols = c(1:15), values_to = "effect", names_to = "region") %>%
  left_join(., mod_reg_key) %>%
  mutate(type = "truth")

# addition of time component ------
genSpline <- function(x, t, r, df = 5, degree) {
  full_design <- array(NA, dim = c(r, t, df + 2))
  for(i in 1:r) {
    basis <- bs(x = x[, i], df = df, degree = degree, 
                Boundary.knots = range(x[, i]), intercept = FALSE)    
    full_design[i, ,] <- cbind(rep(1,t), x[, i], basis)
  }
  return(full_design)
}

X <- matrix(rnorm(t*r), t, r)
X_full <- genSpline(X, t, r, df = 5, degree = 3)
df_kappa <- matrix(NA, r, t)
for(i in 1:r) {
  df_kappa[i,] <- X_full[i, , ] %*% betas_kappa[, i]
}

X_long <- X %>% as_tibble() %>% mutate(time = c(1:t)) %>%
  rename_with(., ~ gsub("V", "reg", .x, fixed = TRUE)) %>%
  pivot_longer(cols = c(1:all_of(r)), values_to = "linear", names_to = "region")

kappa_effects <- t(df_kappa) %>% as_tibble() %>% mutate(time = c(1:t)) %>%
  rename_with(., ~ gsub("V", "reg", .x, fixed = TRUE)) %>%
  pivot_longer(cols = c(1:all_of(r)), values_to = "effect", names_to = "region") %>%
  left_join(., X_long) %>%
  left_join(., mod_reg_key)

truth_kappa <- ggplot(kappa_effects, aes(x=linear, y=effect, group = region)) +
  geom_line(aes(linetype=NA_L1CODE, color = NA_L2CODE)) + labs(title = "Regression on kappa")
truth_kappa

# df_nu <- X_full %*% betas_nu
# nu_effects <- df_nu %>% as_tibble() %>%
#   rename_with(., ~ gsub("V", "reg", .x, fixed = TRUE)) %>%
#   mutate(design = as.vector(X)) %>%
#   # change cols to appropriate number of regions; if 15 regions, make 1:15, etc
#   pivot_longer(cols = c(1:15), values_to = "effect", names_to = "region") %>%
#   left_join(., mod_reg_key) %>%
#   mutate(type = "truth")
# 
# truth_nu <- ggplot(nu_effects, aes(x=design, y=effect, group = region)) +
#   geom_line(aes(linetype=NA_L1CODE, color = NA_L2CODE)) + labs(title = "Regression on nu")
# truth_nu
# 
# df_xi <- X_full %*% betas_xi
# xi_effects <- df_xi %>% as_tibble() %>%
#   rename_with(., ~ gsub("V", "reg", .x, fixed = TRUE)) %>%
#   mutate(design = as.vector(X)) %>%
#   # change cols to appropriate ximber of regions; if 15 regions, make 1:15, etc
#   pivot_longer(cols = c(1:15), values_to = "effect", names_to = "region") %>%
#   left_join(., mod_reg_key) %>%
#   mutate(type = "truth")
# 
# truth_xi <- ggplot(xi_effects, aes(x=design, y=effect, group = region)) +
#   geom_line(aes(linetype=NA_L1CODE, color = NA_L2CODE)) + labs(title = "Regression on xi")
# truth_xi


nb <- read_rds('~/Desktop/csu/research/josephs_paper/data/processed/nb.rds')
ecoregions <- read_rds(file = "ecoregions.RDS")
nb_agg <- aggregate(nb, ecoregions$NA_L3NAME)
nbInfo <- nb2WB(nb_agg)
nb_mat <- nb2mat(nb_agg, style = 'B')

W <- nb_mat # unnormalized weight matrix
D <- diag(nbInfo$num)
smallW <- W[1:r, 1:r]
smallD <- D[1:r, 1:r]
tau <- rexp(1)
eta <- runif(1)
phi_mat <- matrix(NA, t, r)
phi_mat[1,] <- rnorm(rep(0,r), tau*(smallD-smallW)) # one timepoint
for(i in 2:t) {
  phi_mat[i,] <- rnorm(eta * phi_mat[i-1,], tau*(smallD-smallW))
}

### generate neighborhood data for car prior - this is what paper used ------
listw <- nb2listw(nb_agg, style = 'B', zero.policy = TRUE)
B <- as(listw, 'symmetricMatrix')
smallB <- B[1:r, 1:r]
n_edges = length(smallB@i)
node1 = smallB@i+1 # add one to offset zero-based index
node2 = smallB@j+1

XB_kappa <- matrix(NA, r, t)
XB_nu <- matrix(NA, r, t)
XB_xi <- matrix(NA, r, t)

betas_nu <- matrix(runif(p*r, -.25, 0.25), p, r)
betas_xi <- matrix(runif(p*r, -.25, .25), p, r)
for(i in 1:r) {
  XB_kappa[i,] <- X_full[i, , ] %*% betas_kappa[, i]/3 + phi_mat[,i]/6
  XB_nu[i, ] <- X_full[i, , ] %*% betas_nu[, i]
  XB_xi[i, ] <- X_full[i, , ] %*% betas_xi[, i]/4
}
range(exp(XB_kappa))
range(exp(XB_nu))
range(exp(XB_xi))
kappa_true <- c(exp(XB_kappa))
nu_true <- c(exp(XB_nu))
xi_true <- c(exp(XB_xi))

# kappa_true <- c(exp(X_full %*% betas_kappa/5 + phi_mat/10)) # creates a vector with first n elements for region 1, second n for region 2, etc
# range(kappa_true) # if too large, consider changing constant in denominator in previous line
# betas_nu <- matrix(runif(p*r, -.25, 0.25), p, r)
# nu_true <- c(exp(X_full %*% betas_nu))
# range(nu_true)
# betas_xi <- matrix(runif(p*r, -.25, .25), p, r)
# xi_true <- c(exp(X_full %*% betas_xi/5))
# range(xi_true) # if greater than 1, consider changing constant in denominator in previous line

y <- rep(NA, t*r)
sigma_true <- rep(NA, t*r)
for(i in 1:(t*r)) {
  sigma_true[i] <- nu_true[i]/(1 + xi_true[i])
  y[i] <- g1_random(n = 1, sigma = sigma_true[i], xi = xi_true[i], kappa = kappa_true[i])
  # if (y[i] == 0) {
  #   y[i] = y[i] + 1e-10
  # } else {
  #   y[i] = y[i]
  # }
}
range(sigma_true)
range(y)

stan_d <- list(
  t = t,
  p = p,
  r = r,
  
  # indicator matrices for region correlation
  l3 = l3,
  l2 = l2,
  l1 = l1,
  
  # indicator matrices for AR(1) process
  equal = equal,
  bp_lin = bp_lin,
  bp_square = bp_square,
  bp_cube = bp_cube,
  bp_quart = bp_quart,
  
  # data
  X = X_full,
  y = y,
  
  N_edges = n_edges,
  node1 = node1,
  node2 = node2,
  
  # true parameters to use in diagnostics post sampling
  truth = list(betas_kappa = betas_kappa, betas_nu = betas_nu, betas_xi = betas_xi, 
               rho1 = rho1, rho2 = rho2, 
               bp_kappa = bp_kappa,
               phi_mat = phi_mat,
               tau = tau,
               eta = eta)
)
toy_data_ar1_icarphi_st <- stan_d
write_rds(toy_data_ar1_icarphi_st, 'manuscript/scripts/toy_sim/g1/data/toy_data_ar1_icarphi_spatial_time.rds')

# ------ only one regression-----
# full_effects <- basis_df[[2]] %>% as_tibble() %>%
#   rename_with(., ~ gsub("V", "reg", .x, fixed = TRUE)) %>% 
#   mutate(design = as.vector(X)) %>% 
#   # change cols to appropriate number of regions; if 15 regions, make 1:15, etc
#   pivot_longer(cols = c(1:15), values_to = "effect", names_to = "region") %>%
#   left_join(., mod_reg_key)
# 
# truth_plot <- ggplot(full_effects, aes(x=design, y=effect, group = region)) + 
#   geom_line(aes(linetype=NA_L1CODE, color = NA_L2CODE))
# truth_plot

# -----
# phi offset with spatial AR process, static time -------
# source('R/import_data.R') # go to script and run first few lines to get nbInfo
# library(nimble)
# CM <- as.carCM(nbInfo$adj, nbInfo$weights, nbInfo$num)
# W <- nb_mat # unnormalized weight matrix
# C <- sweep(W, 1, rowSums(W), FUN="/") # row-normalized weight matrix (each row sums to 1)
# M <- diag(CM$M) # diagonal matrix of inverses of conditional variances
# mu_phi <- rep(0,r)
# cov_phi <- 0.5^2 * inverse(0.5 * (diag(r) - 0.5 * C[1:r, 1:r])) %*% M[1:r, 1:r]
# gamma_bounds <- carBounds(CM$C, nbInfo$adj, nbInfo$num, CM$M)
# phi_static <- rmnorm(1, mu_phi, varcov = cov_phi)


betas_kappa <- matrix(runif(p*r, -.25, 0.25), p, r)
betas_nu <- matrix(runif(p*r, -.25, 0.25), p, r)
betas_xi <- matrix(runif(p*r, -.25, .25), p, r)

X_final <- as.matrix(basis_df[[1]])
kappa_true <- c(exp(X_final %*% betas_kappa)) # creates a vector with first n elements for region 1, second n for region 2, etc
range(kappa_true) # if too large, consider changing constant in denominator in previous line
nu_true <- c(exp(X_final %*% betas_nu))
range(nu_true)
xi_true <- c(exp(X_final %*% betas_xi/15))
range(xi_true) # if greater than 1, consider changing constant in denominator in previous line

y <- rep(NA, n*r)
sigma_true <- rep(NA, n*r)
for(i in 1:(n*r)) {
  sigma_true[i] <- nu_true[i]/(1 + xi_true[i])
  y[i] <- g1_random(n = 1, sigma = sigma_true[i], xi = xi_true[i], kappa = kappa_true[i])
  if (y[i] == 0) {
    y[i] = y[i] + 1e-10
  } else {
    y[i] = y[i]
  }
}
range(sigma_true)
range(y)

stan_d <- list(
  n = n,
  p = p,
  r = r,
  
  # indicator matrices for region correlation
  l3 = l3,
  l2 = l2,
  l1 = l1,
  
  # indicator matrices for AR(1) process
  equal = equal,
  bp_lin = bp_lin,
  bp_square = bp_square,
  bp_cube = bp_cube,
  bp_quart = bp_quart,
  
  # data
  X = X_final,
  y = y,
  
  # true parameters to use in diagnostics post sampling
  # truth = list(betas_kappa = betas_kappa, betas_nu = betas_nu, betas_xi = betas_xi, rho1 = rho1, rho2 = rho2, bp = bp, alpha = alpha_true, phi = phi_true)
  truth = list(betas_kappa = betas_kappa, betas_nu = betas_nu, betas_xi = betas_xi, 
               rho1 = rho1, rho2 = rho2, 
               bp = bp)
)
toy_data_ar1_xi <- stan_d
write_rds(toy_data_ar1_xi, 'manuscript/scripts/toy_sim/g1/data/toy_data_ar1_xi.rds')

