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

t <- 252 # timepoints
p <- 7 # parameters

# create correlation matrix from 3 levels of relationships using real ecoregions
load(file = "./sim-study/shared-data/region_key.RData")

full_reg_key <- as_tibble(region_key) %>% 
  mutate(region = sprintf("reg%d", 1:84),
         NA_L2CODE = as.factor(NA_L2CODE),
         NA_L1CODE = as.factor(NA_L1CODE),
         NA_L3CODE = as.factor(NA_L3CODE))

# indices of selected regions, pulled from full region key

r <- dim(full_reg_key)[1]
reg_cols <- full_reg_key$region
level3 <- matrix(0, 84, 84)
level2 <- matrix(0, 84, 84)
level1 <- matrix(0, 84, 84)

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
l3 <- level3
l2 <- level2
l1 <- level1
rho1 <- 0.54; rho2 <- 0.45
corr <- l3 + rho2 * l2 + rho1 * l1
chol_corr <- chol(corr)

# add in AR prior to betas -------
bp_kappa <- runif(1)/2
bp_nu <- runif(1)/2
bp_xi <- runif(1)/2

# create indicator matrices to create AR(1) covariance matrix for use in data generation and in stan model
zeroes <- matrix(0, p, p)
equal <- diag(p)
bp_lin <- zeroes
bp_square <- zeroes
bp_cube <- zeroes
bp_quart <- zeroes
cov_vec_idx <- 3:7 # 1st element is global intercept; 2, 8, 14, 20, 26, and 32 are the linear terms of each of the 6 covariates
# cov_vec_idx <- c(3:7, 9:13, 15:19, 21:25, 27:31, 33:37) # 1st element is global intercept; 2, 8, 14, 20, 26, and 32 are the linear terms of each of the 6 covariates

for(i in cov_vec_idx) {
  for(j in cov_vec_idx) {
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
chol_ar1_kappa <- chol(cov_ar1_kappa)
cov_ar1_nu <- equal + bp_nu * bp_lin + bp_nu^2 * bp_square + bp_nu^3 * bp_cube + bp_nu^4 * bp_quart
chol_ar1_nu <- chol(cov_ar1_nu)
cov_ar1_xi <- equal + bp_xi * bp_lin + bp_xi^2 * bp_square + bp_xi^3 * bp_cube + bp_xi^4 * bp_quart
chol_ar1_xi <- chol(cov_ar1_xi)

# normal process --------
Z_kappa <- matrix(rnorm(r * p), p, r)
Z_nu <- matrix(rnorm(r * p), p, r)
Z_xi <- matrix(rnorm(r * p), p, r)

# create betas from std normal with AR(1) covariance matrix and 4 region correlation matrix
betas_kappa <- t(chol_ar1_kappa) %*% Z_kappa %*% chol_corr
betas_nu <- t(chol_ar1_nu) %*% Z_nu %*% chol_corr
betas_xi <- t(chol_ar1_xi) %*% Z_xi %*% chol_corr

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
df_nu <- matrix(NA, r, t)
df_xi <- matrix(NA, r, t)
for(i in 1:r) {
  df_kappa[i,] <- X_full[i, , ] %*% betas_kappa[, i]
  df_nu[i,] <- X_full[i, , ] %*% betas_nu[, i]
  df_xi[i,] <- X_full[i, , ] %*% betas_xi[, i]
}

X_long <- X %>% as_tibble() %>%
  rename_with(., ~ reg_cols) %>% 
  mutate(time = c(1:t)) %>%
  pivot_longer(cols = c(1:all_of(r)), values_to = "linear", names_to = "region")

kappa_effects <- t(df_kappa) %>% as_tibble() %>%
  rename_with(., ~ reg_cols) %>%
  mutate(time = c(1:t)) %>%
  pivot_longer(cols = c(1:all_of(r)), values_to = "effect", names_to = "region") %>%
  left_join(., X_long) %>%
  left_join(., full_reg_key) %>% mutate(type = "truth")

ggplot(kappa_effects, aes(x=linear, y=effect, group = region)) + 
  geom_line(aes(linetype=NA_L1CODE, color = NA_L2CODE))

nu_effects <- t(df_nu) %>% as_tibble() %>% 
  rename_with(., ~ reg_cols) %>% 
  mutate(time = c(1:t)) %>%
  pivot_longer(cols = c(1:all_of(r)), values_to = "effect", names_to = "region") %>%
  left_join(., X_long) %>%
  left_join(., full_reg_key) %>% mutate(type = "truth")

ggplot(nu_effects, aes(x=linear, y=effect, group = region)) + 
  geom_line(aes(linetype=NA_L1CODE, color = NA_L2CODE))

xi_effects <- t(df_xi) %>% as_tibble() %>% 
  rename_with(., ~ reg_cols) %>% 
  mutate(time = c(1:t)) %>%
  pivot_longer(cols = c(1:all_of(r)), values_to = "effect", names_to = "region") %>%
  left_join(., X_long) %>%
  left_join(., full_reg_key) %>% mutate(type = "truth")

ggplot(xi_effects, aes(x=linear, y=effect, group = region)) + 
  geom_line(aes(linetype=NA_L1CODE, color = NA_L2CODE))

nb <- read_rds('./sim-study/shared-data/nb.rds')
ecoregions <- read_rds(file = "./sim-study/shared-data/ecoregions.RDS")
nb_agg <- aggregate(nb, ecoregions$NA_L3NAME)
nbInfo <- nb2WB(nb_agg)
nb_mat <- nb2mat(nb_agg, style = 'B')

W <- nb_mat # unnormalized weight matrix
D <- diag(nbInfo$num)

smallD <- D
smallW <- W
tau <- rexp(1)
eta <- runif(1)
alpha <- 0.99
phi_mat_kappa <- matrix(NA, t, r)
phi_mat_nu <- matrix(NA, t, r)
phi_mat_xi <- matrix(NA, t, r)
Q <- tau * (smallD - alpha * smallW)
phi_mat_kappa[1,] <- rnorm(rep(0,r), Q) # first timepoint
phi_mat_nu[1,] <- rnorm(rep(0,r), Q) # first timepoint
phi_mat_xi[1,] <- rnorm(rep(0,r), Q) # first timepoint
for(i in 2:t) {
  phi_mat_kappa[i,] <- rnorm(eta * phi_mat_kappa[i-1,], Q)
  phi_mat_nu[i,] <- rnorm(eta * phi_mat_nu[i-1,], Q)
  phi_mat_xi[i,] <- rnorm(eta * phi_mat_xi[i-1,], Q)
}

### generate neighborhood data for icar prior in stan
listw <- nb2listw(nb_agg, style = 'B', zero.policy = TRUE)
B <- as(listw, 'symmetricMatrix')
smallB <- B
n_edges = length(smallB@i)
node1 = smallB@i+1 # add one to offset zero-based index
node2 = smallB@j+1

reg_kappa <- matrix(NA, t, r)
reg_nu <- matrix(NA, t, r)
reg_xi <- matrix(NA, t, r)

for(i in 1:r) {
  reg_kappa[,i] <- X_full[i, , ] %*% betas_kappa[, i]/4 + phi_mat_kappa[,i]/10
  reg_nu[,i] <- X_full[i, , ] %*% betas_nu[, i]/4 + phi_mat_nu[,i]/10
  reg_xi[,i] <- X_full[i, , ] %*% betas_xi[, i]/8 + phi_mat_xi[,i]/10
}
range(exp(reg_kappa))
range(exp(reg_nu))
range(exp(reg_xi))
kappa_true <- c(exp(reg_kappa))[idx_burns_all]
nu_true <- c(exp(reg_nu))[idx_burns_all]
xi_true <- c(exp(reg_xi))[idx_burns_all]

y <- rep(NA, length(idx_burns_all))
sigma_true <- rep(NA, length(idx_burns_all))
for(i in 1:length(idx_burns_all)) {
  sigma_true[i] <- nu_true[i]/(1 + xi_true[i])
  y[i] <- g1_random(n = 1, sigma = sigma_true[i], xi = xi_true[i], kappa = kappa_true[i])
  if (y[i] == 0) {
    y[i] = y[i] + 1e-10
  } else {
    y[i] = y[i]
  }
}

y_obs <- y[idx_burns_obs]
# range(sigma_true)

toy_data <- list(
  r = 84, # total number of regions
  t_tb = t,
  p = p,
  N_obs = length(idx_burns_obs),
  N_mis = length(idx_burns_mis),
  N_all = length(idx_burns_all),
  ii_obs = idx_burns_obs,
  ii_mis = idx_burns_mis,
  ii_all = idx_burns_all, # for broadcasting params in likelihood
  # M = 3, # of params based on data in egpd density
  
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
  
  n_edges = length(B@i),
  node1 = B@i + 1, # add one to offset zero-based index
  node2 = B@j + 1,
  
  # training data
  X = X_full,
  y_obs = y_obs
)
