library(readr)
library(LearnBayes)
library(Matrix)
library(splines)
library(tidyverse)


g3_cdf_inv <- function(u, sigma = sigma, xi = xi, delta = delta) {
  (sigma/xi) * ((qbeta((1-u), (1/delta), 2)^(-xi/delta)) - 1)
}

g3_rng <- function(n, sigma = sigma, xi = xi, delta = delta) {
  g3_cdf_inv(runif(n), sigma, xi, delta)
}

n <- 2000
r <- 10 # regions
p <- 7 # parameters

# create correlation matrix from 3 levels of relationships using real ecoregions
load(file = "region_key.RData")
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
bp <- runif(1)/2
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
cov_ar1 <- equal + bp * bp_lin + bp^2 * bp_square + bp^3 * bp_cube + bp^4 * bp_quart
chol_ar1 <- chol(cov_ar1)

# kron_cov <- cov_ar1 %x% corr
# kron_chol <- chol(kron_cov)

# normal process --------
Z_delta <- matrix(rnorm(r * p), r, p)
# Z_nu <- matrix(rmnorm(r * p, mean = rep(0, r*p), varcov = diag(r * p)), p, r)
# Z_xi <- matrix(rmnorm(r * p, mean = rep(0, r*p), varcov = diag(r * p)), p, r)

# create betas from std normal with AR(1) covariance matrix and 4 region correlation matrix
betas_delta <- t(t(chol_corr) %*% Z_delta %*% chol_ar1)

genSpline <- function(x, n, df = 5, degree, theta) {
  basis <- bs(x = x, df = df, degree = degree, 
              Boundary.knots = range(x), intercept = FALSE)
  full_design <- cbind(rep(1, n), x, basis)
  return(list(full_design, full_design %*% theta))
}
X <- matrix(rnorm(n*1), n, 1)
basis_df <- genSpline(X, n, df = 5, degree = 3, betas_delta)

full_effects <- basis_df[[2]] %>% as_tibble() %>%
  rename_with(., ~ gsub("V", "reg", .x, fixed = TRUE)) %>% 
  mutate(design = as.vector(X)) %>% 
  # change cols to appropriate number of regions; if 15 regions, make 1:15, etc
  pivot_longer(cols = c(1:10), values_to = "effect", names_to = "region") %>%
  left_join(., mod_reg_key)

truth_plot <- ggplot(full_effects, aes(x=design, y=effect, group = region)) + 
  geom_line(aes(linetype=NA_L1CODE, color = NA_L2CODE))
truth_plot

X_final <- as.matrix(basis_df[[1]])

# betas_delta <- matrix(runif(p*r, -.25, 0.25), p, r)
betas_nu <- matrix(runif(p*r, -.25, 0.25), p, r)
betas_xi <- matrix(runif(p*r, -.25, .25), p, r)

delta_true <- c(exp(X_final %*% betas_delta/3))
nu_true <- c(exp(X_final %*% betas_nu)) 
xi_true <- c(exp(X_final %*% betas_xi)) 

y <- rep(NA, n*r)
sigma_true <- rep(NA, n*r)
for(i in 1:(n*r)) {
  sigma_true[i] <- nu_true[i]/(1 + xi_true[i])
  y[i] <- g3_rng(1, sigma = sigma_true[i], xi = xi_true[i], delta = delta_true[i])
}

range(y)
range(delta_true)

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
  truth = list(betas_delta = betas_delta, betas_nu = betas_nu, betas_xi = betas_xi,
               rho1 = rho1, rho2 = rho2,
               bp = bp)
)

toy_data_g3_corr_ar1_10reg <- stan_d
write_rds(toy_data_g3_corr_ar1_10reg, 'manuscript/scripts/toy_sim/g3/data/toy_data_g3_corr_ar1_10reg.rds')

