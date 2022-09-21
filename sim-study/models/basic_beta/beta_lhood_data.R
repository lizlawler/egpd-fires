library(readr)
library(Matrix)
library(splines)
library(tidyverse)
library(spdep)
library(spatialreg)
library(assertthat)
library(LearnBayes)

n <- 1000 # timepoints
p <- 5 # parameters

# create correlation matrix from 3 levels of relationships using real ecoregions
load(file = "./sim-study/shared-data/region_key.RData")
level1_all8and9 <- seq(8,9.9,0.1)

mod_reg_key <- as_tibble(region_key) %>% 
  mutate(region = sprintf("reg%d", 1:84),
         NA_L2CODE = as.factor(NA_L2CODE),
         NA_L1CODE = as.factor(NA_L1CODE),
         NA_L3CODE = as.factor(NA_L3CODE)) %>%
  filter(NA_L2CODE %in% level1_all8and9)

full_reg_key <- as_tibble(region_key) %>% 
  mutate(region = sprintf("reg%d", 1:84),
         NA_L2CODE = as.factor(NA_L2CODE),
         NA_L1CODE = as.factor(NA_L1CODE),
         NA_L3CODE = as.factor(NA_L3CODE))

# indices of selected regions, pulled from full region key
idx <- which(full_reg_key$NA_L2CODE %in% mod_reg_key$NA_L2CODE)

r <- dim(mod_reg_key)[1]
reg_cols <- mod_reg_key$region
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
l3 <- level3[idx, idx]
l2 <- level2[idx, idx]
l1 <- level1[idx, idx]
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

# normal process --------
set.seed(123)
Z <- matrix(rnorm(r * p), r, p)
mu <- rep(0, r*p)
set.seed(123)
Z_mn <- matrix(rmnorm(1, mu, diag(r*p)), r, p)

# create betas from std normal with AR(1) covariance matrix and 4 region correlation matrix
assert_that(all(round(t(chol_corr) %*% chol_corr, 2) == round(corr, 2)))
assert_that(all(round(t(chol_ar1) %*% chol_ar1, 2) == round(cov_ar1, 2)))

betas <- t(chol_corr) %*% Z %*% chol_ar1
test <- t(chol_corr) %*% Z_mn %*% chol_ar1
assert_that(all(round(test, 2) == round(betas,2)))

full_cov <- cov_ar1 %x% corr
set.seed(123)
betas_mn <- matrix(rmnorm(1, mu, full_cov), r, p)
assert_that(all(round(betas_mn, 2) == round(betas,2)))

# get back Z
Z_rev <- solve(t(chol_corr)) %*% betas %*% solve(chol_ar1)
assert_that(all(round(Z_rev, 2) == round(Z,2)))

beta_data <- list(
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
  betas = betas
)
