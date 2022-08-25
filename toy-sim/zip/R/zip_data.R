library(readr)
library(Matrix)
library(splines)
library(tidyverse)
library(spdep)
library(spatialreg)
library(extraDistr)


t <- 200 # timepoints
p <- 7 # parameters

# create correlation matrix from 3 levels of relationships using real ecoregions
load(file = "toy-sim/shared-data/region_key.RData")
level1_all8and9 <- seq(8,9.9,0.1)
# level1_all10and11 <- seq(10,11.9,0.1)

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

# %>%
#   filter(NA_L2CODE == 9.4 | NA_L2CODE == 8.3 | NA_L2CODE == 8.4)
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
bp_lambda <- runif(1)/2

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
cov_ar1_lambda <- equal + bp_lambda * bp_lin + bp_lambda^2 * bp_square + bp_lambda^3 * bp_cube + bp_lambda^4 * bp_quart
chol_ar1_lambda <- chol(cov_ar1_lambda)

# normal process --------
Z_lambda <- matrix(rnorm(r * p), r, p)

# create betas from std normal with AR(1) covariance matrix and 4 region correlation matrix
betas_lambda <- t(t(chol_corr) %*% Z_lambda %*% chol_ar1_lambda)

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
df_lambda <- matrix(NA, r, t)
for(i in 1:r) {
  df_lambda[i,] <- X_full[i, , ] %*% betas_lambda[, i]
}

# X_long <- X %>% as_tibble() %>% mutate(time = c(1:t)) %>%
#   rename_with(., ~ gsub("V", "reg", .x, fixed = TRUE)) %>%
#   pivot_longer(cols = c(1:all_of(r)), values_to = "linear", names_to = "region")

X_long <- X %>% as_tibble() %>%
  rename_with(., ~ reg_cols) %>% 
  mutate(time = c(1:t)) %>%
  pivot_longer(cols = c(1:all_of(r)), values_to = "linear", names_to = "region")

# lambda_effects <- t(df_lambda) %>% as_tibble() %>% mutate(time = c(1:t)) %>%
#   rename_with(., ~ gsub("V", "reg", .x, fixed = TRUE)) %>%
#   pivot_longer(cols = c(1:all_of(r)), values_to = "effect", names_to = "region") %>%
#   left_join(., X_long) %>%
#   left_join(., mod_reg_key) %>% mutate(type = "truth")

lambda_effects <- t(df_lambda) %>% as_tibble() %>% 
  rename_with(., ~ reg_cols) %>% 
  mutate(time = c(1:t)) %>%
  pivot_longer(cols = c(1:all_of(r)), values_to = "effect", names_to = "region") %>%
  left_join(., X_long) %>%
  left_join(., mod_reg_key) %>% mutate(type = "truth")

truth_lambda <- ggplot(lambda_effects, aes(x=linear, y=effect, group = region)) +
  geom_line(aes(linetype=NA_L1CODE, color = NA_L2CODE)) + labs(title = "Regression on lambda")
truth_lambda
# ggsave(paste0('manuscript/scripts/toy_sim/g1/rep_summit/truth_lambda_', 
#               format(as.POSIXlt(Sys.time(), "America/Denver"), "%d%b%Y_%H%M"),
#               ".pdf"), 
#        plot=truth_lambda, 
#        device = "pdf")

nb <- read_rds('~/Desktop/research/egpd-fires/data/processed/nb.rds')
ecoregions <- read_rds(file = "toy-sim/shared-data/ecoregions.RDS")
nb_agg <- aggregate(nb, ecoregions$NA_L3NAME)
nbInfo <- nb2WB(nb_agg)
nb_mat <- nb2mat(nb_agg, style = 'B')

W <- nb_mat # unnormalized weight matrix
D <- diag(nbInfo$num)

smallD <- D[idx, idx]
smallW <- W[idx, idx]
tau <- rexp(1)
eta <- runif(1)
alpha <- 0.99
phi_mat <- matrix(NA, t, r)
Q <- tau * (smallD - alpha * smallW)
phi_mat[1,] <- rnorm(rep(0,r), Q) # one timepoint
for(i in 2:t) {
  phi_mat[i,] <- rnorm(eta * phi_mat[i-1,], Q)
}

### generate neighborhood data for icar prior in stan
listw <- nb2listw(nb_agg, style = 'B', zero.policy = TRUE)
B <- as(listw, 'symmetricMatrix')
smallB <- B[idx, idx]
n_edges = length(smallB@i)
node1 = smallB@i+1 # add one to offset zero-based index
node2 = smallB@j+1

reg_lambda <- matrix(NA, r, t)
for(i in 1:r) {
  reg_lambda[i,] <- X_full[i, , ] %*% betas_lambda[, i]/5 + phi_mat[, i]/10
}
range(c(exp(reg_lambda)))
lambda_true <- c(exp(reg_lambda))
# y <- rep(NA, t*r)
# for(i in 1:(t*r)) y[i] <- rpois(1, lambda_true[i])
pi <- 0.1
y <- rep(NA, t*r)
for(i in 1:(t*r)) y[i] <- rzip(1, lambda_true[i], pi)
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
  pi_val = pi,
  
  N_edges = n_edges,
  node1 = node1,
  node2 = node2,
  
  #   true parameters to use in diagnostics post sampling
  truth = list(betas_lambda = betas_lambda, 
               rho1 = rho1, rho2 = rho2, 
               bp_lambda = bp_lambda, 
               phi_mat = phi_mat, tau = tau, eta = eta)
)
# toy_data_icarphi_st <- stan_d
write_rds(stan_d, 'manuscript/scripts/toy_sim/zip/data/zip_reg10and11only_pi-const.rds')
