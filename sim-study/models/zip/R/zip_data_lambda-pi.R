library(readr)
library(Matrix)
library(splines)
library(tidyverse)
library(spdep)
library(spatialreg)
library(extraDistr)


t <- 300 # timepoints
p <- 37 # parameters

# create correlation matrix from 3 levels of relationships using real ecoregions
load(file = "./sim-study/shared-data/region_key.RData")
# level1_all8and9 <- seq(8,9.9,0.1)

# mod_reg_key <- as_tibble(region_key) %>% 
#   mutate(region = sprintf("reg%d", 1:84),
#          NA_L2CODE = as.factor(NA_L2CODE),
#          NA_L1CODE = as.factor(NA_L1CODE),
#          NA_L3CODE = as.factor(NA_L3CODE)) %>%
#   filter(NA_L2CODE %in% level1_all8and9)

full_reg_key <- as_tibble(region_key) %>% 
  mutate(region = sprintf("reg%d", 1:84),
         NA_L2CODE = as.factor(NA_L2CODE),
         NA_L1CODE = as.factor(NA_L1CODE),
         NA_L3CODE = as.factor(NA_L3CODE))

# indices of selected regions, pulled from full region key
idx <- which(full_reg_key$NA_L2CODE %in% mod_reg_key$NA_L2CODE)

# %>%
#   filter(NA_L2CODE == 9.4 | NA_L2CODE == 8.3 | NA_L2CODE == 8.4)
# r <- dim(mod_reg_key)[1]
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
# l3 <- level3[idx, idx]
# l2 <- level2[idx, idx]
# l1 <- level1[idx, idx]

l3 <- level3
l2 <- level2
l1 <- level1
rho1 <- 0.54; rho2 <- 0.45
corr <- l3 + rho2 * l2 + rho1 * l1
chol_corr <- chol(corr)

# add in AR prior to betas -------
bp_lambda <- runif(1)/2
bp_pi <- runif(1)/2

# create indicator matrices to create AR(1) covariance matrix for use in data generation and in stan model
# zeroes <- matrix(0, p, p)
# equal <- diag(p)
# bp_lin <- zeroes
# bp_square <- zeroes
# bp_cube <- zeroes
# bp_quart <- zeroes
# for(i in 3:p) { # indices 1 and 2 are for the intercept column and the linear column
#   for(j in 3:p) {
#     if (i + 1 == j | i - 1 == j) {
#       bp_lin[i, j] = 1
#     } else if (i + 2 == j | i - 2 == j) {
#       bp_square[i, j] = 1
#     } else if (i + 3 == j | i - 3 == j) {
#       bp_cube[i, j] = 1
#     } else if (i + 4 == j | i - 4 == j) {
#       bp_quart[i, j] = 1
#     } else {
#       bp_lin[i, j] = 0
#       bp_square[i, j] = 0
#       bp_cube[i, j] = 0
#       bp_quart[i, j] = 0
#     }
#   }
# }

zeroes <- matrix(0, p, p)
equal <- diag(p)
bp_lin <- zeroes
bp_square <- zeroes
bp_cube <- zeroes
bp_quart <- zeroes
cov_vec_idx <- c(3:7, 9:13, 15:19, 21:25, 27:31, 33:37) 

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

cov_ar1_lambda <- equal + bp_lambda * bp_lin + bp_lambda^2 * bp_square + bp_lambda^3 * bp_cube + bp_lambda^4 * bp_quart
chol_ar1_lambda <- chol(cov_ar1_lambda)
cov_ar1_pi <- equal + bp_pi * bp_lin + bp_pi^2 * bp_square + bp_pi^3 * bp_cube + bp_pi^4 * bp_quart
chol_ar1_pi <- chol(cov_ar1_pi)
# normal process --------
Z_lambda <- matrix(rnorm(r * p), p, r)
Z_pi <- matrix(rnorm(r * p), p, r)

# create betas from std normal with AR(1) covariance matrix and 4 region correlation matrix
betas_lambda <- t(chol_ar1_lambda) %*% Z_lambda %*% chol_corr
betas_pi <- t(chol_ar1_pi) %*% Z_pi %*% chol_corr

# addition of time component ------
genSpline <- function(x, t, r, lin_p, df = 5, degree) {
  full_design <- array(NA, dim = c(r, t, (lin_p*df+lin_p+1))) # lin_p*df = number of total basis splines; lin_p = linear param; +1 is for global intercept
  for(i in 1:r) {
    full_design[i,,1] <- rep(1,t)
    idx <- 2:7
    for(j in 1:lin_p) {
      basis <- bs(x = x[i, ,j], df = df, degree = degree, 
                  Boundary.knots = range(x[i, ,j]), intercept = FALSE)    
      full_design[i, ,idx] <- cbind(x[i, ,j], basis)
      idx <- idx + 6
    }
  }
  return(full_design)
} 
lin_p <- 6
X <- array(matrix(rnorm(t*lin_p), t, lin_p), dim = c(r, t, lin_p))
X_full <- genSpline(X, t, r, lin_p = lin_p, df = 5, degree = 3)
df_lambda_v1 <- matrix(NA, r, t)
df_lambda_v2 <- matrix(NA, r, t)
df_lambda_v3 <- matrix(NA, r, t)
df_lambda_v4 <- matrix(NA, r, t)
df_lambda_v5 <- matrix(NA, r, t)
df_lambda_v6 <- matrix(NA, r, t)
for(i in 1:r) {
  df_lambda_v1[i,] <- X_full[i, ,1:7 ] %*% betas_lambda[1:7, i]
  df_lambda_v2[i,] <- X_full[i, ,c(1,8:13) ] %*% betas_lambda[c(1,8:13), i]
  df_lambda_v3[i,] <- X_full[i, ,c(1,14:19) ] %*% betas_lambda[c(1,14:19), i]
  df_lambda_v4[i,] <- X_full[i, ,c(1,20:25) ] %*% betas_lambda[c(1,20:25), i]
  df_lambda_v5[i,] <- X_full[i, ,c(1,26:31) ] %*% betas_lambda[c(1,26:31), i]
  df_lambda_v6[i,] <- X_full[i, ,c(1,32:37) ] %*% betas_lambda[c(1,32:37), i]
}

X_new <- list()
for(i in 1:r) {
  X_new[[i]] <- X[i,,] %>% as_tibble() %>% 
    mutate(region = reg_cols[i], 
           time = c(1:t))
}

X_long_v1 <- X_new %>% bind_rows() %>% select(c(V1, region, time)) %>% rename(linear = V1)
X_long_v2 <- X_new %>% bind_rows() %>% select(c(V2, region, time)) %>% rename(linear = V2)
X_long_v3 <- X_new %>% bind_rows() %>% select(c(V3, region, time)) %>% rename(linear = V3)
X_long_v4 <- X_new %>% bind_rows() %>% select(c(V4, region, time)) %>% rename(linear = V4)
X_long_v5 <- X_new %>% bind_rows() %>% select(c(V5, region, time)) %>% rename(linear = V5)
X_long_v6 <- X_new %>% bind_rows() %>% select(c(V6, region, time)) %>% rename(linear = V6)

# %>%
#   rename_with(., ~ reg_cols) %>% 
#   mutate(time = c(1:t)) %>%
#   pivot_longer(cols = c(1:all_of(r)), values_to = "linear", names_to = "region")

lambda_effects_v1 <- t(df_lambda_v1) %>% as_tibble() %>% 
  rename_with(., ~ reg_cols) %>% 
  mutate(time = c(1:t)) %>%
  pivot_longer(cols = c(1:all_of(r)), values_to = "effect", names_to = "region") %>%
  left_join(., X_long_v1) %>%
  left_join(., full_reg_key) %>% mutate(type = "truth")

lambda_effects_v2 <- t(df_lambda_v2) %>% as_tibble() %>% 
  rename_with(., ~ reg_cols) %>% 
  mutate(time = c(1:t)) %>%
  pivot_longer(cols = c(1:all_of(r)), values_to = "effect", names_to = "region") %>%
  left_join(., X_long_v2) %>%
  left_join(., mod_reg_key) %>% mutate(type = "truth")

lambda_effects_v3 <- t(df_lambda_v3) %>% as_tibble() %>% 
  rename_with(., ~ reg_cols) %>% 
  mutate(time = c(1:t)) %>%
  pivot_longer(cols = c(1:all_of(r)), values_to = "effect", names_to = "region") %>%
  left_join(., X_long_v3) %>%
  left_join(., mod_reg_key) %>% mutate(type = "truth")

# ggplot(lambda_effects_v1, aes(x=linear, y=effect, group = region)) +
#   geom_line(aes(linetype=NA_L1CODE, color = NA_L2CODE))
# 
# ggplot(lambda_effects_v2, aes(x=linear, y=effect, group = region)) +
#   geom_line(aes(linetype=NA_L1CODE, color = NA_L2CODE))
# 
# ggplot(lambda_effects_v3, aes(x=linear, y=effect, group = region)) +
#   geom_line(aes(linetype=NA_L1CODE, color = NA_L2CODE))

# df_pi <- matrix(NA, r, t)
# for(i in 1:r) {
#   df_pi[i,] <- X_full[i, , ] %*% betas_pi[, i]
# }
# 
# pi_effects <- t(df_pi) %>% as_tibble() %>% 
#   rename_with(., ~ reg_cols) %>% 
#   mutate(time = c(1:t)) %>%
#   pivot_longer(cols = c(1:all_of(r)), values_to = "effect", names_to = "region") %>%
#   left_join(., X_long) %>%
#   left_join(., mod_reg_key) %>% mutate(type = "truth")


# truth_lambda <- ggplot(lambda_effects, aes(x=linear, y=effect, group = region)) +
#   geom_line(aes(linetype=NA_L1CODE, color = NA_L2CODE)) + labs(title = "Regression on lambda")
# truth_lambda

nb <- read_rds('./sim-study/shared-data/nb.rds')
ecoregions <- read_rds(file = "./sim-study/shared-data/ecoregions.RDS")
nb_agg <- aggregate(nb, ecoregions$NA_L3NAME)
nbInfo <- nb2WB(nb_agg)
nb_mat <- nb2mat(nb_agg, style = 'B')

W <- nb_mat # unnormalized weight matrix
D <- diag(nbInfo$num)

# smallD <- D[idx, idx]
# smallW <- W[idx, idx]
smallD <- D
smallW <- W
tau <- rexp(1)
eta <- runif(1)
alpha <- 0.99
phi_mat_lambda <- matrix(NA, t, r)
phi_mat_pi <- matrix(NA, t, r)
Q <- tau * (smallD - alpha * smallW)
phi_mat_lambda[1,] <- rnorm(rep(0,r), Q) # one timepoint
phi_mat_pi[1,] <- rnorm(rep(0,r), Q) # one timepoint
for(i in 2:t) {
  phi_mat_lambda[i,] <- rnorm(eta * phi_mat_lambda[i-1,], Q)
  phi_mat_pi[i,] <- rnorm(eta * phi_mat_pi[i-1,], Q)
}

### generate neighborhood data for icar prior in stan
listw <- nb2listw(nb_agg, style = 'B', zero.policy = TRUE)
B <- as(listw, 'symmetricMatrix')
# smallB <- B[idx, idx]
smallB <- B
n_edges = length(smallB@i)
node1 = smallB@i+1 # add one to offset zero-based index
node2 = smallB@j+1

reg_lambda <- matrix(NA, t, r)
reg_pi <- matrix(NA, t, r)
for(i in 1:r) {
  reg_lambda[,i] <- X_full[i, , ] %*% betas_lambda[, i]/10 + phi_mat_lambda[, i]/15
  reg_pi[,i] <- X_full[i, , ] %*% betas_pi[, i]/5 + phi_mat_pi[, i]/10
}
expit <- function(x) exp(x)/(1 + exp(x))
range(exp(reg_lambda))
range(expit(reg_pi))
lambda_true <- exp(reg_lambda)
pi_true <- expit(reg_pi)
# use empirical values by region for zero prob; generated in "make-stan-data" script, line 77
class(reg_zeroes)

y <- matrix(rep(NA, t*r), nrow = t, ncol = r)
for(i in 1:r) {
  for(j in 1:t) {
    y[j, i] <- rzip(1, lambda_true[j, i], reg_zeroes[i])
  }
}
# for(i in 1:(t*r)) y[i] <- rzip(1, lambda_true[i], pi_true[i])
range(y)

toy_data <- list(
  T = t,
  p = p,
  R = r,
  N = t*r,
  
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
  # pi_val = pi,
  
  n_edges = n_edges,
  node1 = node1,
  node2 = node2,
  
  # true parameters to use in diagnostics post sampling
  truth = list(betas_lambda = betas_lambda,
               betas_pi = betas_pi,
               rho1 = rho1, rho2 = rho2,
               bp_lambda = bp_lambda,
               bp_pi = bp_pi,
               phi_mat_lambda = phi_mat_lambda, 
               phi_mat_pi = phi_mat_pi, 
               tau = tau, eta = eta)
)
