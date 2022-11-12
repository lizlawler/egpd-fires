library(tidyverse)
library(lubridate)
library(rstan)
library(rgdal)
library(assertthat)
library(sf)
library(splines)
library(spdep)
library(spatialreg)

# Albers equal area (AEA) conic projection of North America
aea_proj <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
source("./full-model/process-data/merge-data.R")

# generate spatial neighbors
if (!file.exists('./full-model/data/processed/nb.rds')) {
  nb <- poly2nb(as(ecoregions, 'Spatial'))
  write_rds(nb, './full-model/data/processed/nb.rds')
} else {
  nb <- read_rds('./full-model/data/processed/nb.rds')
}

nb_agg <- aggregate(nb, ecoregions$NA_L3NAME)
nb_mat <- nb2mat(nb_agg, style = 'B')

# generate neighborhood data for car prior
listw <- nb2listw(nb_agg, style = 'B', zero.policy = TRUE)
listw$style
B <- as(listw, 'symmetricMatrix')
# B is suitable for building N, N_edges, node1, and node2
# following http://mc-stan.org/users/documentation/case-studies/icar_stan.html

# generate correlation indicato matrix matrices
# create correlation matrix from 3 levels of relationships using real ecoregions
load(file = "./sim-study/shared-data/region_key.RData")
reg_key <- as_tibble(region_key) %>% 
  mutate(region = c(1:84),
         NA_L2CODE = as.factor(NA_L2CODE),
         NA_L1CODE = as.factor(NA_L1CODE),
         NA_L3CODE = as.factor(NA_L3CODE))

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

test <- level3 + 0.45*level2 + 0.54*level1
choltest <- chol(test)
assert_that(all(test == round(t(choltest) %*% choltest, 2)))

# generate AR(1) indicator matrices for use in covariance matrix in stan
p <- 37
zeroes <- matrix(0, p, p)
equal <- diag(p)
bp_lin <- zeroes
bp_square <- zeroes
bp_cube <- zeroes
bp_quart <- zeroes
cov_vec_idx <- setdiff(c(1:37), c(1, seq(2,37, 6))) # pull indices related to splines only

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

## FINAL DATA PROCESSING ## -----------
ecoregion_df <- as(ecoregions, "Spatial") %>%
  data.frame

# get areas for each L3 ecoregion
area_df <- ecoregion_df %>%
  as_tibble() %>%
  group_by(NA_L3NAME) %>%
  summarize(area = sum(Shape_Area)) %>% arrange(NA_L3NAME)

burn_df <- mtbs %>% arrange(NA_L3NAME, ym)

count_df <- count_df %>%
  left_join(area_df) %>%
  arrange(NA_L3NAME, ym)

er_df <- dplyr::distinct(data.frame(ecoregions),
                       NA_L3NAME, NA_L2NAME, NA_L1NAME,
                       NA_L2CODE, NA_L1CODE) %>%
  as_tibble %>%
  filter(NA_L2NAME != 'UPPER GILA MOUNTAINS (?)')

er_covs <- ecoregion_summaries %>%
  left_join(er_df) %>%
  filter(!NA_L2NAME == "UPPER GILA MOUNTAINS (?)",
         year > 1983,
         ym <= max(mtbs$ym)) %>%
  mutate(log_housing_density = log(housing_density),
         pr = ifelse(pr < 0 , 0, pr)) %>%
  left_join(area_df) %>%
  droplevels() %>%
  mutate(er_ym = paste(NA_L3NAME, ym, sep = "_")) %>%
  arrange(NA_L3NAME, ym)

assert_that(length(setdiff(er_covs$NA_L3NAME, count_df$NA_L3NAME)) == 0)
assert_that(!anyDuplicated(er_covs))
er_covs$id <- 1:nrow(er_covs)


# Create training sets, including years from 1984 to cutoff_year + 1 (training: 1984-2011)
cutoff_year <- 2012

# count split
train_counts <- count_df %>%
  filter(year < cutoff_year) %>%
  left_join(er_covs)

holdout_counts <- count_df %>%
  filter(year >= cutoff_year) %>%
  left_join(er_covs)

idx_count_hold <- match(holdout_counts$er_ym, er_covs$er_ym)
assert_that(all(holdout_counts$er_ym == er_covs[idx_count_hold, ]$er_ym))

# burn split; including missing data
train_burns_full <- burn_df %>%
  filter(FIRE_YEAR < cutoff_year) %>%
  right_join(er_covs %>% filter(year < cutoff_year)) %>% arrange(NA_L3NAME, ym)

holdout_burns_full <-  burn_df %>%
  filter(FIRE_YEAR >= cutoff_year) %>%
  right_join(er_covs %>% filter(year >= cutoff_year)) %>% arrange(NA_L3NAME, ym)


# this data frame has no duplicate ecoregion X timestep combos
N <- length(unique(er_covs$NA_L3NAME))
T <- length(unique(er_covs$ym))

assert_that(identical(nrow(er_covs), N * T))

# Create b-splines for climate vars
# normalization attempt ------
# library(cellWise)
# vars <- c('log_housing_density', 'vs',
#           'pr', 'prev_12mo_precip', 'tmmx',
#           'rmin')
# 
# df_each <- 5
# deg_each <- 3
# X_bs <- list()
# X_lin <- list()
# X_bs_df <- list()
# for (i in seq_along(vars)) {
#   varname <- paste("lin", vars[i], sep = "_")
#   # incorporate data normalization HERE and then create splines
#   X_lin[[i]] <- transfo(er_covs[[vars[i]]], type = "YJ")$Zt
#   X_bs[[i]] <- bs(x = X_lin[[i]], df = df_each, degree = deg_each, 
#                   Boundary.knots = range(X_lin[[i]]), intercept = FALSE)
#   
#   X_bs_df[[i]] <- X_bs[[i]] %>% as_tibble()
#   names(X_bs_df[[i]]) <- paste('bs', vars[[i]], 1:df_each, sep = '_')
#   X_bs_df[[i]] <- X_bs_df[[i]] %>%
#     mutate(!!varname := X_lin[[i]]) %>%
#     relocate(!!varname, before = where(is.character))
# }
# X_bs_df <- bind_cols(X_bs_df)
# assert_that(!any(is.na(X_bs_df)))
# 
# X_full <- X_bs_df %>% mutate(er_ym = er_covs$er_ym, NA_L3NAME = er_covs$NA_L3NAME, year = er_covs$year) %>% 
#   mutate(intercept = 1, .before = lin_log_housing_density)
# X_train <- X_full %>% filter(year < cutoff_year) %>% select(-year)

# standardization attempt -------
vars <- c('log_housing_density', 'vs',
          'pr', 'prev_12mo_precip', 'tmmx',
          'rmin')

df_each <- 5
deg_each <- 3
X_bs <- list()
X_lin <- list()
X_bs_df <- list()

for (i in seq_along(vars)) {
  varname <- paste("lin", vars[i], sep = "_")
  # incorporate data normalization HERE and then create splines
  X_lin[[i]] <- (er_covs[[vars[i]]] - mean(er_covs[[vars[i]]]))/sd(er_covs[[vars[i]]])
  X_bs[[i]] <- bs(x = X_lin[[i]], df = df_each, degree = deg_each,
                  Boundary.knots = range(X_lin[[i]]), intercept = FALSE)
  
  X_bs_df[[i]] <- X_bs[[i]] %>% as_tibble()
  names(X_bs_df[[i]]) <- paste('bs', vars[[i]], 1:df_each, sep = '_')
  X_bs_df[[i]] <- X_bs_df[[i]] %>%
    mutate(!!varname := X_lin[[i]]) %>%
    relocate(!!varname, before = where(is.character))
}

X_bs_df <- bind_cols(X_bs_df)
assert_that(!any(is.na(X_bs_df)))

X_full <- X_bs_df %>% mutate(er_ym = er_covs$er_ym, NA_L3NAME = er_covs$NA_L3NAME, year = er_covs$year) %>% 
  mutate(intercept = 1, .before = lin_log_housing_density)
X_train <- X_full %>% filter(year < cutoff_year) %>% select(-year)

# split X_full and X_train into list of 84 design matrices, then reshape to an array for stan model
X_list_full <- lapply(split(X_full, X_full$NA_L3NAME), function(x) select(x, -NA_L3NAME))
assert_that(all(bind_rows(X_list_full)$er_ym == X_full$er_ym))
X_list_train <- lapply(split(X_train, X_train$NA_L3NAME), function(x) select(x, -NA_L3NAME))
assert_that(all(bind_rows(X_list_train)$er_ym == X_train$er_ym))

## BURN COUNTS ## ----------
# pull matrix of training counts and make sure it matches training covariates
nfire_list_train <- lapply(split(train_counts, train_counts$NA_L3NAME), function(x) select(x, c(n_fire, er_ym)))
assert_that(all(mapply(function(x, y) identical(x$er_ym, y$er_ym), X_list_train, nfire_list_train)) == TRUE)
nfire_matrix_train <- matrix(unlist(lapply(nfire_list_train, function(x) select(x, n_fire))), nrow(nfire_list_train[[1]]), 84)
iden_vec <- c()
for(i in 1:84) {
  iden_vec[i] <- all(nfire_matrix_train[, i] == nfire_list_train[[i]]$n_fire)
}
assert_that(all(iden_vec) == TRUE)

# pull matrix of holdout counts
nfire_list_hold <- lapply(split(holdout_counts, holdout_counts$NA_L3NAME), function(x) select(x, c(n_fire, er_ym)))
idx_count_hold <- match(nfire_list_hold[[1]]$er_ym, X_list_full[[1]]$er_ym)
assert_that(all(mapply(function(x, y) identical(x[idx_count_hold,]$er_ym, y$er_ym), X_list_full, nfire_list_hold)) == TRUE)
nfire_matrix_hold <- matrix(unlist(lapply(nfire_list_hold, function(x) select(x, n_fire))), nrow(nfire_list_hold[[1]]), 84)
iden_vec <- c()
for(i in 1:84) {
  iden_vec[i] <- all(nfire_matrix_hold[, i] == nfire_list_hold[[i]]$n_fire)
}
assert_that(all(iden_vec) == TRUE)

## BURN AREA ## ----------
# pull indices from train_burns_full for use in stan model
idx_tb_mis <- which(is.na(train_burns_full$BurnBndAc)) # indices of missing y
idx_tb_obs <- which(!is.na(train_burns_full$BurnBndAc)) # indices of observed y; subset of all of the rows in the training dataframe
idx_tb_all <- match(train_burns_full$er_ym, X_train$er_ym) # indices to broadcast kappa, sigma, and xi in the model for observed and missing y
assert_that(all(X_train$er_ym[idx_tb_all] == train_burns_full$er_ym)) # check broadcasting works
assert_that(all(X_train$er_ym[idx_tb_all][idx_tb_obs] == train_burns_full$er_ym[idx_tb_obs]))

# pull indices from holdout_burns_full for use in log scores
idx_hold_obs <- which(!is.na(holdout_burns_full$BurnBndAc)) # pulls indices wrt holdout dataset
idx_hold_all <- match(holdout_burns_full$er_ym, X_full$er_ym) # indices of holdout dataset wrt full X
assert_that(all(X_full$er_ym[idx_hold_all] == holdout_burns_full$er_ym)) 
assert_that(all(X_full$er_ym[idx_hold_all][idx_hold_obs] == holdout_burns_full$er_ym[idx_hold_obs]))

# original burn area
burn_train_obs_og <- train_burns_full$BurnBndAc[idx_tb_obs] - 1000
assert_that(all(!is.na(burn_train_obs_og)))
hist(burn_train_obs_og)
burn_hold_obs_og <- holdout_burns_full$BurnBndAc[idx_hold_obs] - 1000
assert_that(all(!is.na(burn_hold_obs_og)))
hist(burn_hold_obs_og)

# use square root of burn area 
burn_train_obs_sqrt <- sqrt(burn_train_obs_og)
assert_that(all(!is.na(burn_train_obs_sqrt)))
hist(burn_train_obs_sqrt)
burn_hold_obs_sqrt <- sqrt(burn_hold_obs_og)
assert_that(all(!is.na(burn_hold_obs_sqrt)))
hist(burn_hold_obs_sqrt)

# reshape design matrices into arrays for stan model
X_list_full <- lapply(X_list_full, function(x) select(x, -c(year, er_ym)))
X_list_train <- lapply(X_list_train, function(x) select(x, -er_ym))
t_train <- nrow(X_list_train[[1]])
t_all <- nrow(X_list_full[[1]])
X_array_train <- array(NA, dim = c(84, t_train, 37))
X_array_full <- array(NA, dim = c(84, t_all, 37))
for(i in 1:84) {
  X_array_train[i, ,] <- as.matrix(X_list_train[[i]])
  X_array_full[i, ,] <- as.matrix(X_list_full[[i]])
}

# Bundle up data into a list too pass to Stan -----------------------------
stan_data_og <- list(
  r = 84, # total number of regions
  p = p,
  t_all = t_all,
  t_train = t_train,
  t_hold = t_all - t_train,
  
  # covariate data
  X_full = X_array_full,
  X_train = X_array_train,
  
  # count data
  area_offset = log(area_df$area * 1e-11) / 2,
  y_train_count = nfire_matrix_train,
  y_hold_count = nfire_matrix_hold,
  idx_train_er = 1:t_train,
  idx_hold_er = idx_count_hold,

  # burn data 
  # training data
  y_train_obs = burn_train_obs_og,
  ii_tb_obs = idx_tb_obs,
  ii_tb_mis = idx_tb_mis,
  ii_tb_all = idx_tb_all, # for broadcasting params in likelihood
  N_tb_obs = length(idx_tb_obs),
  N_tb_mis = length(idx_tb_mis),
  N_tb_all = length(idx_tb_all),
  
  # holdout data
  y_hold_obs = burn_hold_obs_og,
  N_hold_obs = length(idx_hold_obs),
  N_hold_all = length(idx_hold_all),
  ii_hold_obs = idx_hold_obs,
  ii_hold_all = idx_hold_all, 
  
  # indicator matrices for region correlation
  l3 = level3,
  l2 = level2,
  l1 = level1,
  
  # indicator matrices for AR(1) process
  equal = equal,
  bp_lin = bp_lin,
  bp_square = bp_square,
  bp_cube = bp_cube,
  bp_quart = bp_quart,
  
  n_edges = length(B@i),
  node1 = B@i + 1, # add one to offset zero-based index
  node2 = B@j + 1,
  
  effects = list(reg_key = reg_key, vars = vars, X_lin = X_lin)
)

stan_data_sqrt <- list(
  r = 84, # total number of regions
  p = p,
  t_all = t_all,
  t_train = t_train,
  t_hold = t_all - t_train,
  
  # covariate data
  X_full = X_array_full,
  X_train = X_array_train,
  
  # count data
  area_offset = log(area_df$area * 1e-11) / 2,
  y_train_count = nfire_matrix_train,
  y_hold_count = nfire_matrix_hold,
  idx_train_er = 1:t_train,
  idx_hold_er = idx_count_hold,
  
  # burn data 
  # training data
  y_train_obs = burn_train_obs_sqrt,
  ii_tb_obs = idx_tb_obs,
  ii_tb_mis = idx_tb_mis,
  ii_tb_all = idx_tb_all, # for broadcasting params in likelihood
  N_tb_obs = length(idx_tb_obs),
  N_tb_mis = length(idx_tb_mis),
  N_tb_all = length(idx_tb_all),
  
  # holdout data
  y_hold_obs = burn_hold_obs_sqrt,
  N_hold_obs = length(idx_hold_obs),
  N_hold_all = length(idx_hold_all),
  ii_hold_obs = idx_hold_obs,
  ii_hold_all = idx_hold_all, 
  
  # indicator matrices for region correlation
  l3 = level3,
  l2 = level2,
  l1 = level1,
  
  # indicator matrices for AR(1) process
  equal = equal,
  bp_lin = bp_lin,
  bp_square = bp_square,
  bp_cube = bp_cube,
  bp_quart = bp_quart,
  
  n_edges = length(B@i),
  node1 = B@i + 1, # add one to offset zero-based index
  node2 = B@j + 1,
  
  effects = list(reg_key = reg_key, vars = vars, X_lin = X_lin)
)

# assert that there are no missing values in stan_d
assert_that(!any(lapply(stan_data_og, function(x) any(is.na(x))) %>% unlist))
assert_that(!any(lapply(stan_data_sqrt, function(x) any(is.na(x))) %>% unlist))

saveRDS(stan_data_og, file = './full-model/data/stan_data_og.RDS')
saveRDS(stan_data_sqrt, file = './full-model/data/stan_data_sqrt.RDS')
print('stan_data.rds written!')

