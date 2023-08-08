library(tidyverse)
library(lubridate)
library(cmdstanr)
library(assertthat)
library(raster)
library(sp)
library(spdep)
library(spatialreg)
library(splines)

source("./full-model/data/process-data/merge-data.R")

# generate spatial neighbors
if (!file.exists('./full-model/data/processed/nb.rds')) {
  sf_use_s2(FALSE)
  nb <- poly2nb(as(ecoregion_shp, 'Spatial'))
  write_rds(nb, './full-model/data/processed/nb.rds')
} else {
  nb <- read_rds('./full-model/data/processed/nb.rds')
}

nb_agg <- aggregate(nb, ecoregion_shp$NA_L3NAME)
nb_mat <- nb2mat(nb_agg, style = 'B')

# generate neighborhood data for car prior
listw <- nb2listw(nb_agg, style = 'B', zero.policy = TRUE)
listw$style
B <- as(listw, 'symmetricMatrix')
# B is suitable for building N, N_edges, node1, and node2
# following http://mc-stan.org/users/documentation/case-studies/icar_stan.html

# generate correlation indicator matrix matrices
# create correlation matrix from 3 levels of relationships using real ecoregions
if (!file.exists('./full-model/data/processed/region_key.rds')) {
  region_key <- ecoregion_shp %>% 
    as_tibble() %>% 
    dplyr::select(NA_L3NAME, NA_L3CODE, NA_L2CODE, NA_L1CODE) %>% 
    unique() %>% arrange(NA_L3NAME, .locale = "en")
  write_rds(region_key, './full-model/data/processed/region_key.rds')
} else {
  region_key <- read_rds('./full-model/data/processed/region_key.rds')
}

assert_that(all(region_key$NA_L3NAME == attributes(nb_agg)$region.id))

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

## FINAL DATA PROCESSING ## -----------
ecoregion_df <- ecoregion_shp %>% as_tibble()

# get areas for each L3 ecoregion
area_df <- ecoregion_df %>%
  group_by(NA_L3NAME) %>%
  summarize(area = sum(Shape_Area)) %>% ungroup() %>%
  arrange(NA_L3NAME, .locale = "en")
assert_that(all(region_key$NA_L3NAME == area_df$NA_L3NAME))

burn_df <- mtbs_er %>% as_tibble() %>% 
  dplyr::select(BurnBndAc, fire_yr, fire_mon, ym, 
         NA_L3CODE, NA_L3NAME, NA_L2CODE, NA_L2NAME, NA_L1CODE, NA_L1NAME, 
         Shape_Leng, Shape_Area) %>%
  arrange(NA_L3NAME, ym, .locale = "en")

count_df <- count_df_climate %>% 
  left_join(count_df_erc) %>% left_join(count_df_fwi) %>%
  left_join(area_df) %>%
  arrange(NA_L3NAME, ym, .locale = "en")

er_df <- ecoregion_df %>% 
  dplyr::select(NA_L3NAME, NA_L2NAME, NA_L1NAME, NA_L2CODE, NA_L1CODE) %>% 
  distinct() %>%
  filter(NA_L2NAME != 'UPPER GILA MOUNTAINS (?)')

ecoregion_summaries <- ecoregion_summaries %>% 
  left_join(ecoregion_summaries_erc) %>%
  left_join(ecoregion_summaries_fwi)

er_covs <- ecoregion_summaries %>%
  left_join(er_df) %>%
  filter(NA_L2NAME != "UPPER GILA MOUNTAINS (?)",
         ym >= min(mtbs_er$ym),
         ym <= max(mtbs_er$ym)) %>%
  mutate(log_housing_density = log(housing_density),
         pr = ifelse(pr < 0 , 0, pr)) %>%
  left_join(area_df) %>%
  droplevels() %>%
  mutate(er_ym = paste(NA_L3NAME, ym, sep = "_")) %>%
  arrange(NA_L3NAME, ym, .locale = "en")

assert_that(length(setdiff(er_covs$NA_L3NAME, count_df$NA_L3NAME)) == 0)
assert_that(!anyDuplicated(er_covs))
assert_that(all(region_key$NA_L3NAME == unique(er_covs$NA_L3NAME)))
er_covs$id <- 1:nrow(er_covs)

# Create training and test sets
# test data: first five and last five years
# training data: middle 20 years
all_years <- 1990:2020
first_five <- 1990:1994
last_five <- 2020:2016
test_years <- sort(c(first_five, last_five))
train_years <- setdiff(all_years, test_years)

# count split
train_counts <- count_df %>%
  filter(year %in% train_years) %>%
  left_join(er_covs)
holdout_counts <- count_df %>%
  filter(year %in% test_years) %>%
  left_join(er_covs)
assert_that(sum(nrow(holdout_counts), nrow(train_counts)) == nrow(count_df))
assert_that(all(region_key$NA_L3NAME == unique(train_counts$NA_L3NAME)))
assert_that(all(region_key$NA_L3NAME == unique(holdout_counts$NA_L3NAME)))

idx_count_hold <- match(holdout_counts$er_ym, er_covs$er_ym)
assert_that(all(holdout_counts$er_ym == er_covs[idx_count_hold, ]$er_ym))

# burn split; including missing data
train_burns_full <- burn_df %>%
  filter(fire_yr %in% train_years) %>%
  right_join(er_covs %>% filter(year %in% train_years), 
             by = join_by("fire_yr" == "year", "fire_mon" == "month", NA_L3NAME)) %>% 
  arrange(NA_L3NAME, ym.y, .locale = "en")

holdout_burns_full <-  burn_df %>%
  filter(fire_yr %in% test_years) %>%
  right_join(er_covs %>% filter(year %in% test_years), 
             by = join_by("fire_yr" == "year", "fire_mon" == "month", NA_L3NAME)) %>% 
  arrange(NA_L3NAME, ym.y, .locale = "en")
assert_that(all(unique(holdout_burns_full$fire_yr) == test_years))
assert_that(all(unique(train_burns_full$fire_yr) == train_years))
assert_that(all(region_key$NA_L3NAME == unique(train_burns_full$NA_L3NAME)))
assert_that(all(region_key$NA_L3NAME == unique(holdout_burns_full$NA_L3NAME)))

# this data frame has no duplicate ecoregion X timestep combos
N <- as.numeric(length(unique(er_covs$NA_L3NAME)))
T <- as.numeric(length(unique(er_covs$ym)))
assert_that(nrow(er_covs) == N * T)

# standardization, then B-spline creation -------
vars <- c('log_housing_density', 'vs',
          'pr', 'prev_12mo_precip', 'tmmx',
          'rmin', 'erc', 'fwi')

df_each <- 5
deg_each <- 3
X_bs <- list()
X_lin <- list()
X_bs_df <- list()
un_std <- list()

for (i in seq_along(vars)) {
  varname <- paste("lin", vars[i], sep = "_")
  # data standardization
  re_scale_var <- t(c(mean(er_covs[[vars[i]]]), sd(er_covs[[vars[i]]])))
  X_lin[[i]] <- (er_covs[[vars[i]]] - mean(er_covs[[vars[i]]]))/sd(er_covs[[vars[i]]])
  
  # spline creation for each covariate
  X_bs[[i]] <- bs(x = X_lin[[i]], df = df_each, degree = deg_each,
                  Boundary.knots = range(X_lin[[i]]), intercept = FALSE)
  
  X_bs_df[[i]] <- X_bs[[i]] %>% as_tibble()
  names(X_bs_df[[i]]) <- paste('bs', vars[[i]], 1:df_each, sep = '_')
  X_bs_df[[i]] <- X_bs_df[[i]] %>%
    mutate(!!varname := X_lin[[i]]) %>%
    relocate(!!varname, before = where(is.character))
  un_std[[i]] <- re_scale_var %>% as_tibble() %>% rename(mean = V1, sd = V2) %>% mutate(variable = vars[[i]])
}

X_bs_df <- bind_cols(X_bs_df)
assert_that(!any(is.na(X_bs_df)))

un_std <- bind_rows(un_std)

X_full <- X_bs_df %>% mutate(er_ym = er_covs$er_ym, NA_L3NAME = er_covs$NA_L3NAME, year = er_covs$year) %>% 
  mutate(intercept = 1, .before = lin_log_housing_density)
X_train <- X_full %>% filter(year %in% train_years) %>% dplyr::select(-year)
erc_idx <- c(1, grep('housing', colnames(X_train)), grep('_erc', colnames(X_train)))
fwi_idx <- c(1, grep('housing', colnames(X_train)), grep('_fwi', colnames(X_train)))
erc_fwi_idx <- unique(c(erc_idx, fwi_idx))
climate_idx <- 1:37
assert_that(all(region_key$NA_L3NAME == unique(X_full$NA_L3NAME)))
assert_that(all(region_key$NA_L3NAME == unique(X_train$NA_L3NAME)))

# split X_full and X_train into list of 84 design matrices, then reshape to an array for stan model
X_list_full <- lapply(split(X_full, X_full$NA_L3NAME), function(x) dplyr::select(x, -NA_L3NAME))
assert_that(all(bind_rows(X_list_full)$er_ym == X_full$er_ym))
X_list_train <- lapply(split(X_train, X_train$NA_L3NAME), function(x) dplyr::select(x, -NA_L3NAME))
assert_that(all(bind_rows(X_list_train)$er_ym == X_train$er_ym))

## BURN COUNTS ## ----------
# pull matrix of training counts and make sure it matches training covariates
nfire_list_train <- lapply(split(train_counts, train_counts$NA_L3NAME), function(x) dplyr::select(x, c(n_fire, er_ym)))
assert_that(all(mapply(function(x, y) identical(x$er_ym, y$er_ym), X_list_train, nfire_list_train)) == TRUE)
nfire_matrix_train <- matrix(unlist(lapply(nfire_list_train, function(x) dplyr::select(x, n_fire))), nrow(nfire_list_train[[1]]), 84)
iden_vec <- c()
for(i in 1:84) {
  iden_vec[i] <- all(nfire_matrix_train[, i] == nfire_list_train[[i]]$n_fire)
}
assert_that(all(iden_vec) == TRUE)

# pull matrix of holdout counts
nfire_list_hold <- lapply(split(holdout_counts, holdout_counts$NA_L3NAME), function(x) dplyr::select(x, c(n_fire, er_ym)))
idx_count_hold <- match(nfire_list_hold[[1]]$er_ym, X_list_full[[1]]$er_ym)
assert_that(all(mapply(function(x, y) identical(x[idx_count_hold,]$er_ym, y$er_ym), X_list_full, nfire_list_hold)) == TRUE)
nfire_matrix_hold <- matrix(unlist(lapply(nfire_list_hold, function(x) dplyr::select(x, n_fire))), nrow(nfire_list_hold[[1]]), 84)
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
burn_train_obs_og <- train_burns_full$BurnBndAc[idx_tb_obs]/1000
assert_that(all(!is.na(burn_train_obs_og)))
hist(burn_train_obs_og)
burn_hold_obs_og <- holdout_burns_full$BurnBndAc[idx_hold_obs]/1000
assert_that(all(!is.na(burn_hold_obs_og)))
hist(burn_hold_obs_og)

# reshape design matrices into arrays for stan model
X_list_full <- lapply(X_list_full, function(x) dplyr::select(x, -c(year, er_ym)))
X_list_train <- lapply(X_list_train, function(x) dplyr::select(x, -er_ym))
t_train <- nrow(X_list_train[[1]])
t_all <- nrow(X_list_full[[1]])
X_array_train_climate <- array(NA, dim = c(84, t_train, 37))
X_array_full_climate <- array(NA, dim = c(84, t_all, 37))
X_array_train_erc <- array(NA, dim = c(84, t_train, 13))
X_array_full_erc <- array(NA, dim = c(84, t_all, 13))
X_array_train_fwi <- array(NA, dim = c(84, t_train, 13))
X_array_full_fwi <- array(NA, dim = c(84, t_all, 13))
X_array_train_erc_fwi <- array(NA, dim = c(84, t_train, 19))
X_array_full_erc_fwi <- array(NA, dim = c(84, t_all, 19))
for(i in 1:84) {
  X_array_train_climate[i, ,] <- as.matrix(X_list_train[[i]][,climate_idx])
  X_array_full_climate[i, ,] <- as.matrix(X_list_full[[i]][,climate_idx])
  X_array_train_erc[i, ,] <- as.matrix(X_list_train[[i]][,erc_idx])
  X_array_full_erc[i, ,] <- as.matrix(X_list_full[[i]][,erc_idx])
  X_array_train_fwi[i, ,] <- as.matrix(X_list_train[[i]][,fwi_idx])
  X_array_full_fwi[i, ,] <- as.matrix(X_list_full[[i]][,fwi_idx])
  X_array_train_erc_fwi[i, ,] <- as.matrix(X_list_train[[i]][,erc_fwi_idx])
  X_array_full_erc_fwi[i, ,] <- as.matrix(X_list_full[[i]][,erc_fwi_idx])
}

## MODEL WITH ALL CLIMATE COVARIATES ## -----------------------------------------
# generate AR(1) indicator matrices for use in covariance matrix
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

# points for twCRPS calculation
n_int <- 5000
int_holdout <- max(burn_hold_obs_og) - min(burn_hold_obs_og)
int_train <- max(burn_train_obs_og) - min(burn_train_obs_og)
int_pts_holdout <- min(burn_hold_obs_og) + (1:n_int)*(int_holdout/n_int)
int_pts_train <- min(burn_train_obs_og) + (1:n_int)*(int_train/n_int)

# Bundle up data into a list to pass to Stan -----------------------------
stan_data_climate <- list(
  R = 84, # total number of regions
  p = p,
  T_all = t_all,
  T_train = t_train,
  T_hold = t_all - t_train,
  
  # covariate data
  X_full = X_array_full_climate,
  X_train = X_array_train_climate,
  
  # count data
  area_offset = log(area_df$area * 1e-11) / 2,
  y_train_count = nfire_matrix_train,
  y_hold_count = nfire_matrix_hold,
  idx_train_er = 1:t_train,
  idx_hold_er = idx_count_hold,

  # burn data 
  # training data
  y_min = min(burn_hold_obs_og, burn_train_obs_og),
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
  
  # for twCRPS
  n_int = n_int,
  int_train = int_train,
  int_pts_train = int_pts_train,
  int_holdout = int_holdout,
  int_pts_holdout = int_pts_holdout
)

## JOINT MODEL ----------------------------
# climate covariates for count, ERC for burns
stan_data_joint <- list(
  R = 84, # total number of regions
  p = p,
  p_burn = 13,
  T_all = t_all,
  T_train = t_train,
  T_hold = t_all - t_train,
  
  # covariate data
  X_full_count = X_array_full_climate,
  X_train_count = X_array_train_climate,
  X_full_burn = X_array_full_erc,
  X_train_burn = X_array_train_erc,
  
  # count data
  area_offset = log(area_df$area * 1e-11) / 2,
  y_train_count = nfire_matrix_train,
  y_hold_count = nfire_matrix_hold,
  idx_train_er = 1:t_train,
  idx_hold_er = idx_count_hold,
  
  # burn data 
  # training data
  y_min = min(burn_hold_obs_og, burn_train_obs_og),
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
  
  # for twCRPS
  n_int = n_int,
  int_train = int_train,
  int_pts_train = int_pts_train,
  int_holdout = int_holdout,
  int_pts_holdout = int_pts_holdout
)

## MODEL WITH ERC and HOUSING COVARIATES ## -----------------------------------------
# use the first 13 x 13 block from the AR(1) covariance indicator matrices generated above
p <- 13

# Bundle up data into a list to pass to Stan -----------------------------
stan_data_erc <- list(
  R = 84, # total number of regions
  p = p,
  T_all = t_all,
  T_train = t_train,
  T_hold = t_all - t_train,
  
  # covariate data
  X_full = X_array_full_erc,
  X_train = X_array_train_erc,
  
  # count data
  area_offset = log(area_df$area * 1e-11) / 2,
  y_train_count = nfire_matrix_train,
  y_hold_count = nfire_matrix_hold,
  idx_train_er = 1:t_train,
  idx_hold_er = idx_count_hold,
  
  # burn data 
  # training data
  y_min = min(burn_hold_obs_og, burn_train_obs_og),
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
  equal = equal[1:p, 1:p],
  bp_lin = bp_lin[1:p, 1:p],
  bp_square = bp_square[1:p, 1:p],
  bp_cube = bp_cube[1:p, 1:p],
  bp_quart = bp_quart[1:p, 1:p],
  
  n_edges = length(B@i),
  node1 = B@i + 1, # add one to offset zero-based index
  node2 = B@j + 1,
  
  # for twCRPS
  n_int = n_int,
  int_train = int_train,
  int_pts_train = int_pts_train,
  int_holdout = int_holdout,
  int_pts_holdout = int_pts_holdout
)

## MODEL WITH FWI and HOUSING COVARIATES ## -----------------------------------------
# Bundle up data into a list to pass to Stan -----------------------------
stan_data_fwi <- list(
  R = 84, # total number of regions
  p = p,
  T_all = t_all,
  T_train = t_train,
  T_hold = t_all - t_train,
  
  # covariate data
  X_full = X_array_full_fwi,
  X_train = X_array_train_fwi,
  
  # count data
  area_offset = log(area_df$area * 1e-11) / 2,
  y_train_count = nfire_matrix_train,
  y_hold_count = nfire_matrix_hold,
  idx_train_er = 1:t_train,
  idx_hold_er = idx_count_hold,
  
  # burn data 
  # training data
  y_min = min(burn_hold_obs_og, burn_train_obs_og),
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
  equal = equal[1:p, 1:p],
  bp_lin = bp_lin[1:p, 1:p],
  bp_square = bp_square[1:p, 1:p],
  bp_cube = bp_cube[1:p, 1:p],
  bp_quart = bp_quart[1:p, 1:p],
  
  n_edges = length(B@i),
  node1 = B@i + 1, # add one to offset zero-based index
  node2 = B@j + 1,
  
  # for twCRPS
  n_int = n_int,
  int_train = int_train,
  int_pts_train = int_pts_train,
  int_holdout = int_holdout,
  int_pts_holdout = int_pts_holdout
)

## MODEL WITH ERC, FWI and HOUSING COVARIATES ## -----------------------------------------
p <- 19
# Bundle up data into a list to pass to Stan -----------------------------
stan_data_erc_fwi <- list(
  R = 84, # total number of regions
  p = p,
  T_all = t_all,
  T_train = t_train,
  T_hold = t_all - t_train,
  
  # covariate data
  X_full = X_array_full_erc_fwi,
  X_train = X_array_train_erc_fwi,
  
  # count data
  area_offset = log(area_df$area * 1e-11) / 2,
  y_train_count = nfire_matrix_train,
  y_hold_count = nfire_matrix_hold,
  idx_train_er = 1:t_train,
  idx_hold_er = idx_count_hold,
  
  # burn data 
  # training data
  y_min = min(burn_hold_obs_og, burn_train_obs_og),
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
  equal = equal[1:p, 1:p],
  bp_lin = bp_lin[1:p, 1:p],
  bp_square = bp_square[1:p, 1:p],
  bp_cube = bp_cube[1:p, 1:p],
  bp_quart = bp_quart[1:p, 1:p],
  
  n_edges = length(B@i),
  node1 = B@i + 1, # add one to offset zero-based index
  node2 = B@j + 1,
  
  # for twCRPS
  n_int = n_int,
  int_train = int_train,
  int_pts_train = int_pts_train,
  int_holdout = int_holdout,
  int_pts_holdout = int_pts_holdout
)

# assert that there are no missing values in stan_d
assert_that(!any(lapply(stan_data_climate, function(x) any(is.na(x))) %>% unlist))
assert_that(!any(lapply(stan_data_erc, function(x) any(is.na(x))) %>% unlist))
assert_that(!any(lapply(stan_data_fwi, function(x) any(is.na(x))) %>% unlist))
assert_that(!any(lapply(stan_data_erc_fwi, function(x) any(is.na(x))) %>% unlist))
assert_that(!any(lapply(stan_data_joint, function(x) any(is.na(x))) %>% unlist))
saveRDS(stan_data_climate, file = './full-model/data/stan_data_climate.RDS')
saveRDS(stan_data_erc, file = './full-model/data/stan_data_erc.RDS')
saveRDS(stan_data_fwi, file = './full-model/data/stan_data_fwi.RDS')
saveRDS(stan_data_erc_fwi, file = './full-model/data/stan_data_erc_fwi.RDS')
saveRDS(stan_data_joint, file = './full-model/data/stan_data_joint.RDS')
write_stan_json(data = stan_data_climate, file = './full-model/data/stan_data_climate.json')
write_stan_json(data = stan_data_erc, file = './full-model/data/stan_data_erc.json')
write_stan_json(data = stan_data_fwi, file = './full-model/data/stan_data_fwi.json')
write_stan_json(data = stan_data_erc_fwi, file = './full-model/data/stan_data_erc_fwi.json')
write_stan_json(data = stan_data_joint, file = './full-model/data/stan_data_joint.json')
saveRDS(un_std, file = './full-model/data/un_std_all.RDS')
