##############################################################################
##  This script does the following:                                         ##
##  (1) Standardizes the covariates and generates B spline bases for each   ##
##  (2) Splits the count and size datasets into training and test sets      ##
##  (3) Creates list of necessary data elements to be imported into         ##
##      Stan model.                                                         ##
##############################################################################

source("./data/03_merge_data.R")
library(tidyverse)
library(lubridate)
library(cmdstanr)
library(assertthat)
library(raster)
library(sp)
library(spdep)
library(spatialreg)
library(splines)

###################################################
##  SPATIAL NEIGHBORS AND HIERACHICAL NESTING   ##
##################################################

if (!file.exists('./data/processed/nb.rds')) {
  sf_use_s2(FALSE)
  nb <- poly2nb(as(ecoregion_shp, 'Spatial'))
  write_rds(nb, './data/processed/nb.rds')
} else {
  nb <- read_rds('./data/processed/nb.rds')
}

nb_agg <- aggregate(nb, ecoregion_shp$NA_L3NAME)
nb_mat <- nb2mat(nb_agg, style = 'B')

# generate neighborhood data for car prior
listw <- nb2listw(nb_agg, style = 'B', zero.policy = TRUE)
listw$style
B <- as(listw, 'symmetricMatrix')
# B is suitable for building N, N_edges, node1, and node2
# see: http://mc-stan.org/users/documentation/case-studies/icar_stan.html

# generate correlation indicator matrix matrices
# create correlation matrix from 3 levels of relationships using real ecoregions
if (!file.exists('./data/processed/region_key.rds')) {
  region_key <- ecoregion_shp %>% 
    as_tibble() %>% 
    dplyr::select(NA_L1NAME, NA_L3NAME, NA_L3CODE, NA_L2CODE, NA_L1CODE) %>% 
    unique() %>% arrange(NA_L3NAME, .locale = "en") %>%
    mutate(region = c(1:84),
           NA_L2CODE = as.factor(NA_L2CODE),
           NA_L1CODE = as.factor(NA_L1CODE),
           NA_L3CODE = as.factor(NA_L3CODE))
  write_rds(region_key, './data/processed/region_key.rds')
} else {
  region_key <- read_rds('./data/processed/region_key.rds')
}

assert_that(all(region_key$NA_L3NAME == attributes(nb_agg)$region.id))

# initialize empty matrices for each level of ecoregion
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

###############################
##  FINAL DATA PROCESSING    ##
###############################
ecoregion_df <- ecoregion_shp %>% as_tibble() 
# aggregate areas for each L3 ecoregion
area_df <- ecoregion_df %>%
  group_by(NA_L3NAME) %>%
  summarize(area = sum(Shape_Area)) %>% ungroup() %>%
  arrange(NA_L3NAME, .locale = "en")
assert_that(all(region_key$NA_L3NAME == area_df$NA_L3NAME))

# aggregate burn areas for each L3 ecoregion
size_df <- mtbs_er %>% as_tibble() %>% 
  dplyr::select(BurnBndAc, fire_yr, fire_mon, ym, 
         NA_L3CODE, NA_L3NAME, NA_L2CODE, NA_L2NAME, NA_L1CODE, NA_L1NAME, 
         Shape_Leng, Shape_Area) %>%
  arrange(NA_L3NAME, ym, .locale = "en")

size_df_agg <- size_df %>%
  mutate(across(where(is.character), as.factor)) %>%
  group_by(NA_L3CODE, fire_yr, NA_L3NAME, NA_L1NAME, NA_L1CODE, NA_L2NAME, NA_L2CODE) %>%
  summarise(total_burns = sum(BurnBndAc)) %>%
  ungroup()
saveRDS(size_df_agg, file = "data/processed/obs_burned_areas.RDS")

count_df <- count_df %>%
  arrange(NA_L3NAME, ym, .locale = "en")

er_df <- ecoregion_df %>% 
  dplyr::select(NA_L3NAME, NA_L2NAME, NA_L1NAME, NA_L2CODE, NA_L1CODE) %>% 
  distinct()

er_covs <- ecoregion_summaries %>%
  left_join(er_df) %>%
  filter(ym >= min(mtbs_er$ym),
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
# check that er_covs df has no duplicate ecoregion X timestep combos
N <- as.numeric(length(unique(er_covs$NA_L3NAME)))
T <- as.numeric(length(unique(er_covs$ym)))
assert_that(nrow(er_covs) == N * T)

##  Create training and test sets --------------------------------------------
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

# burn split; including 'missing' data
train_sizes_full <- size_df %>%
  filter(fire_yr %in% train_years) %>%
  right_join(er_covs %>% filter(year %in% train_years), 
             by = join_by("fire_yr" == "year", "fire_mon" == "month", NA_L3NAME)) %>% 
  arrange(NA_L3NAME, ym.y, .locale = "en")
holdout_sizes_full <-  size_df %>%
  filter(fire_yr %in% test_years) %>%
  right_join(er_covs %>% filter(year %in% test_years), 
             by = join_by("fire_yr" == "year", "fire_mon" == "month", NA_L3NAME)) %>% 
  arrange(NA_L3NAME, ym.y, .locale = "en")
assert_that(all(unique(holdout_sizes_full$fire_yr) == test_years))
assert_that(all(unique(train_sizes_full$fire_yr) == train_years))
assert_that(all(region_key$NA_L3NAME == unique(train_sizes_full$NA_L3NAME)))
assert_that(all(region_key$NA_L3NAME == unique(holdout_sizes_full$NA_L3NAME)))

# data standardization, then B-spline creation -------
vars <- c('log_housing_density', 'vs',
          'pr', 'prev_12mo_precip', 'tmmx',
          'rmin', 'erc', 'fwi')

df_each <- 5
deg_each <- 3
X_bs <- list()
X_lin <- list()
X_bs_df <- list()
# store means and sd of data to back transform for results later on
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
saveRDS(un_std, file = './data/processed/un_std_all.RDS')

X_full <- X_bs_df %>% mutate(er_ym = er_covs$er_ym, NA_L3NAME = er_covs$NA_L3NAME, year = er_covs$year) %>% 
  mutate(intercept = 1, .before = lin_log_housing_density)
X_train <- X_full %>% filter(year %in% train_years) %>% dplyr::select(-year)
erc_idx <- c(1, grep('housing', colnames(X_train)), grep('_erc', colnames(X_train)))
fwi_idx <- c(1, grep('housing', colnames(X_train)), grep('_fwi', colnames(X_train)))
erc_fwi_idx <- unique(c(erc_idx, fwi_idx))
climate_idx <- 1:37
assert_that(all(region_key$NA_L3NAME == unique(X_full$NA_L3NAME)))
assert_that(all(region_key$NA_L3NAME == unique(X_train$NA_L3NAME)))

# split X_full and X_train into list of 84 design matrices, then reshape into an array
X_list_full <- lapply(split(X_full, X_full$NA_L3NAME), function(x) dplyr::select(x, -NA_L3NAME))
assert_that(all(bind_rows(X_list_full)$er_ym == X_full$er_ym))
X_list_train <- lapply(split(X_train, X_train$NA_L3NAME), function(x) dplyr::select(x, -NA_L3NAME))
assert_that(all(bind_rows(X_list_train)$er_ym == X_train$er_ym))

## Wildfire counts ## ----------
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

## Wildfire burned areas ## ----------
# pull indices from train_sizes_full for use in stan model
idx_ts_mis <- which(is.na(train_sizes_full$BurnBndAc)) # indices of 'missing' burned areas
idx_ts_obs <- which(!is.na(train_sizes_full$BurnBndAc)) # indices of observed burned areas; subset of all of the rows in the training dataframe
idx_ts_all <- match(train_sizes_full$er_ym, X_train$er_ym) # indices to broadcast kappa, sigma, and xi in the model for observed and missing sizes
assert_that(all(X_train$er_ym[idx_ts_all] == train_sizes_full$er_ym)) # check broadcasting works
assert_that(all(X_train$er_ym[idx_ts_all][idx_ts_obs] == train_sizes_full$er_ym[idx_ts_obs]))

# pull indices from holdout_sizes_full for use in log scores
idx_hold_obs <- which(!is.na(holdout_sizes_full$BurnBndAc)) # pulls indices wrt holdout dataset
idx_hold_all <- match(holdout_sizes_full$er_ym, X_full$er_ym) # indices of holdout dataset wrt full X
assert_that(all(X_full$er_ym[idx_hold_all] == holdout_sizes_full$er_ym)) 
assert_that(all(X_full$er_ym[idx_hold_all][idx_hold_obs] == holdout_sizes_full$er_ym[idx_hold_obs]))

# scale burned areas down by 1000
size_train_obs <- train_sizes_full$BurnBndAc[idx_ts_obs]/1000
assert_that(all(!is.na(size_train_obs)))
size_hold_obs <- holdout_sizes_full$BurnBndAc[idx_hold_obs]/1000
assert_that(all(!is.na(size_hold_obs)))

# reshape design matrices into arrays for use in Stan
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

# generate points for twCRPS calculation
n_int <- 5000
int_range <- max(size_hold_obs) - min(size_hold_obs)
int_pts <- min(size_hold_obs) + (1:n_int)*(int_range/n_int)

###########################################
##  GENERATE DATA LISTS FOR USE IN STAN  ##
###########################################
##  Model using climate covariates and housing density ------------------
stan_data_climate <- list(
  R = 84,
  p = p,
  T_all = t_all,
  T_train = t_train,
  T_hold = t_all - t_train,
  
  # covariate data
  X_full = X_array_full_climate,
  X_train = X_array_train_climate,
  
  # count data
  area_offset = log(area_df$area * 1e-11) / 2,
  # the coefficient of an offset variable in Poisson regression is always 1, so we don't want the offset itself
  # to be huge since it can't be rescaled in the regression; but we still want the magnitude of the areas
  # to be reflected in the offset. Log of area/10^11 accounts for this, and we further minimize the scale of the offset by
  # dividing by a factor of 2
  y_train_count = nfire_matrix_train,
  y_hold_count = nfire_matrix_hold,
  idx_train_er = setdiff(1:372, idx_count_hold),
  idx_hold_er = idx_count_hold,

  # burn data 
  # training data
  y_min = min(size_hold_obs, size_train_obs),
  y_train_obs = size_train_obs,
  ii_ts_obs = idx_ts_obs,
  ii_ts_mis = idx_ts_mis,
  ii_ts_all = idx_ts_all, # for broadcasting params in likelihood
  N_ts_obs = length(idx_ts_obs),
  N_ts_mis = length(idx_ts_mis),
  N_ts_all = length(idx_ts_all),
  
  # holdout data
  y_hold_obs = size_hold_obs,
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
  
  # points for twCRPS
  n_int = n_int,
  int_range = int_range,
  int_pts = int_pts
)

##  Model using ERC and housing density ----------------------------------
p <- as.numeric(dim(X_array_full_erc)[3])
# use the first 13 x 13 block from the AR(1) covariance indicator matrices generated above
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
  idx_train_er = setdiff(1:372, idx_count_hold),
  idx_hold_er = idx_count_hold,
  
  # burn data 
  # training data
  y_min = min(size_hold_obs, size_train_obs),
  y_train_obs = size_train_obs,
  ii_ts_obs = idx_ts_obs,
  ii_ts_mis = idx_ts_mis,
  ii_ts_all = idx_ts_all, # for broadcasting params in likelihood
  N_ts_obs = length(idx_ts_obs),
  N_ts_mis = length(idx_ts_mis),
  N_ts_all = length(idx_ts_all),
  
  # holdout data
  y_hold_obs = size_hold_obs,
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
  
  # points for twCRPS
  n_int = n_int,
  int_range = int_range,
  int_pts = int_pts
)

##  Model using FWI and housing density -------------------------------------
p <- as.numeric(dim(X_array_full_fwi)[3])
stan_data_fwi <- list(
  R = 84, 
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
  idx_train_er = setdiff(1:372, idx_count_hold),
  idx_hold_er = idx_count_hold,
  
  # burn data 
  # training data
  y_min = min(size_hold_obs, size_train_obs),
  y_train_obs = size_train_obs,
  ii_ts_obs = idx_ts_obs,
  ii_ts_mis = idx_ts_mis,
  ii_ts_all = idx_ts_all, # for broadcasting params in likelihood
  N_ts_obs = length(idx_ts_obs),
  N_ts_mis = length(idx_ts_mis),
  N_ts_all = length(idx_ts_all),
  
  # holdout data
  y_hold_obs = size_hold_obs,
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
  
  # points for twCRPS
  n_int = n_int,
  int_range = int_range,
  int_pts = int_pts
)

##  Model using ERC, FWI, and housing density -------------------------------
p <- as.numeric(dim(X_array_full_erc_fwi)[3])
stan_data_erc_fwi <- list(
  R = 84,
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
  idx_train_er = setdiff(1:372, idx_count_hold),
  idx_hold_er = idx_count_hold,
  
  # burn data 
  # training data
  y_min = min(size_hold_obs, size_train_obs),
  y_train_obs = size_train_obs,
  ii_ts_obs = idx_ts_obs,
  ii_ts_mis = idx_ts_mis,
  ii_ts_all = idx_ts_all, # for broadcasting params in likelihood
  N_ts_obs = length(idx_ts_obs),
  N_ts_mis = length(idx_ts_mis),
  N_ts_all = length(idx_ts_all),
  
  # holdout data
  y_hold_obs = size_hold_obs,
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
  
  # points for twCRPS
  n_int = n_int,
  int_range = int_range,
  int_pts = int_pts
)

##  Joint model, which uses climate covariates and housing density for 
##  count sub-model and ERC, FWI, and housing density for sizes
p <- 37
stan_data_joint <- list(
  R = 84,
  p = p,
  p_size = as.numeric(dim(X_array_full_erc_fwi)[3]),
  T_all = t_all,
  T_train = t_train,
  T_hold = t_all - t_train,
  
  # covariate data
  X_full_count = X_array_full_climate,
  X_train_count = X_array_train_climate,
  X_full_size = X_array_full_erc_fwi,
  X_train_size = X_array_train_erc_fwi,
  
  # count data
  area_offset = log(area_df$area * 1e-11) / 2,
  y_train_count = nfire_matrix_train,
  y_hold_count = nfire_matrix_hold,
  idx_train_er = setdiff(1:372, idx_count_hold),
  idx_hold_er = idx_count_hold,
  
  # sizes data 
  # training data
  y_min = min(size_hold_obs, size_train_obs),
  y_train_size_obs = size_train_obs,
  ii_ts_obs = idx_ts_obs,
  ii_ts_mis = idx_ts_mis,
  ii_ts_all = idx_ts_all, # for broadcasting params in likelihood
  N_ts_obs = length(idx_ts_obs),
  N_ts_mis = length(idx_ts_mis),
  N_ts_all = length(idx_ts_all),
  
  # holdout data
  y_hold_size_obs = size_hold_obs,
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
  
  # points for twCRPS
  n_int = n_int,
  int_range = int_range,
  int_pts = int_pts
)

# check that there are no missing values in any of these lists
assert_that(!any(lapply(stan_data_climate, function(x) any(is.na(x))) %>% unlist))
assert_that(!any(lapply(stan_data_erc, function(x) any(is.na(x))) %>% unlist))
assert_that(!any(lapply(stan_data_fwi, function(x) any(is.na(x))) %>% unlist))
assert_that(!any(lapply(stan_data_erc_fwi, function(x) any(is.na(x))) %>% unlist))
assert_that(!any(lapply(stan_data_joint, function(x) any(is.na(x))) %>% unlist))
# save data lists as json files for use in cmdstan for model execution
write_stan_json(data = stan_data_climate, file = './data/stan_lists/data_climate.json')
write_stan_json(data = stan_data_erc, file = './data/stan_lists/data_erc.json')
write_stan_json(data = stan_data_fwi, file = './data/stan_lists/data_fwi.json')
write_stan_json(data = stan_data_erc_fwi, file = './data/stan_lists/data_erc_fwi.json')
write_stan_json(data = stan_data_joint, file = './data/stan_lists/data_joint.json')
# save data lists as RDS files for use in result processing
saveRDS(stan_data_climate, file = './data/stan_lists/data_climate.RDS')
saveRDS(stan_data_erc, file = './data/stan_lists/data_erc.RDS')
saveRDS(stan_data_fwi, file = './data/stan_lists/data_fwi.RDS')
saveRDS(stan_data_erc_fwi, file = './data/stan_lists/data_erc_fwi.RDS')
saveRDS(stan_data_joint, file = './data/stan_lists/data_joint.RDS')
