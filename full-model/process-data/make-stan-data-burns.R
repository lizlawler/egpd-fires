library(tidyverse)
library(lubridate)
library(rstan)
library(splines)
library(spdep)
library(spatialreg)

source('./full-model/process-data/merge-data.R')

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

st_covs <- ecoregion_summaries %>%
  left_join(er_df) %>%
  filter(!NA_L2NAME == "UPPER GILA MOUNTAINS (?)",
         year > 1983,
         ym <= max(mtbs$ym)) %>%
  mutate(log_housing_density = log(housing_density),
         pr = ifelse(pr < 0 , 0, pr)) %>%
  left_join(area_df) %>%
  droplevels %>%
  mutate(er_ym = paste(NA_L3NAME, ym, sep = "_")) %>%
  arrange(NA_L3NAME, ym)

# %>%
#   filter(year < 2005)
# REMOVE THIS LINE WHEN RUNNING ON FULL DATASET

assert_that(length(setdiff(st_covs$NA_L3NAME, count_df$NA_L3NAME)) == 0)
assert_that(!anyDuplicated(st_covs))
st_covs$id <- 1:nrow(st_covs)


# Create training sets, including years from 1984 to cutoff_year - 1
cutoff_year <- 2005

# train_counts <- count_df %>%
#   filter(year < cutoff_year) %>%
#   left_join(st_covs)


# including missing data
train_burns_full <- burn_df %>%
  filter(FIRE_YEAR < cutoff_year) %>%
  right_join(st_covs %>% filter(year < cutoff_year)) %>% arrange(NA_L3NAME, ym)
# 
# train_burns_obs <- burn_df %>%
#   filter(FIRE_YEAR < cutoff_year) %>%
#   left_join(st_covs) %>% arrange(NA_L3NAME, ym)

holdout_burns_full <-  burn_df %>%
  filter(FIRE_YEAR >= cutoff_year) %>%
  right_join(st_covs %>% filter(year >= cutoff_year)) %>% arrange(NA_L3NAME, ym)
# holdout_burns_miss <-  burn_df %>%
#   filter(FIRE_YEAR >= cutoff_year & FIRE_YEAR <= upper_cutoff) %>%
#   left_join(st_covs)
# 

# 
# holdout_b_idx <- match(holdout_burns$er_ym, holdout_counts$er_ym)

# this data frame has no duplicate ecoregion X timestep combos
N <- length(unique(st_covs$NA_L3NAME))
T <- length(unique(st_covs$ym))

assert_that(identical(nrow(st_covs), N * T))

# Create b-splines for climate vars
vars <- c('log_housing_density', 'vs',
          'pr', 'prev_12mo_precip', 'tmmx',
          'rmin')

df_each <- 5
deg_each <- 3
X_bs <- list()
X_bs_df <- list()
for (i in seq_along(vars)) {
  varname <- paste("lin", vars[i], sep = "_")
  X_bs[[i]] <- bs(x = st_covs[[vars[i]]], df = df_each, degree = deg_each, 
                  Boundary.knots = range(st_covs[[vars[i]]]), intercept = FALSE)

  X_bs_df[[i]] <- X_bs[[i]] %>% as_tibble()
  names(X_bs_df[[i]]) <- paste('bs', vars[[i]], 1:df_each, sep = '_')
  X_bs_df[[i]] <- X_bs_df[[i]] %>%
    mutate(!!varname := st_covs[[vars[i]]]) %>%
    relocate(!!varname, before = where(is.character))
}
X_bs_df <- bind_cols(X_bs_df)
assert_that(!any(is.na(X_bs_df)))

# X_full <- X_bs_df %>% mutate(NA_L3NAME = st_covs$NA_L3NAME) %>% mutate(intercept = 1, .before = lin_log_housing_density)
X_full <- X_bs_df %>% mutate(er_ym = st_covs$er_ym, NA_L3NAME = st_covs$NA_L3NAME, year = st_covs$year) %>% 
  mutate(intercept = 1, .before = lin_log_housing_density)
X_tb <- X_full %>% filter(year < 2005) %>% select(-year)

# design matrix for training burn areas -------
# pull indices from train_burns_full for use in stan model
idx_burns_tb_mis <- which(is.na(train_burns_full$BurnBndAc)) # indices of missing y
idx_burns_tb_obs <- which(!is.na(train_burns_full$BurnBndAc)) # indices of observed y; subset of all of the rows in the training dataframe
idx_burns_tb_all <- match(train_burns_full$er_ym, X_tb$er_ym) # indices to broadcast kappa, sigma, and xi in the model for observed and missing y
assert_that(all(X_tb[idx_burns_tb_all, ]$er_ym == train_burns_full$er_ym)) # check broadcasting works
assert_that(all(X_tb[idx_burns_tb_all, ][idx_burns_tb_obs, ]$er_ym == train_burns_full[idx_burns_tb_obs,]$er_ym))

# pull indices for holdout_burns_full for use in log scores
idx_burns_hold_obs <- which(!is.na(holdout_burns_full$BurnBndAc))
idx_burns_hold_all <- match(holdout_burns_full$er_ym, X_full$er_ym) # indices to broadcast kappa, sigma, and xi in the model for observed and missing y
assert_that(all(X_full[idx_burns_hold_all, ]$er_ym == holdout_burns_full$er_ym)) # check broadcasting works
assert_that(all(X_full[idx_burns_hold_all, ][idx_burns_hold_obs, ]$er_ym == holdout_burns_full[idx_burns_hold_obs,]$er_ym))

# use square root of burn area so MCMC chains mix 
burn_train_obs <- sqrt(train_burns_full$BurnBndAc[idx_burns_tb_obs])
assert_that(all(!is.na(burn_train_obs)))
hist(burn_train_obs)
burn_hold_obs <- sqrt(holdout_burns_full$BurnBndAc[idx_burns_hold_obs])
assert_that(all(!is.na(burn_hold_obs)))
hist(burn_hold_obs)

# function to standardize design matrices
std_data <- function(x) {
  if(sd(x) != 0) {
    return((x-mean(x))/sd(x)) # don't want to divide by zero for any columns where all values are the same
  } else {
    return(x)
  }
}

# first, split X_full into 84 design matrices
# then, standardize each regions design matrix
X_list_full <- lapply(split(X_full, X_full$NA_L3NAME), function(x) select(x, -NA_L3NAME))
assert_that(all(bind_rows(X_list_full)$er_ym == X_full$er_ym))
design_train_idx <- which(X_list_full[[1]]$year < 2005)
X_list_full <- lapply(X_list_full, function(x) select(x, -c(year, er_ym)))
X_list_full_std <- lapply(X_list_full, function(df) apply(df, 2, std_data))

# last, split into training and also keep full (for predictive effect in log scores)
X_list_tb <- lapply(X_list_full_std, function(x) x[design_train_idx,])

# reshape each for use in the stan model
t_tb <- nrow(X_list_tb[[1]])
t_all <- nrow(X_list_full_std[[1]])
X_array_tb <- array(NA, dim = c(84, t_tb, 37))
X_array_full <- array(NA, dim = c(84, t_all, 37))
for(i in 1:84) {
  X_array_tb[i, ,] <- as.matrix(X_list_tb[[i]])
  X_array_full[i, ,] <- as.matrix(X_list_full_std[[i]])
}

# generate correlation matrix indicator matrices (three total, one for each level of ecoregion)
# create correlation matrix from 3 levels of relationships using real ecoregions
load(file = "./sim-study/shared-data/region_key.RData")

full_reg_key <- as_tibble(region_key) %>% 
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

# generate AR(1) indicator matrices for use in covariance matrix in stan
p <- ncol(X_array_tb[1, ,])
zeroes <- matrix(0, p, p)
equal <- diag(p)
bp_lin <- zeroes
bp_square <- zeroes
bp_cube <- zeroes
bp_quart <- zeroes
cov_vec_idx <- c(3:7, 9:13, 15:19, 21:25, 27:31, 33:37) # 1st element is global intercept; 2, 8, 14, 20, 26, and 32 are the linear terms of each of the 6 covariates

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

# Bundle up data into a list too pass to Stan -----------------------------
stan_data <- list(
  r = 84, # total number of regions
  t_all = t_all,
  p = p,
  
  # training data
  X_tb = X_array_tb,
  y_tb_obs = burn_train_obs,
  N_tb_obs = length(idx_burns_tb_obs),
  N_tb_mis = length(idx_burns_tb_mis),
  N_tb_all = length(idx_burns_tb_all),
  ii_tb_obs = idx_burns_tb_obs,
  ii_tb_mis = idx_burns_tb_mis,
  ii_tb_all = idx_burns_tb_all, # for broadcasting params in likelihood

  # for predicting effects and assessing holdout fit
  X_full = X_array_full,
  y_hold_obs = burn_hold_obs,
  N_hold_obs = length(idx_burns_hold_obs),
  N_tb_all = length(idx_burns_tb_all),
  ii_tb_obs = idx_burns_tb_obs,
  ii_tb_mis = idx_burns_tb_mis,
  ii_tb_all = idx_burns_tb_all, # for broadcasting params in likelihood
  
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
  
  # training data
  X = X_std_array,
  y_obs = burn_vec
)

# assert that there are no missing values in stan_d
assert_that(!any(lapply(stan_data, function(x) any(is.na(x))) %>% unlist))

saveRDS(stan_data, file = "full-model/fire-sims/burns/data/burns_sliced-index_Xstd_ystd.RDS")

zi_d <- stan_d
zi_d$M <- 2
write_rds(zi_d, 'data/processed/zi_d.rds')

write_rds(stan_d, file = 'data/processed/stan_d.rds')
print('stan_d.rds written!')


write_rds(st_covs, 'data/processed/st_covs.rds')
write_rds(cutoff_year, 'data/processed/cutoff_year.rds')
write_rds(train_counts, 'data/processed/train_counts.rds')
write_rds(holdout_counts, 'data/processed/holdout_counts.rds')
write_rds(train_burns, 'data/processed/train_burns.rds')
write_rds(holdout_burns, 'data/processed/holdout_burns.rds')
write_rds(colnamesX, 'data/processed/colnamesX.rds')
write_rds(X, 'data/processed/X.rds')
write_rds(ecoregions, 'data/processed/ecoregions.rds')
write_rds(vars, 'data/processed/vars.rds')
write_rds(mtbs, 'data/processed/mtbs.rds')
