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
  arrange(NA_L3NAME, ym) %>% filter(year < 2005)

# %>%
#   filter(year < 2005)
# REMOVE THIS LINE WHEN RUNNING ON FULL DATASET

assert_that(length(setdiff(st_covs$NA_L3NAME, count_df$NA_L3NAME)) == 0)
assert_that(!anyDuplicated(st_covs))
st_covs$id <- 1:nrow(st_covs)


# Create training sets, including years from 1984 to cutoff_year - 1
cutoff_year <- 2005
# remember to change "2016" (or "upper_cutoff" to whatever data is most current

# train_counts <- count_df %>%
#   filter(year < cutoff_year) %>%
#   left_join(st_covs)


# including missing data
train_burns_full <- burn_df %>%
  filter(FIRE_YEAR < cutoff_year) %>%
  right_join(st_covs) %>% arrange(NA_L3NAME, ym)
# 
# train_burns_obs <- burn_df %>%
#   filter(FIRE_YEAR < cutoff_year) %>%
#   left_join(st_covs) %>% arrange(NA_L3NAME, ym)

# upper_cutoff <- 2005 # REMOVE THIS LINE WHEN RUNNING ON FULL DATASET
# holdout_counts <- count_df %>%
#   filter(year >= cutoff_year & year <= upper_cutoff)
# 
# holdout_c_idx <- match(holdout_counts$er_ym, st_covs$er_ym)
# 
# holdout_burns_miss <-  burn_df %>%
#   filter(FIRE_YEAR >= cutoff_year & FIRE_YEAR <= upper_cutoff) %>%
#   left_join(st_covs)
# 
# holdout_burns_obs <-  burn_df %>%
#   filter(FIRE_YEAR >= cutoff_year & FIRE_YEAR <= upper_cutoff) %>%
#   left_join(st_covs)
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
X_full <- X_bs_df %>% mutate(er_ym = st_covs$er_ym, NA_L3NAME = st_covs$NA_L3NAME) %>% 
  mutate(intercept = 1, .before = lin_log_housing_density)

# design matrix for training counts ------
# is a subset of the rows of X, based on which rows show up in train_counts
count_idx_train <- match(train_counts$er_ym, X_full$er_ym)
assert_that(all(diff(count_idx_train)) == 1)
X_tc <- X_full[count_idx_train, ]
assert_that(identical(nrow(X_tc), nrow(train_counts)))
assert_that(all(X_tc$er_ym == train_counts$er_ym))

# split X into 84 matrices that are 192 x 37
X_list_tc <- lapply(split(X_tc, X_tc$NA_L3NAME), function(x) select(x, -NA_L3NAME))
nfire_list <- lapply(split(train_counts, train_counts$NA_L3NAME), function(x) select(x, c(n_fire, er_ym)))
assert_that(all(mapply(function(x, y) identical(x$er_ym, y$er_ym), X_list_tc, nfire_list)) == TRUE)
T_tc <- nrow(X_list_tc[[1]])
X_list_tc <- lapply(X_list_tc, function(x) select(x, -er_ym))
X_array_tc <- array(NA, dim = c(84, T_tc, 37))
for(i in 1:84) {
  X_array_tc[i, ,] <- as.matrix(X_list_tc[[i]])
}

nfire_matrix <- matrix(unlist(lapply(nfire_list, function(x) select(x, n_fire))), nrow(nfire_list[[1]]), 84)
iden_vec <- c()
for(i in 1:84) {
  iden_vec[i] <- all(nfire_matrix[, i] == nfire_list[[i]]$n_fire)
}
assert_that(all(iden_vec) == TRUE)

# ensure that split data frames are still in correct order
# assert_that(all(bind_rows(X_list_tc)$lin_log_housing_density == log(train_counts$housing_density)))
# assert_that(all(bind_rows(X_list_tc)$lin_vs == train_counts$vs))
# assert_that(all(bind_rows(X_list_tc)$lin_pr == train_counts$pr))
# assert_that(all(bind_rows(X_list_tc)$lin_vs == train_counts$vs))
# assert_that(all(bind_rows(X_list_tc)$lin_prev_12mo_precip == train_counts$prev_12mo_precip))
# assert_that(all(bind_rows(X_list_tc)$lin_tmmx == train_counts$tmmx))
# assert_that(all(bind_rows(X_list_tc)$lin_rmin == train_counts$rmin))
assert_that(all(bind_rows(X_list_tc) == X_full[count_idx_train,-c(38:39)]))
iden_vec <- c()
for(i in 1:84) {
  iden_vec[i] <- all(X_array_tc[i,,] == X_list_tc[[i]])
}
assert_that(all(iden_vec) == TRUE)

# st_covs_tc <- st_covs[count_idx_train, ]
# st_covs_tc_er <- split(st_covs_tc, st_covs_tc$NA_L3NAME)
# 
# count_idx_future <- setdiff(1:nrow(st_covs), count_idx_train)
# assert_that(all(st_covs_tc$er_ym == train_counts$er_ym))

# design matrix for training burn areas -------
# pull indices from train_burns_full for use in stan model
idx_burns_mis <- which(is.na(train_burns_full$BurnBndAc)) # indices of missing y
idx_burns_obs <- which(!is.na(train_burns_full$BurnBndAc)) # indices of observed y
idx_burns_all <- match(train_burns_full$er_ym, X_full$er_ym) # indices to broadcast kappa, sigma, and xi in the model for observed and missing y
assert_that(all(X_full[idx_burns_all,]$er_ym == train_burns_full$er_ym)) # check broadcasting works

burns_obs <- train_burns_full %>% filter(!is.na(BurnBndAc)) # dataframe of observed valeus

X_list_tb <- lapply(split(X_full, X_full$NA_L3NAME), function(x) select(x, -NA_L3NAME))
assert_that(all(bind_rows(X_list_tb)$er_ym == X_full$er_ym))

# reshape list of X's into array that is 84 x T_tb x 37
T_tb <- nrow(X_list_tb[[1]])
X_list_tb <- lapply(X_list_tb, function(x) select(x, -er_ym))
X_array_tb <- array(NA, dim = c(84, T_tb, 37))
for(i in 1:84) {
  X_array_tb[i, ,] <- as.matrix(X_list_tb[[i]])
}

compare_X <- X_full %>% select(-c(NA_L3NAME, er_ym))

# generate correlation matrix indicators
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
l3 <- level3
l2 <- level2
l1 <- level1

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

# standardize design matrix
std_data <- function(x) {
  if(sd(x) != 0) {
    return((x-mean(x))/sd(x)) # don't want to divide by zero for any columns where all values are the same
  } else {
    return(x)
  }
}

X_std_array <- array(NA, dim = c(84, T_tb, 37))
for(i in 1:84) {
  X_std_array[i, ,] <- apply(X_array_tb[i, ,], 2, std_data)
}


# Bundle up data into a list too pass to Stan -----------------------------
min_size <- 1e3
stan_data <- list(
  R = 84, # total number of regions
  T = T_tb,
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
  X = X_std_array,
  y_obs = burns_obs$BurnBndAc - min_size
)

# assert that there are no missing values in stan_d
assert_that(!any(lapply(stan_data, function(x) any(is.na(x))) %>% unlist))

saveRDS(stan_data, file = "full-model/simulations/g1/data/burns_sliced-index_std.RDS")

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
