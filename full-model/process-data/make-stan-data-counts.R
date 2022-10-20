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

assert_that(length(setdiff(st_covs$NA_L3NAME, count_df$NA_L3NAME)) == 0)
assert_that(!anyDuplicated(st_covs))
st_covs$id <- 1:nrow(st_covs)


# Create training sets, including years from 1984 to cutoff_year - 1
cutoff_year <- 2006
# remember to change "2016" (or "upper_cutoff" to whatever data is most current

train_counts <- count_df %>%
  filter(year < cutoff_year) %>%
  left_join(st_covs)

holdout_counts <- count_df %>%
  filter(year >= cutoff_year)

idx_count_ho <- match(holdout_counts$er_ym, st_covs$er_ym)
assert_that(all(holdout_counts$er_ym == st_covs[holdout_c_idx, ]$er_ym))

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
idx_count_train <- match(train_counts$er_ym, X_full$er_ym)
assert_that(all(diff(idx_count_train)) == 1)
X_tc <- X_full[idx_count_train, ]
assert_that(identical(nrow(X_tc), nrow(train_counts)))
assert_that(all(X_tc$er_ym == train_counts$er_ym))
assert_that(all(X_full[holdout_c_idx,]$er_ym == holdout_counts$er_ym))

# split training X into 84 matrices that are t_tc x 37
X_list_tc <- lapply(split(X_tc, X_tc$NA_L3NAME), function(x) select(x, -NA_L3NAME))
nfire_list_tc <- lapply(split(train_counts, train_counts$NA_L3NAME), function(x) select(x, c(n_fire, er_ym)))
assert_that(all(mapply(function(x, y) identical(x$er_ym, y$er_ym), X_list_tc, nfire_list_tc)) == TRUE)
t_tc <- nrow(X_list_tc[[1]])
X_list_tc <- lapply(X_list_tc, function(x) select(x, -er_ym))
X_array_tc <- array(NA, dim = c(84, t_tc, 37))
for(i in 1:84) {
  X_array_tc[i, ,] <- as.matrix(X_list_tc[[i]])
}
assert_that(all(bind_rows(X_list_tc) == X_full[idx_count_train,-c(38:39)]))

nfire_matrix_tc <- matrix(unlist(lapply(nfire_list_tc, function(x) select(x, n_fire))), nrow(nfire_list_tc[[1]]), 84)
iden_vec <- c()
for(i in 1:84) {
  iden_vec[i] <- all(nfire_matrix_tc[, i] == nfire_list_tc[[i]]$n_fire)
}
assert_that(all(iden_vec) == TRUE)

iden_vec <- c()
for(i in 1:84) {
  iden_vec[i] <- all(X_array_tc[i,,] == X_list_tc[[i]])
}
assert_that(all(iden_vec) == TRUE)

# split full X matrix into 84 matrices that are t_all x 37
X_list_full <- lapply(split(X_full, X_full$NA_L3NAME), function(x) select(x, -NA_L3NAME))
t_all <- nrow(X_list_full[[1]])
X_list_full <- lapply(X_list_full, function(x) select(x, -er_ym))
X_array_full <- array(NA, dim = c(84, t_all, 37))
for(i in 1:84) {
  X_array_full[i, ,] <- as.matrix(X_list_full[[i]])
}
assert_that(all(bind_rows(X_list_full) == X_full[,-c(38:39)]))
idx_tc_by_er <- 1:t_tc
assert_that(all(X_list_tc[[1]] == X_list_full[[1]][idx_tc_by_er,]))
idx_holdout_by_er <- setdiff(1:t_all, idx_tc_by_er)

# grab holdout counts
X_list_ho <- lapply(split(X_full[idx_count_ho,], X_full[idx_count_ho,]$NA_L3NAME), function(x) select(x, -NA_L3NAME))
nfire_list_ho <- lapply(split(holdout_counts, holdout_counts$NA_L3NAME), function(x) select(x, c(n_fire, er_ym)))
assert_that(all(mapply(function(x, y) identical(x$er_ym, y$er_ym), X_list_ho, nfire_list_ho)) == TRUE)
nfire_matrix_ho <- matrix(unlist(lapply(nfire_list_ho, function(x) select(x, n_fire))), nrow(nfire_list_ho[[1]]), 84)

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
p <- ncol(X_array_tc[1, ,])
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

X_std_array_tc <- array(NA, dim = c(84, t_tc, 37))
X_std_array_full <- array(NA, dim = c(84, t_all, 37))
for(i in 1:84) {
  X_std_array_tc[i, ,] <- apply(X_array_tc[i, ,], 2, std_data)
  X_std_array_full[i, ,] <- apply(X_array_full[i, ,], 2, std_data)
}



# Bundle up data into a list too pass to Stan -----------------------------
# min_size <- 1e3
stan_data <- list(
  r = 84, # total number of regions
  p = p,
  
  # all data
  t_all = t_all,
  X_all_tmpt = X_std_array_full,
  area_offset = log(area_df$area * 1e-11) / 2,
  
  # training data
  t_tc = t_tc,
  X_tc = X_std_array,
  y_tc = nfire_matrix_tc,
  idx_tc_er = idx_tc_by_er,

  # holdout data
  t_ho = t_all - t_tc,
  y_hold = nfire_matrix_ho,
  idx_hold_er = idx_holdout_by_er,
  
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
  node2 = B@j + 1
)

# assert that there are no missing values in stan_d
assert_that(!any(lapply(stan_data, function(x) any(is.na(x))) %>% unlist))

saveRDS(stan_data, file = "full-model/fire-sims/counts/data/stan_data_train-hold_counts.RDS")

zi_d <- stan_d
zi_d$M <- 2
write_rds(zi_d, 'data/processed/zi_d.rds')

write_rds(stan_d, file = 'data/processed/stan_d.rds')
print('stan_d.rds written!')

