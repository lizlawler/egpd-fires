library(tidyverse)
library(cowplot)
library(tidyverse)
library(patchwork)
library(ggrepel)
library(rstan)
# 
# model_fits <- list.files(pattern = "*.RDS", recursive = TRUE)
# post_params_files <- grep(model_fits, pattern = 'post_params', value = TRUE, invert = FALSE)
# models <- c("g1_ogY", "g1_sqrtY", "zinb_byER", "zinb", "zip")
# post_params_list <- vector("list", 5)
# names(post_params_list) <- models
# for(i in seq_along(post_params_files)) {
#   post_params_list[[i]] <- readRDS(post_params_files[i])
# }
# counts <- models[3:5]
# burns <- models[1:2]
# 
# holdout_loglik_counts <- vector("list", 3)
# train_loglik_counts <- vector("list", 3)
# for(i in seq_along(counts)) {
#   holdout_loglik_counts[[i]] <- post_params_list[[counts[i]]]$holdout_loglik %>%
#     apply(., 1, c) %>%
#     as_tibble() %>%
#     pivot_longer(cols = everything(), names_to = "iter") %>%
#     mutate(iter = as.numeric(gsub("V", "", iter))) %>%
#     group_by(iter) %>%
#     summarize(value = sum(value)) %>%
#     mutate(model = counts[i], train = FALSE)
#   train_loglik_counts[[i]] <- post_params_list[[counts[i]]]$train_loglik %>%
#     apply(., 1, c) %>%
#     as_tibble() %>%
#     pivot_longer(cols = everything(), names_to = "iter") %>%
#     mutate(iter = as.numeric(gsub("V", "", iter))) %>%
#     group_by(iter) %>%
#     summarize(value = sum(value)) %>%
#     mutate(model = counts[i], train = TRUE)
# }
# 
# holdout_loglik_c <- bind_rows(holdout_loglik_counts)
# train_loglik_c <- bind_rows(train_loglik_counts)
# 
# ll_full <- holdout_loglik_c %>%
#   full_join(train_loglik_c) %>%
#   mutate(train = ifelse(train == TRUE, 'train', 'test')) %>%
#   pivot_wider(names_from = train, values_from = value)
# 
# ll_test <- ll_full %>%
#   group_by(model) %>%
#   summarize(mean_test = mean(test),
#             sd_test = sd(test)) %>%
#   arrange(-mean_test) 
# ll_train <- ll_full %>%
#   group_by(model) %>%
#   summarize(mean_train = mean(train),
#             sd_train = sd(train)) %>%
#   arrange(-mean_train) 

stan_data <- readRDS("full-model/data/stan_data_sqrt.RDS")
X <- stan_data$X_train
vars <- c('log_housing_density', 'vs',
          'pr', 'prev_12mo_precip', 'tmmx',
          'rmin')
X_covar <- c()
X_cols <- vector("list", length(vars))
start <- 2
for(i in seq_along(vars)) {
  X_covar[i] <- paste0("X_", vars[i])
  X_cols[[i]] <- c(1, start:(start+5))
  start = start + 6
}

load(file = "./full-model/data/processed/region_key.RData")
full_reg_key <- as_tibble(region_key) %>% 
  mutate(region = c(1:84),
         NA_L2CODE = as.factor(NA_L2CODE),
         NA_L1CODE = as.factor(NA_L1CODE),
         NA_L3CODE = as.factor(NA_L3CODE))
reg_cols <- full_reg_key$region
g1_nocens_kappa <- extract(egpd_fit, pars = "beta_kappa")[[1]] %>% apply(., c(2,3), median)
g1_nocens_nu <- extract(egpd_fit, pars = "beta_nu")[[1]] %>% apply(., c(2,3), median)
g1_nocens_xi <- extract(egpd_fit, pars = "beta_xi")[[1]] %>% apply(., c(2,3), median)

# zinb_byER_lambda <- apply(post_params_list$zinb_byER$beta_lambda, c(2,3), median)
# zinb_lambda <- apply(post_params_list$zinb$beta_lambda, c(2,3), median)

r <- 84
t <- nrow(X[1,,])
# zip_lambda_effects_list <- vector("list", 6)
# for(i in seq_along(vars)) {
#   zip_lambda_effects_df <- matrix(NA, t, r)
#   for(j in 1:r) {
#     zip_lambda_effects_df[, j] <- X[j, , X_cols[[i]]] %*% zip_lambda[X_cols[[i]], j]
#   }
#   zip_lambda_effects_list[[i]] <- zip_lambda_effects_df %>% 
#     as_tibble() %>% 
#     rename_with(., ~ reg_cols) %>% 
#     mutate(time = c(1:t)) %>% 
#     pivot_longer(cols = c(1:all_of(r)), values_to = "effect", names_to = "region") %>%
#         mutate(region = as.numeric(region), covar = vars[i])
# }
# X_long <- list()
# for(i in 1:r) {
#   X_long[[i]] <- X[i,,] %>% as_tibble() %>% 
#     mutate(region = reg_cols[i], 
#            time = c(1:t))
# }

g1_nocens_kappa_effects_list <- vector("list", 6)
for(i in seq_along(vars)) {
  g1_nocens_kappa_effects_df <- matrix(NA, t, r)
  for(j in 1:r) {
    g1_nocens_kappa_effects_df[, j] <- X[j, , X_cols[[i]]] %*% g1_nocens_kappa[X_cols[[i]], j]
  }
  g1_nocens_kappa_effects_list[[i]] <- g1_nocens_kappa_effects_df %>% 
    as_tibble() %>% 
    rename_with(., ~ reg_cols) %>% 
    mutate(time = c(1:t)) %>% 
    pivot_longer(cols = c(1:all_of(r)), values_to = "effect", names_to = "region") %>%
    mutate(region = as.numeric(region), covar = vars[i])
}
X_long <- list()
for(i in 1:r) {
  X_long[[i]] <- X[i,,] %>% as_tibble() %>% 
    mutate(region = reg_cols[i], 
           time = c(1:t))
}

X_long_v1 <- X_long %>% bind_rows() %>% select(c(V2, region, time)) %>% rename(linear = V2) %>% mutate(covar = vars[1])
X_long_v2 <- X_long %>% bind_rows() %>% select(c(V8, region, time)) %>% rename(linear = V8) %>% mutate(covar = vars[2])
X_long_v3 <- X_long %>% bind_rows() %>% select(c(V14, region, time)) %>% rename(linear = V14) %>% mutate(covar = vars[3])
X_long_v4 <- X_long %>% bind_rows() %>% select(c(V20, region, time)) %>% rename(linear = V20) %>% mutate(covar = vars[4])
X_long_v5 <- X_long %>% bind_rows() %>% select(c(V26, region, time)) %>% rename(linear = V26) %>% mutate(covar = vars[5])
X_long_v6 <- X_long %>% bind_rows() %>% select(c(V32, region, time)) %>% rename(linear = V32) %>% mutate(covar = vars[6])

kappa_effects_v1 <- g1_nocens_kappa_effects_list[[1]] %>% left_join(., X_long_v1) %>% left_join(., full_reg_key)
# kappa_housing <- ggplot(kappa_effects_v1, aes(x=linear, y = effect, group = region)) + 
#   geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE))
# ggsave("full-model/figures/g1/effects/nocens_kappa_housing.png", dpi = 320, type = "cairo")

kappa_effects_v2 <- g1_nocens_kappa_effects_list[[2]] %>% left_join(., X_long_v2) %>% left_join(., full_reg_key)
# kappa_wind <- ggplot(kappa_effects_v2, aes(x=linear, y = effect, group = region)) + 
#   geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE))
# ggsave("full-model/figures/g1/effects/nocens_kappa_wind.png", dpi = 320, type = "cairo")

kappa_effects_v3 <- g1_nocens_kappa_effects_list[[3]] %>% left_join(., X_long_v3) %>% left_join(., full_reg_key)
# kappa_rain <- ggplot(kappa_effects_v3, aes(x=linear, y = effect, group = region)) + 
#   geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE))
# ggsave("full-model/figures/g1/effects/nocens_kappa_rain.png", dpi = 320, type = "cairo")

kappa_effects_v4 <- g1_nocens_kappa_effects_list[[4]] %>% left_join(., X_long_v4) %>% left_join(., full_reg_key)
# kappa_rain12mo <- ggplot(kappa_effects_v4, aes(x=linear, y = effect, group = region)) + 
#   geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE))
# ggsave("full-model/figures/g1/effects/nocens_kappa_rain12mo.png", dpi = 320, type = "cairo")

kappa_effects_v5 <- g1_nocens_kappa_effects_list[[5]] %>% left_join(., X_long_v5) %>% left_join(., full_reg_key)
# kappa_maxtemp <- ggplot(kappa_effects_v5, aes(x=linear, y = effect, group = region)) + 
#   geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE))
# ggsave("full-model/figures/g1/effects/nocens_kappa_maxtemp.png", dpi = 320, type = "cairo")

kappa_effects_v6 <- g1_nocens_kappa_effects_list[[6]] %>% left_join(., X_long_v6) %>% left_join(., full_reg_key)
# kappa_relhumid <- ggplot(kappa_effects_v6, aes(x=linear, y = effect, group = region)) + 
#   geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE))
# ggsave("full-model/figures/g1/effects/nocens_kappa_relhumid.png", dpi = 320, type = "cairo")

full_nocens_kappa <- rbind(kappa_effects_v1, kappa_effects_v2, kappa_effects_v3, kappa_effects_v4, kappa_effects_v5, kappa_effects_v6)
nocens_kappa_grid <- ggplot(full_nocens_kappa, aes(x = linear, y = effect, group = region)) + geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE)) + facet_wrap(. ~ covar, scales = "free_x")
ggsave("full-model/figures/g1/effects/nocens_kappa_allvars.png", dpi = 320, type = "cairo")


g1_nocens_nu_effects_list <- vector("list", 6)
for(i in seq_along(vars)) {
  g1_nocens_nu_effects_df <- matrix(NA, t, r)
  for(j in 1:r) {
    g1_nocens_nu_effects_df[, j] <- X[j, , X_cols[[i]]] %*% g1_nocens_nu[X_cols[[i]], j]
  }
  g1_nocens_nu_effects_list[[i]] <- g1_nocens_nu_effects_df %>% 
    as_tibble() %>% 
    rename_with(., ~ reg_cols) %>% 
    mutate(time = c(1:t)) %>% 
    pivot_longer(cols = c(1:all_of(r)), values_to = "effect", names_to = "region") %>%
    mutate(region = as.numeric(region), covar = vars[i])
}

nu_effects_v1 <- g1_nocens_nu_effects_list[[1]] %>% left_join(., X_long_v1) %>% left_join(., full_reg_key)
nu_effects_v2 <- g1_nocens_nu_effects_list[[2]] %>% left_join(., X_long_v2) %>% left_join(., full_reg_key)
nu_effects_v3 <- g1_nocens_nu_effects_list[[3]] %>% left_join(., X_long_v3) %>% left_join(., full_reg_key)
nu_effects_v4 <- g1_nocens_nu_effects_list[[4]] %>% left_join(., X_long_v4) %>% left_join(., full_reg_key)
nu_effects_v5 <- g1_nocens_nu_effects_list[[5]] %>% left_join(., X_long_v5) %>% left_join(., full_reg_key)
nu_effects_v6 <- g1_nocens_nu_effects_list[[6]] %>% left_join(., X_long_v6) %>% left_join(., full_reg_key)

full_nocens_nu <- rbind(nu_effects_v1, nu_effects_v2, nu_effects_v3, nu_effects_v4, nu_effects_v5, nu_effects_v6)
nocens_nu_grid <- ggplot(full_nocens_nu, aes(x = linear, y = effect, group = region)) + geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE)) + facet_wrap(. ~ covar, scales = "free_x")
ggsave("full-model/figures/g1/effects/nocens_nu_allvars.png", dpi = 320, type = "cairo")

g1_nocens_xi_effects_list <- vector("list", 6)
for(i in seq_along(vars)) {
  g1_nocens_xi_effects_df <- matrix(NA, t, r)
  for(j in 1:r) {
    g1_nocens_xi_effects_df[, j] <- X[j, , X_cols[[i]]] %*% g1_nocens_xi[X_cols[[i]], j]
  }
  g1_nocens_xi_effects_list[[i]] <- g1_nocens_xi_effects_df %>% 
    as_tibble() %>% 
    rename_with(., ~ reg_cols) %>% 
    mutate(time = c(1:t)) %>% 
    pivot_longer(cols = c(1:all_of(r)), values_to = "effect", names_to = "region") %>%
    mutate(region = as.numeric(region), covar = vars[i])
}

xi_effects_v1 <- g1_nocens_xi_effects_list[[1]] %>% left_join(., X_long_v1) %>% left_join(., full_reg_key)
xi_effects_v2 <- g1_nocens_xi_effects_list[[2]] %>% left_join(., X_long_v2) %>% left_join(., full_reg_key)
xi_effects_v3 <- g1_nocens_xi_effects_list[[3]] %>% left_join(., X_long_v3) %>% left_join(., full_reg_key)
xi_effects_v4 <- g1_nocens_xi_effects_list[[4]] %>% left_join(., X_long_v4) %>% left_join(., full_reg_key)
xi_effects_v5 <- g1_nocens_xi_effects_list[[5]] %>% left_join(., X_long_v5) %>% left_join(., full_reg_key)
xi_effects_v6 <- g1_nocens_xi_effects_list[[6]] %>% left_join(., X_long_v6) %>% left_join(., full_reg_key)

full_nocens_xi <- rbind(xi_effects_v1, xi_effects_v2, xi_effects_v3, xi_effects_v4, xi_effects_v5, xi_effects_v6)
nocens_xi_grid <- ggplot(full_nocens_xi, aes(x = linear, y = effect, group = region)) + geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE)) + facet_wrap(. ~ covar, scales = "free_x")
ggsave("full-model/figures/g1/effects/nocens_xi_allvars.png", dpi = 320, type = "cairo")

# grab data and pull relevant parameters
areas <- stan_data$y_train_obs
areas_test <- stan_data$y_hold_obs
areas_idx_obs <- stan_data$ii_tb_obs
areas_idx_all <- stan_data$ii_tb_all
kappa_train_obs <- extract(egpd_fit, pars = "kappa_train")[[1]]
sigma_train_obs <- extract(egpd_fit, pars = "sigma_train")[[1]]
xi_train_obs <- extract(egpd_fit, pars = "xi_train")[[1]]

f1_cdf <- function(x, sigma = sigma, xi = xi, kappa = kappa) {
  (1 - (1 + xi * (x/sigma))^(-1/xi))^kappa
}

f1_pdf <- function(x, sigma = sigma, xi = xi, kappa = kappa) {
  lpdf <- log(kappa) - log(sigma) - (1/xi + 1) * log(1 + xi * (x/sigma)) + 
    (kappa-1) * log(1 - (1 + xi * (x/sigma))^(-1/xi))
  return(exp(lpdf))
}

g1_pdf <- function(x, sigma = sigma, xi = xi, kappa = kappa) {
  lower <- f1_cdf(1.001, sigma, xi, kappa)
  return(f1_pdf(x, sigma, xi, kappa)/(1-lower))
}

g1_cdf <- function(x, sigma = sigma, xi = xi, kappa = kappa) {
  lower <- f1_cdf(1.001, sigma, xi, kappa)
  return((f1_cdf(x, sigma, xi, kappa) - lower)/(1-lower))
}

g1_cdf_inv <- function(u, sigma = sigma, xi = xi, kappa = kappa) {
  lower <- f1_cdf(1.001, sigma, xi, kappa)
  u_adj <- u * (1-lower) + lower
  (sigma/xi) * ((1-u_adj^(1/kappa))^-xi - 1)
}

g1_rng <- function(n, sigma, xi, kappa) {
  u = runif(n)
  return(g1_cdf_inv(u, sigma, xi, kappa))
}

max_space <- burn_df %>% group_by(NA_L3NAME) %>% summarize(burn = sqrt(max(BurnBndAc)/1000))

wt_fcn <- function(pred) {
  lpred <- length(pred)
  if(lpred > 1) {
    wt_vect <- rep(NA, lpred)
    for(i in 1:lpred) {
      wt_vect[i] <- pnorm(pred[i], mean = 9.8, sd = 7.3)
    } 
    return(wt_vect)
  } else {
    return(pnorm(pred, mean = 9.8, sd = 7.3))
  }
}
preds <- matrix(NA, nrow = 8018, ncol = 1000)
ecdf <- matrix(NA, nrow = 8018, ncol = 1000)
for(i in 1:8018) {
  preds[i,] <- g1_rng(1000, sigma_train_obs[i], xi_train_obs[i], kappa_train_obs[i])
  for(j in 1:1000) {
    ecdf[i,j] <- g1_cdf(preds[i,j], sigma_train_obs[i], xi_train_obs[i], kappa_train_obs[i])
  }
}
g1_cdf(preds[1,1], sigma_train_obs[1], xi_train_obs[1], kappa_train_obs[1])
wts <- wt_fcn(preds[1,])
sum((ecdf[1,] - (areas[1] <= preds[1,]))^2*wts)
sum((preds[1,] - (areas[1] <= preds[1,]))^2 * wts)/ 1000 
k <- 1:1000
sum((ecdf[1,k] - (areas[1] <= preds[1,k]))^2 * wts[k])
# zinb_byER_lambda_effects_list <- vector("list", 6)
# for(i in seq_along(vars)) {
#   zinb_byER_lambda_effects_df <- matrix(NA, t, r)
#   for(j in 1:r) {
#     zinb_byER_lambda_effects_df[, j] <- X[j, , X_cols[[i]]] %*% zinb_byER_lambda[X_cols[[i]], j]
#   }
#   zinb_byER_lambda_effects_list[[i]] <- zinb_byER_lambda_effects_df %>% 
#     as_tibble() %>% 
#     rename_with(., ~ reg_cols) %>% 
#     mutate(time = c(1:t)) %>% 
#     pivot_longer(cols = c(1:all_of(r)), values_to = "effect", names_to = "region") %>%
#     mutate(region = as.numeric(region), covar = vars[i])
# }
# 
# lambda_effects_zinbER_v1 <- zinb_byER_lambda_effects_list[[1]] %>% left_join(., X_long_v1) %>% left_join(., full_reg_key)
# ggplot(lambda_effects_zinbER_v1, aes(x=linear, y = effect, group = region)) +
#   geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE))
# 
# lambda_effects_zinbER_v2 <- zinb_byER_lambda_effects_list[[2]] %>% left_join(., X_long_v2) %>% left_join(., full_reg_key)
# ggplot(lambda_effects_zinbER_v2, aes(x=linear, y = effect, group = region)) +
#   geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE))
# 
# lambda_effects_zinbER_v3 <- zinb_byER_lambda_effects_list[[3]] %>% left_join(., X_long_v3) %>% left_join(., full_reg_key)
# ggplot(lambda_effects_zinbER_v3, aes(x=linear, y = effect, group = region)) +
#   geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE))
# 
# lambda_effects_zinbER_v4 <- zinb_byER_lambda_effects_list[[4]] %>% left_join(., X_long_v4) %>% left_join(., full_reg_key)
# ggplot(lambda_effects_zinbER_v4, aes(x=linear, y = effect, group = region)) +
#   geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE))
# 
# lambda_effects_zinbER_v5 <- zinb_byER_lambda_effects_list[[5]] %>% left_join(., X_long_v5) %>% left_join(., full_reg_key)
# ggplot(lambda_effects_zinbER_v5, aes(x=linear, y = effect, group = region)) +
#   geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE))
# 
# lambda_effects_zinbER_v6 <- zinb_byER_lambda_effects_list[[6]] %>% left_join(., X_long_v6) %>% left_join(., full_reg_key)
# ggplot(lambda_effects_zinbER_v6, aes(x=linear, y = effect, group = region)) +
#   geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE))
# 
# 
# zinb_lambda_effects_list <- vector("list", 6)
# for(i in seq_along(vars)) {
#   zinb_lambda_effects_df <- matrix(NA, t, r)
#   for(j in 1:r) {
#     zinb_lambda_effects_df[, j] <- X[j, , X_cols[[i]]] %*% zinb_lambda[X_cols[[i]], j]
#   }
#   zinb_lambda_effects_list[[i]] <- zinb_lambda_effects_df %>% 
#     as_tibble() %>% 
#     rename_with(., ~ reg_cols) %>% 
#     mutate(time = c(1:t)) %>% 
#     pivot_longer(cols = c(1:all_of(r)), values_to = "effect", names_to = "region") %>%
#     mutate(region = as.numeric(region), covar = vars[i])
# }
# 
# lambda_effects_zinb_v1 <- zinb_lambda_effects_list[[1]] %>% left_join(., X_long_v1) %>% left_join(., full_reg_key)
# ggplot(lambda_effects_zinb_v1, aes(x=linear, y = effect, group = region)) +
#   geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE))
# 
# lambda_effects_zinb_v2 <- zinb_lambda_effects_list[[2]] %>% left_join(., X_long_v2) %>% left_join(., full_reg_key)
# ggplot(lambda_effects_zinb_v2, aes(x=linear, y = effect, group = region)) +
#   geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE))
# 
# lambda_effects_zinb_v3 <- zinb_lambda_effects_list[[3]] %>% left_join(., X_long_v3) %>% left_join(., full_reg_key)
# ggplot(lambda_effects_zinb_v3, aes(x=linear, y = effect, group = region)) +
#   geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE))
# 
# lambda_effects_zinb_v4 <- zinb_lambda_effects_list[[4]] %>% left_join(., X_long_v4) %>% left_join(., full_reg_key)
# ggplot(lambda_effects_zinb_v4, aes(x=linear, y = effect, group = region)) +
#   geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE))
# 
# lambda_effects_zinb_v5 <- zinb_lambda_effects_list[[5]] %>% left_join(., X_long_v5) %>% left_join(., full_reg_key)
# ggplot(lambda_effects_zinb_v5, aes(x=linear, y = effect, group = region)) +
#   geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE))
# 
# lambda_effects_zinb_v6 <- zinb_lambda_effects_list[[6]] %>% left_join(., X_long_v6) %>% left_join(., full_reg_key)
# ggplot(lambda_effects_zinb_v6, aes(x=linear, y = effect, group = region)) +
#   geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE))
# 
# g1_sqrtY_kappa <- apply(post_params_list$g1_sqrtY$beta_kappa, c(2,3), median)
# g1_ogY_kappa <- apply(post_params_list$g1_ogY$beta_kappa, c(2,3), median)
# 
# g1_sqrtY_effects_list <- vector("list", 6)
# for(i in seq_along(vars)) {
#   g1_sqrtY_effects_df <- matrix(NA, t, r)
#   for(j in 1:r) {
#     g1_sqrtY_effects_df[, j] <- X[j, , X_cols[[i]]] %*% g1_sqrtY_kappa[X_cols[[i]], j]
#   }
#   g1_sqrtY_effects_list[[i]] <- g1_sqrtY_effects_df %>% 
#     as_tibble() %>% 
#     rename_with(., ~ reg_cols) %>% 
#     mutate(time = c(1:t)) %>% 
#     pivot_longer(cols = c(1:all_of(r)), values_to = "effect", names_to = "region") %>%
#     mutate(region = as.numeric(region), covar = vars[i])
# }
# 
# kappa_effects_sqrtY_v1 <- g1_sqrtY_effects_list[[1]] %>% left_join(., X_long_v1) %>% left_join(., full_reg_key)
# ggplot(kappa_effects_sqrtY_v1, aes(x=linear, y = effect, group = region)) +
#   geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE))
# 
# kappa_effects_sqrtY_v2 <- g1_sqrtY_effects_list[[2]] %>% left_join(., X_long_v2) %>% left_join(., full_reg_key)
# ggplot(kappa_effects_sqrtY_v2, aes(x=linear, y = effect, group = region)) +
#   geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE))
# 
# kappa_effects_sqrtY_v3 <- g1_sqrtY_effects_list[[3]] %>% left_join(., X_long_v3) %>% left_join(., full_reg_key)
# ggplot(kappa_effects_sqrtY_v3, aes(x=linear, y = effect, group = region)) +
#   geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE))
# 
# kappa_effects_sqrtY_v4 <- g1_sqrtY_effects_list[[4]] %>% left_join(., X_long_v4) %>% left_join(., full_reg_key)
# ggplot(kappa_effects_sqrtY_v4, aes(x=linear, y = effect, group = region)) +
#   geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE))
# 
# kappa_effects_sqrtY_v5 <- g1_sqrtY_effects_list[[5]] %>% left_join(., X_long_v5) %>% left_join(., full_reg_key)
# ggplot(kappa_effects_sqrtY_v5, aes(x=linear, y = effect, group = region)) +
#   geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE))
# 
# kappa_effects_sqrtY_v6 <- g1_sqrtY_effects_list[[6]] %>% left_join(., X_long_v6) %>% left_join(., full_reg_key)
# ggplot(kappa_effects_sqrtY_v6, aes(x=linear, y = effect, group = region)) +
#   geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE))
# 
# 
# g1_ogY_effects_list <- vector("list", 6)
# for(i in seq_along(vars)) {
#   g1_ogY_effects_df <- matrix(NA, t, r)
#   for(j in 1:r) {
#     g1_ogY_effects_df[, j] <- X[j, , X_cols[[i]]] %*% g1_ogY_kappa[X_cols[[i]], j]
#   }
#   g1_ogY_effects_list[[i]] <- g1_ogY_effects_df %>% 
#     as_tibble() %>% 
#     rename_with(., ~ reg_cols) %>% 
#     mutate(time = c(1:t)) %>% 
#     pivot_longer(cols = c(1:all_of(r)), values_to = "effect", names_to = "region") %>%
#     mutate(region = as.numeric(region), covar = vars[i])
# }
# 
# kappa_effects_ogY_v1 <- g1_ogY_effects_list[[1]] %>% left_join(., X_long_v1) %>% left_join(., full_reg_key)
# ggplot(kappa_effects_ogY_v1, aes(x=linear, y = effect, group = region)) +
#   geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE))
# 
# kappa_effects_ogY_v2 <- g1_ogY_effects_list[[2]] %>% left_join(., X_long_v2) %>% left_join(., full_reg_key)
# ggplot(kappa_effects_ogY_v2, aes(x=linear, y = effect, group = region)) +
#   geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE))
# 
# kappa_effects_ogY_v3 <- g1_ogY_effects_list[[3]] %>% left_join(., X_long_v3) %>% left_join(., full_reg_key)
# ggplot(kappa_effects_ogY_v3, aes(x=linear, y = effect, group = region)) +
#   geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE))
# 
# kappa_effects_ogY_v4 <- g1_ogY_effects_list[[4]] %>% left_join(., X_long_v4) %>% left_join(., full_reg_key)
# ggplot(kappa_effects_ogY_v4, aes(x=linear, y = effect, group = region)) +
#   geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE))
# 
# kappa_effects_ogY_v5 <- g1_ogY_effects_list[[5]] %>% left_join(., X_long_v5) %>% left_join(., full_reg_key)
# ggplot(kappa_effects_ogY_v5, aes(x=linear, y = effect, group = region)) +
#   geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE))
# 
# kappa_effects_ogY_v6 <- g1_ogY_effects_list[[6]] %>% left_join(., X_long_v6) %>% left_join(., full_reg_key)
# ggplot(kappa_effects_ogY_v6, aes(x=linear, y = effect, group = region)) +
#   geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE))
# 
