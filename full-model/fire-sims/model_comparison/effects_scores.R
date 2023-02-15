library(tidyverse)
library(cowplot)
library(tidyverse)
library(patchwork)
library(ggrepel)
library(rstan)
library(stringr)


# following code is for the counts -------
count_fits <- paste0("full-model/fire-sims/counts/", list.files(path = "full-model/fire-sims/counts/", pattern = "*.RDS", recursive = TRUE))
extraction <- function(file_path) {
  object <- readRDS(file_path)
  file <- basename(file_path)
  model <- str_remove(file, "_\\d{2}-Jan-2023_\\d{4}_\\d{4}.RDS")
  loglik_chr <- str_subset(object@model_pars, "log")
  betas_chr <- str_subset(object@model_pars, "beta")
  loglik <- extract(object, pars = loglik_chr)
  betas <- extract(object, pars = betas_chr)
  temp <- list(betas, loglik)
  assign(model, temp, parent.frame())
  rm(object)
}

for(i in 1:length(count_fits)) {
  extraction(count_fits[i])
}

count_names <- lapply(count_fits, function(x) str_remove(basename(x), "_\\d{2}-Jan-2023_\\d{4}_\\d{4}.RDS")) %>% unlist()
holdout_loglik_counts <- vector("list", length(count_names))
train_loglik_counts <- vector("list", length(count_names))
for(i in seq_along(count_names)) {
  count_loglik <- get(count_names[[i]])[[2]]
  holdout_loglik_counts[[i]] <- count_loglik$holdout_loglik %>%
    apply(., 1, c) %>%
    as_tibble() %>%
    pivot_longer(cols = everything(), names_to = "iter") %>%
    mutate(iter = as.numeric(gsub("V", "", iter))) %>%
    group_by(iter) %>%
    summarize(value = sum(value)) %>%
    mutate(model = count_names[i], train = FALSE)
  train_loglik_counts[[i]] <- count_loglik$train_loglik %>%
    apply(., 1, c) %>%
    as_tibble() %>%
    pivot_longer(cols = everything(), names_to = "iter") %>%
    mutate(iter = as.numeric(gsub("V", "", iter))) %>%
    group_by(iter) %>%
    summarize(value = sum(value)) %>%
    mutate(model = count_names[i], train = TRUE)
}
# 
holdout_loglik_c <- bind_rows(holdout_loglik_counts)
train_loglik_c <- bind_rows(train_loglik_counts)
# 
ll_full <- holdout_loglik_c %>%
  full_join(train_loglik_c) %>%
  mutate(train = ifelse(train == TRUE, 'train', 'test')) %>%
  pivot_wider(names_from = train, values_from = value)

ll_full %>% pivot_longer(cols = c("test", "train"), names_to = "dataset") %>%
  ggplot(aes(model, value, color = dataset)) + geom_boxplot() + theme_minimal()
ggsave("full-model/figures/model-comp/logscores_counts.png", dpi = 320, type = "cairo", bg = "white")

ll_full %>%
  group_by(model) %>%
  summarize(mean_test = mean(test),
            sd_test = sd(test)) %>%
  arrange(-mean_test)

ll_full %>% select(c(1:3)) %>% pivot_wider(names_from = model, values_from = test) %>% 
  mutate(across(.cols = c(2:7), ~ .x - `zinb_er_pi-ri`)) %>% pivot_longer(cols = c(2:7), names_to = "model") %>%
  group_by(model) %>%
  summarize(mean_diff = mean(value),
            sd_diff = sd(value)) %>%
  arrange(-mean_diff)

ll_full %>%
  group_by(model) %>%
  summarize(mean_train = mean(train),
            sd_train = sd(train)) %>%
  arrange(-mean_train)

ll_full %>% select(c(1,2,4)) %>% pivot_wider(names_from = model, values_from = train) %>% 
  mutate(across(.cols = c(2:7), ~ .x - `zip_pi-reg`)) %>% pivot_longer(cols = c(2:7), names_to = "model") %>%
  group_by(model) %>%
  summarize(mean_diff = mean(value),
            sd_diff = sd(value)) %>%
  arrange(-mean_diff)

saveRDS(ll_full, file = "full-model/figures/model-comp/ll_full_counts.RDS")

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
r <- 84
t <- 252
lambda_counts <- vector("list", length(count_names))
for(i in seq_along(count_names)) {
  count_beta <- get(count_names[[i]])[[1]]
  temp_df <- count_beta$beta_lambda %>%
    apply(., c(2,3), median)
  coef_df_list <- list()
  for(k in seq_along(vars)) {
    stored_df <- matrix(NA, t, r)
    for(j in 1:r) {
      stored_df[, j] <- X[j, , X_cols[[k]]] %*% temp_df[X_cols[[k]], j]
    }
    coef_df_list[[k]] <- stored_df %>% 
      as_tibble() %>% 
      rename_with(., ~ reg_cols) %>% 
      mutate(time = c(1:t)) %>% 
      pivot_longer(cols = c(1:all_of(r)), values_to = "effect", names_to = "region") %>%
      mutate(region = as.numeric(region), covar = vars[k], linear = c(X[,,X_cols[[k]][2]]))
  }
  lambda_counts[[i]] <- bind_rows(coef_df_list) %>% as_tibble() %>% mutate(model = count_names[i]) %>% left_join(., full_reg_key)
}

for(i in 1:length(count_names)) {
  p <- ggplot(lambda_counts[[i]], aes(x = linear, y = effect, group = region)) + 
    geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE), show.legend = FALSE) +
    facet_wrap(. ~ covar, scales = "free_x") + theme_minimal() + ggtitle(count_names[i])
  file_name <- paste0("full-model/figures/model-comp/", count_names[i], ".png")
  ggsave(file_name, p, dpi = 320, type = "cairo", bg = "white")
}



# following code is for the burns log scores and plots -------
burn_fits_g1 <- paste0("full-model/fire-sims/burns/g1/", 
                       list.files(path = "full-model/fire-sims/burns/g1/", pattern = "*.RDS", recursive = TRUE))
no_mix_g1 <- which(str_detect(burn_fits_g1, "nu-ri"))
burn_fits_g1 <- burn_fits_g1[-no_mix_g1]

extraction <- function(file_path) {
  object <- readRDS(file_path)
  file <- basename(file_path)
  model <- str_remove(file, "\\d{2}-Jan-2023_\\d{4}_\\d{4}.RDS")
  loglik_chr <- str_subset(object@model_pars, "log")
  betas_chr <- str_subset(object@model_pars, "beta")
  loglik <- extract(object, pars = loglik_chr)
  betas <- extract(object, pars = betas_chr)
  kappa_reg <- extract(object, pars = c("kappa_train", "kappa_hold")) # for use in twCRPS calculations
  sigma_reg <- extract(object, pars = c("sigma_train", "sigma_hold"))
  xi_reg <- extract(object, pars = c("xi_train", "xi_hold")) 
  temp <- list(betas, loglik, kappa_reg, sigma_reg, xi_reg)
  assign(model, temp, parent.frame())
  rm(object)
}

for(i in 1:length(burn_fits_g1)) {
  extraction(burn_fits_g1[i])
}

burn_names <- lapply(burn_fits_g1, function(x) str_remove(basename(x), "\\d{2}-Jan-2023_\\d{4}_\\d{4}.RDS")) %>% unlist()
holdout_loglik_burns <- vector("list", length(burn_names))
train_loglik_burns <- vector("list", length(burn_names))
for(i in seq_along(burn_names)) {
  burns_loglik <- get(burn_names[[i]])[[2]]
  holdout_loglik_burns[[i]] <- burns_loglik$holdout_loglik %>%
    apply(., 1, c) %>%
    as_tibble() %>%
    pivot_longer(cols = everything(), names_to = "iter") %>%
    mutate(iter = as.numeric(gsub("V", "", iter))) %>%
    group_by(iter) %>%
    summarize(value = sum(value)) %>%
    mutate(model = burn_names[i], train = FALSE)
  train_loglik_burns[[i]] <- burns_loglik$train_loglik %>%
    apply(., 1, c) %>%
    as_tibble() %>%
    pivot_longer(cols = everything(), names_to = "iter") %>%
    mutate(iter = as.numeric(gsub("V", "", iter))) %>%
    group_by(iter) %>%
    summarize(value = sum(value)) %>%
    mutate(model = burn_names[i], train = TRUE)
}
# 
holdout_loglik_b <- bind_rows(holdout_loglik_burns)
train_loglik_b <- bind_rows(train_loglik_burns)
# 
ll_full <- holdout_loglik_b %>%
  full_join(train_loglik_b) %>%
  mutate(train = ifelse(train == TRUE, 'train', 'test')) %>%
  pivot_wider(names_from = train, values_from = value)

ll_full %>% filter(test != Inf | train != Inf) %>% 
  pivot_longer(cols = c("test", "train"), names_to = "dataset") %>%
  ggplot(aes(model, value, color = dataset)) + geom_boxplot() + theme_minimal()
ggsave("full-model/figures/model-comp/logscores_burns.png", dpi = 320, type = "cairo", bg = "white")

ll_full %>%
  group_by(model) %>%
  summarize(mean_test = mean(test),
            sd_test = sd(test)) %>%
  arrange(-mean_test)

ll_full %>% select(c(1:3)) %>% pivot_wider(names_from = model, values_from = test) %>% 
  mutate(across(.cols = c(2:7), ~ .x - `zinb_er_pi-ri`)) %>% pivot_longer(cols = c(2:7), names_to = "model") %>%
  group_by(model) %>%
  summarize(mean_diff = mean(value),
            sd_diff = sd(value)) %>%
  arrange(-mean_diff)

ll_full %>%
  group_by(model) %>%
  summarize(mean_train = mean(train),
            sd_train = sd(train)) %>%
  arrange(-mean_train)

ll_full %>% select(c(1,2,4)) %>% pivot_wider(names_from = model, values_from = train) %>% 
  mutate(across(.cols = c(2:7), ~ .x - `zip_pi-reg`)) %>% pivot_longer(cols = c(2:7), names_to = "model") %>%
  group_by(model) %>%
  summarize(mean_diff = mean(value),
            sd_diff = sd(value)) %>%
  arrange(-mean_diff)

saveRDS(ll_full, file = "full-model/figures/model-comp/ll_full_burns.RDS")

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
r <- 84
t <- 252

nokappa <-  which(str_detect(burn_fits_g1, "kappa-ri"))
burn_names_kappa <- burn_names[-nokappa]
burn_names_nu <- burn_names
kappa_burns <- vector("list", length(burn_names_kappa))
nu_burns <- vector("list", length(burn_names_nu))
for(i in seq_along(burn_names_kappa)) {
  burns_beta <- get(burn_names_kappa[[i]])[[1]]
  temp_df <- burns_beta$beta_kappa %>%
    apply(., c(2,3), median)
  coef_df_list <- list()
  for(k in seq_along(vars)) {
    stored_df <- matrix(NA, t, r)
    for(j in 1:r) {
      stored_df[, j] <- X[j, , X_cols[[k]]] %*% temp_df[X_cols[[k]], j]
    }
    coef_df_list[[k]] <- stored_df %>% 
      as_tibble() %>% 
      rename_with(., ~ reg_cols) %>% 
      mutate(time = c(1:t)) %>% 
      pivot_longer(cols = c(1:all_of(r)), values_to = "effect", names_to = "region") %>%
      mutate(region = as.numeric(region), covar = vars[k], linear = c(X[,,X_cols[[k]][2]]))
  }
  kappa_burns[[i]] <- bind_rows(coef_df_list) %>% as_tibble() %>% mutate(model = burn_names_kappa[i]) %>% left_join(., full_reg_key)
}

for(i in 1:length(burn_names_kappa)) {
  p <- ggplot(kappa_burns[[i]], aes(x = linear, y = effect, group = region)) + 
    geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE), show.legend = FALSE) +
    facet_wrap(. ~ covar, scales = "free_x") + theme_minimal() + ggtitle(burn_names[i])
  file_name <- paste0("full-model/figures/model-comp/kappa_", burn_names_kappa[i], ".png")
  ggsave(file_name, p, dpi = 320, type = "cairo", bg = "white")
}

for(i in seq_along(burn_names_nu)) {
  burns_beta <- get(burn_names_nu[[i]])[[1]]
  temp_df <- burns_beta$beta_nu %>%
    apply(., c(2,3), median)
  coef_df_list <- list()
  for(k in seq_along(vars)) {
    stored_df <- matrix(NA, t, r)
    for(j in 1:r) {
      stored_df[, j] <- X[j, , X_cols[[k]]] %*% temp_df[X_cols[[k]], j]
    }
    coef_df_list[[k]] <- stored_df %>% 
      as_tibble() %>% 
      rename_with(., ~ reg_cols) %>% 
      mutate(time = c(1:t)) %>% 
      pivot_longer(cols = c(1:all_of(r)), values_to = "effect", names_to = "region") %>%
      mutate(region = as.numeric(region), covar = vars[k], linear = c(X[,,X_cols[[k]][2]]))
  }
  nu_burns[[i]] <- bind_rows(coef_df_list) %>% as_tibble() %>% mutate(model = burn_names_nu[i]) %>% left_join(., full_reg_key)
}

for(i in 1:length(burn_names_nu)) {
  p <- ggplot(nu_burns[[i]], aes(x = linear, y = effect, group = region)) + 
    geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE), show.legend = FALSE) +
    facet_wrap(. ~ covar, scales = "free_x") + theme_minimal() + ggtitle(burn_names[i])
  file_name <- paste0("full-model/figures/model-comp/nu_", burn_names_nu[i], ".png")
  ggsave(file_name, p, dpi = 320, type = "cairo", bg = "white")
}


## twCRPS calculations ---------
extraction(burn_fits_g1[3]) # g1_nu-reg_xi-regog
data_og <- readRDS("full-model/data/stan_data_og.RDS")
data_sqrt <- readRDS("full-model/data/stan_data_sqrt.RDS")
yhold <- data_og$y_hold_obs
ytrain <- data_og$y_train_obs

yholdsqrt <- data_sqrt$y_hold_obs
ytrain <- data_og$y_train_obs
kappa_vals <- `g1_nu-reg_xi-regog`[[3]]
sigma_vals <- `g1_nu-reg_xi-regog`[[4]]
xi_vals <- `g1_nu-reg_xi-regog`[[5]]
kappa_hold <- kappa_vals[[2]]
sigma_hold <- sigma_vals[[2]]
xi_hold <- xi_vals[[2]]

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
test <- g1_rng(400, sigma_hold[1,1], xi_hold[1,1], kappa_hold[1,1])

# NEED TO FIGURE OUT APPROPRIATE CHAIN SNAKING ----

twCRPS <- function(obs, n, sigma, xi, kappa) {
  int <- max(obs) - min(obs)
  yn <- length(obs)
  delta <- int/n
  i <- 1:n
  score_chn <- matrix(rep(NA, yn * 3), ncol = 3)
  all_scores <- rep(NA, yn)
  for(elem in 1:yn) {
    for(chn in 1:3) {
      for(iter in 1:1000) {
        forecast <- g1_rng(n, sigma[iter*chn, elem], xi[iter*chn, elem], kappa[iter*chn, elem])
        temp_ecdf <- ecdf(forecast)
        score_chn[elem, chn] <- delta * sum((temp_ecdf(forecast[i]) - (obs[elem] <= forecast[i]))^2 *
                                                pnorm(forecast[i], mean = 20, sd = 4))
      }
    }
  all_scores[elem] <- apply(score_chn, 1, mean) 
  }
  return(mean(all_scores))
}

twCRPS(yhold, 10000, sigma_hold, xi_hold, kappa_hold)


sum((emp(test) - (yhold[1] <= test))^2)/400
plot(ecdf(yhold), xlim = c(0,100), col = 'blue')
curve(pnorm(x, mean = 50, sd = 15), add = TRUE, col = 'red')
