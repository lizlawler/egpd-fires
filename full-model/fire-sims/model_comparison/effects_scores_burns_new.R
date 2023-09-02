library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(tidyverse)
library(stringr)
library(posterior)

# load in extracted gen quants done through Alpine
load("~/Desktop/research/egpd-fires/full-model/fire-sims/model_comparison/gq_newdata_burns.RData")

## log score calculations ---------
nfits <- length(burn_names)
holdout_loglik_burns <- vector("list", nfits)
train_loglik_burns <- vector("list", nfits)
train_twcrps <- vector("list", nfits)
holdout_twcrps <- vector("list", nfits)
for(i in seq_along(burn_names)) {
  model_string <- str_split(burn_names[i], pattern = "_")[[1]]
  if(length(model_string) > 3) {
    model <- model_string[1]
    params <- paste(model_string[2], model_string[3], sep = "_")
    dataset <- model_string[4]
  } else {
    model <- model_string[1]
    params <- model_string[2]
    dataset <- model_string[3]
  }
  holdout_loglik_burns[[i]] <- get(burn_names[i])[["holdout_loglik"]] %>%
    as_draws_df() %>%
    select(-c(".iteration", ".chain")) %>% 
    pivot_longer(cols = !".draw") %>%
    rename(draw = ".draw") %>%
    group_by(draw) %>% 
    summarize(loglik = sum(value)) %>%
    mutate(model = model,
           dataset = dataset,
           params = params,
           train = FALSE)
  train_loglik_burns[[i]] <- get(burn_names[i])[["train_loglik"]] %>%
    as_draws_df() %>%
    select(-c(".iteration", ".chain")) %>% 
    pivot_longer(cols = !".draw") %>%
    rename(draw = ".draw") %>%
    group_by(draw) %>% 
    summarize(loglik = sum(value)) %>%
    mutate(model = model,
           dataset = dataset,
           params = params,
           train = TRUE)
  train_twcrps[[i]] <- get(burn_names[i])[["train_twcrps"]] %>%
    as_draws_df() %>%
    select(-c(".iteration", ".chain")) %>% 
    pivot_longer(cols = !".draw") %>%
    rename(draw = ".draw") %>%
    group_by(draw) %>% 
    summarize(twcrps = mean(value)) %>%
    mutate(model = model,
           dataset = dataset,
           params = params,
           train = TRUE)
  holdout_twcrps[[i]] <- get(burn_names[i])[["holdout_twcrps"]] %>%
    as_draws_df() %>%
    select(-c(".iteration", ".chain")) %>% 
    pivot_longer(cols = !".draw") %>%
    rename(draw = ".draw") %>%
    group_by(draw) %>% 
    summarize(twcrps = mean(value)) %>%
    mutate(model = model,
           dataset = dataset,
           params = params,
           train = FALSE)
}

holdout_loglik_c <- bind_rows(holdout_loglik_burns)
train_loglik_c <- bind_rows(train_loglik_burns)
holdout_twcrps_c <- bind_rows(holdout_twcrps)
train_twcrps_c <- bind_rows(train_twcrps)

ll_full <- holdout_loglik_c %>%
  full_join(train_loglik_c) %>% mutate(full_name = paste(model, params, dataset, sep = "_")) %>%
  filter(full_name != "lognorm_all-reg_climate")

# ll_boxplot_train <- ll_full %>% 
#   ggplot(aes(full_name, loglik, color = train)) + geom_boxplot() + theme_minimal()
# ggsave("full-model/figures/model-comp/logscores_burns_train_31jul2023.png", plot = ll_boxplot_train,
#        dpi = 320, bg = "white")

train_ll_ranked <- ll_full %>% filter(train == TRUE) %>% 
  group_by(model, params, dataset, full_name) %>% 
  summarize(med_train_ll = median(loglik)) %>% arrange(-med_train_ll)
top_mod_train <- as.character(train_ll_ranked$full_name[1])

train_ll_comp <- ll_full %>% filter(train == TRUE) %>% 
  select(c(draw, full_name, loglik)) %>% 
  pivot_wider(names_from = full_name, values_from = loglik, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_train))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(med_diff = median(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% arrange(-med_diff)

test_ll_ranked <- ll_full %>% filter(train == FALSE) %>% 
  group_by(model, params, dataset, full_name) %>% 
  summarize(med_train_ll = median(loglik)) %>% arrange(-med_train_ll)
top_mod_test <- as.character(test_ll_ranked$full_name[1])
test_ll_comp <- ll_full %>% filter(train == FALSE) %>% 
  select(c(draw, full_name, loglik)) %>% 
  pivot_wider(names_from = full_name, values_from = loglik, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_test))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(med_diff = median(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% arrange(-med_diff)
test_ll_comp

test_ll_ranked_climate <- ll_full %>% filter(train == FALSE, dataset == "climate") %>% 
  group_by(model, params, full_name) %>% 
  summarize(med_ll = median(loglik)) %>% arrange(-med_ll)
top_mod_test_climate <- as.character(test_ll_ranked_climate$full_name[1])
test_ll_comp_climate <- ll_full %>% filter(train == FALSE, dataset == "climate") %>% 
  select(c(draw, full_name, loglik)) %>% 
  pivot_wider(names_from = full_name, values_from = loglik, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_test_climate))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(med_diff = median(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% arrange(-med_diff)
test_ll_comp_climate

test_ll_ranked_fwi <- ll_full %>% filter(train == FALSE, dataset == "fwi") %>% 
  group_by(model, params, full_name) %>% 
  summarize(med_ll = median(loglik)) %>% arrange(-med_ll)
top_mod_test_fwi <- as.character(test_ll_ranked_fwi$full_name[1])
test_ll_comp_fwi <- ll_full %>% filter(train == FALSE, dataset == "fwi") %>% 
  select(c(draw, full_name, loglik)) %>% 
  pivot_wider(names_from = full_name, values_from = loglik, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_test_fwi))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(med_diff = median(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% arrange(-med_diff)
test_ll_comp_fwi

test_ll_ranked_erc <- ll_full %>% filter(train == FALSE, dataset == "erc") %>% 
  group_by(model, params, full_name) %>% 
  summarize(med_ll = median(loglik)) %>% arrange(-med_ll)
top_mod_test_erc <- as.character(test_ll_ranked_erc$full_name[1])
test_ll_comp_erc <- ll_full %>% filter(train == FALSE, dataset == "erc") %>% 
  select(c(draw, full_name, loglik)) %>% 
  pivot_wider(names_from = full_name, values_from = loglik, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_test_erc))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(med_diff = median(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% arrange(-med_diff)
test_ll_comp_erc

test_ll_ranked_g1 <- ll_full %>% filter(train == FALSE, model == "g1") %>% 
  group_by(dataset, params, full_name) %>% 
  summarize(med_ll = median(loglik)) %>% arrange(-med_ll)
top_mod_test_g1 <- as.character(test_ll_ranked_g1$full_name[1])
test_ll_comp_g1 <- ll_full %>% filter(train == FALSE, model == "g1") %>% 
  select(c(draw, full_name, loglik)) %>% 
  pivot_wider(names_from = full_name, values_from = loglik, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_test_g1))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(med_diff = median(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% arrange(-med_diff)
test_ll_comp_g1

## twCRPS calculations ---------
twcrps_full <- holdout_twcrps_c %>%
  full_join(train_twcrps_c) %>% mutate(full_name = paste(model, params, dataset, sep = "_")) %>%
  filter(full_name != "lognorm_all-reg_climate")

train_twcrps_ranked <- twcrps_full %>% filter(train == TRUE) %>% 
  group_by(model, params, dataset, full_name) %>% 
  summarize(mean_twcrps = mean(twcrps, na.rm = TRUE)) %>% arrange(mean_twcrps)
top_mod_train_twcrps <- as.character(train_twcrps_ranked$full_name[1])

train_twcrps_comp <- twcrps_full %>% filter(train == TRUE) %>% 
  select(c(draw, full_name, twcrps)) %>% 
  pivot_wider(names_from = full_name, values_from = twcrps, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_train_twcrps))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(mean_diff = mean(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% arrange(mean_diff)

test_twcrps_ranked <- twcrps_full %>% filter(train == FALSE) %>% 
  group_by(model, params, dataset, full_name) %>% 
  summarize(mean_twcrps = mean(twcrps, na.rm = TRUE)) %>% arrange(mean_twcrps)
top_mod_test_twcrps <- as.character(test_twcrps_ranked$full_name[1])
test_twcrps_comp <- twcrps_full %>% filter(train == FALSE) %>% 
  select(c(draw, full_name, twcrps)) %>% 
  pivot_wider(names_from = full_name, values_from = twcrps, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_test_twcrps))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(mean_diff = mean(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% arrange(mean_diff)

test_twcrps_ranked_climate <- twcrps_full %>% filter(train == FALSE, dataset == "climate") %>% 
  group_by(model, params, full_name) %>% 
  summarize(mean_twcrps = mean(twcrps, na.rm = TRUE)) %>% arrange(mean_twcrps)
top_mod_test_twcrps_climate <- as.character(test_twcrps_ranked_climate$full_name[1])
test_twcrps_comp_climate <- twcrps_full %>% filter(train == FALSE, dataset == "climate") %>% 
  select(c(draw, full_name, twcrps)) %>% 
  pivot_wider(names_from = full_name, values_from = twcrps, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_test_twcrps_climate))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(mean_diff = mean(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% arrange(mean_diff)
test_twcrps_comp_climate

test_twcrps_ranked_erc <- twcrps_full %>% filter(train == FALSE, dataset == "erc") %>% 
  group_by(model, params, full_name) %>% 
  summarize(mean_twcrps = mean(twcrps, na.rm = TRUE)) %>% arrange(mean_twcrps)
top_mod_test_twcrps_erc <- as.character(test_twcrps_ranked_erc$full_name[1])
test_twcrps_comp_erc <- twcrps_full %>% filter(train == FALSE, dataset == "erc") %>% 
  select(c(draw, full_name, twcrps)) %>% 
  pivot_wider(names_from = full_name, values_from = twcrps, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_test_twcrps_erc))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(mean_diff = mean(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% arrange(mean_diff)
test_twcrps_comp_erc

test_twcrps_ranked_fwi <- twcrps_full %>% filter(train == FALSE, dataset == "fwi") %>% 
  group_by(model, params, full_name) %>% 
  summarize(mean_twcrps = mean(twcrps, na.rm = TRUE)) %>% arrange(mean_twcrps)
top_mod_test_twcrps_fwi <- as.character(test_twcrps_ranked_fwi$full_name[1])
test_twcrps_comp_fwi <- twcrps_full %>% filter(train == FALSE, dataset == "fwi") %>% 
  select(c(draw, full_name, twcrps)) %>% 
  pivot_wider(names_from = full_name, values_from = twcrps, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_test_twcrps_fwi))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(mean_diff = mean(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% arrange(mean_diff)
test_twcrps_comp_fwi

test_twcrps_ranked_g1 <- twcrps_full %>% filter(train == FALSE, model == "g1") %>% 
  group_by(dataset, params, full_name) %>% 
  summarize(mean_twcrps = mean(twcrps, na.rm = TRUE)) %>% arrange(mean_twcrps)
top_mod_test_twcrps_g1 <- as.character(test_twcrps_ranked_g1$full_name[1])
test_twcrps_comp_g1 <- twcrps_full %>% filter(train == FALSE, model == "g1") %>% 
  select(c(draw, full_name, twcrps)) %>% 
  pivot_wider(names_from = full_name, values_from = twcrps, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_test_twcrps_g1))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(mean_diff = mean(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% arrange(mean_diff)
test_twcrps_comp_g1


## effects plots ----------
files <- paste0("full-model/fire-sims/burns/g1/csv-fits/",
                      list.files("full-model/fire-sims/burns/g1/csv-fits/", pattern = "erc_fwi"))
files <- files[grepl("sigma-ri", files)]
best_erc_fwi <- as_cmdstan_fit(files)
betas <- best_erc_fwi$draws(variables = "beta")

stan_data_climate <- readRDS("full-model/data/stan_data_climate.RDS")
X <- stan_data_climate$X_train
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

region_key <- readRDS(file = "./full-model/data/processed/region_key.rds")
full_reg_key <- as_tibble(region_key) %>% 
  mutate(region = c(1:84),
         NA_L2CODE = as.factor(NA_L2CODE),
         NA_L1CODE = as.factor(NA_L1CODE),
         NA_L3CODE = as.factor(NA_L3CODE))
reg_cols <- full_reg_key$region
r <- 84
t <- stan_data_climate$T_train

best_climate <- `g1_xi-ri_climate`$betas %>%
  as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>% 
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  separate_wider_delim(cols = "name", delim = ",", names = c("param", "coef", "region"))
best_climate <- best_climate %>% 
  mutate(param = as.numeric(gsub("beta\\[", "", param)),
         coef = as.numeric(coef),
         region = as.numeric(gsub("\\]", "", region)))
kappa_climate <- best_climate %>% filter(param == 1) %>% select(-param) %>% 
  group_by(region, coef) %>% summarize(med_val = median(value)) %>% ungroup() %>%
  pivot_wider(names_from = "region", values_from = "med_val") %>% select(-coef) %>% as.matrix()
sigma_climate <- best_climate %>% filter(param == 2) %>% select(-param) %>% 
  group_by(region, coef) %>% summarize(med_val = median(value)) %>% ungroup() %>%
  pivot_wider(names_from = "region", values_from = "med_val") %>% select(-coef) %>% as.matrix()

covar_effect <- function(egpd_param_df, covar_term, linear_term) {
  return(
    egpd_param_df %>% as_tibble() %>% rename_with(., ~ as.character(reg_cols)) %>%
      mutate(time = c(1:t)) %>%
      pivot_longer(cols = c(1:all_of(r)), values_to = "effect", names_to = "region") %>%
      mutate(region = as.numeric(region), covar = covar_term, linear = linear_term)
  )
}

coef_df_list_kappa_climate <- list()
coef_df_list_sigma_climate <- list()
for(k in seq_along(vars)) {
  stored_df_kappa_climate <- matrix(NA, t, r)
  stored_df_sigma_climate <- matrix(NA, t, r)
  for(j in 1:r) {
    stored_df_kappa_climate[, j] <- X[j, , X_cols[[k]]] %*% kappa_climate[X_cols[[k]], j]
    stored_df_sigma_climate[, j] <- X[j, , X_cols[[k]]] %*% sigma_climate[X_cols[[k]], j]
  }
  coef_df_list_kappa_climate[[k]] <- covar_effect(stored_df_kappa_climate, vars[k], c(X[,,X_cols[[k]][2]]))
  coef_df_list_sigma_climate[[k]] <- covar_effect(stored_df_sigma_climate, vars[k], c(X[,,X_cols[[k]][2]]))
}
kappa_climate_effects <- bind_rows(coef_df_list_kappa_climate) %>% as_tibble() %>% left_join(., full_reg_key)
sigma_climate_effects <- bind_rows(coef_df_list_sigma_climate) %>% as_tibble() %>% left_join(., full_reg_key)
p <- ggplot(kappa_climate_effects, aes(x = linear, y = effect, group = region)) + 
  geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE), show.legend = FALSE) +
  facet_wrap(. ~ covar, scales = "free_x") + theme_minimal() + ggtitle("kappa_climate")
file_name <- paste0("full-model/figures/model-comp/g1_kappa_climate", ".png")
ggsave(file_name, p, dpi = 320, bg = "white")
p <- ggplot(sigma_climate_effects, aes(x = linear, y = effect, group = region)) + 
  geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE), show.legend = FALSE) +
  facet_wrap(. ~ covar, scales = "free_x") + theme_minimal() + ggtitle("sigma_climate")
file_name <- paste0("full-model/figures/model-comp/g1_sigma_climate", ".png")
ggsave(file_name, p, dpi = 320, bg = "white")

stan_data_erc_fwi <- readRDS("full-model/data/stan_data_erc_fwi.RDS")
X <- stan_data_erc_fwi$X_train
vars <- c('log_housing_density', 'erc', 'fwi')
X_covar <- c()
X_cols <- vector("list", length(vars))
start <- 2
for(i in seq_along(vars)) {
  X_covar[i] <- paste0("X_", vars[i])
  X_cols[[i]] <- c(1, start:(start+5))
  start = start + 6
}

best_fit <- betas %>%
  as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>% 
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  separate_wider_delim(cols = "name", delim = ",", names = c("param", "coef", "region"))
best_fit <- best_fit %>% 
  mutate(param = as.numeric(gsub("beta\\[", "", param)),
         coef = as.numeric(coef),
         region = as.numeric(gsub("\\]", "", region)))
kappas <- best_fit %>% filter(param == 1) %>% select(-param) %>% 
  group_by(region, coef) %>% summarize(med_val = median(value)) %>% ungroup() %>%
  pivot_wider(names_from = "region", values_from = "med_val") %>% select(-coef) %>% as.matrix()
# sigma_erc <- best_erc %>% filter(param == 2) %>% select(-param) %>% 
#   group_by(region, coef) %>% summarize(med_val = median(value)) %>% ungroup() %>%
#   pivot_wider(names_from = "region", values_from = "med_val") %>% select(-coef) %>% as.matrix()

coef_df_list_kappa <- list()
# coef_df_list_sigma_erc <- list()
for(k in seq_along(vars)) {
  stored_df_kappa <- matrix(NA, t, r)
  # stored_df_sigma <- matrix(NA, t, r)
  for(j in 1:r) {
    stored_df_kappa[, j] <- X[j, , X_cols[[k]]] %*% kappas[X_cols[[k]], j]
    # stored_df_sigma[, j] <- X[j, , X_cols[[k]]] %*% sigma[X_cols[[k]], j]
  }
  coef_df_list_kappa[[k]] <- covar_effect(stored_df_kappa, vars[k], c(X[,,X_cols[[k]][2]]))
  # coef_df_list_sigma[[k]] <- covar_effect(stored_df_sigma, vars[k], c(X[,,X_cols[[k]][2]]))
}
kappa_effects <- bind_rows(coef_df_list_kappa) %>% as_tibble() %>% left_join(., full_reg_key)
# sigma_effects <- bind_rows(coef_df_list_sigma) %>% as_tibble() %>% left_join(., full_reg_key)
p <- ggplot(kappa_effects, aes(x = linear, y = effect, group = region)) + 
  geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE), show.legend = FALSE) +
  facet_wrap(. ~ covar, scales = "free_x") + theme_minimal() + ggtitle("kappa_erc-fwi")
file_name <- paste0("full-model/figures/model-comp/g1_kappa_erc-fwi", ".png")
ggsave(file_name, p, dpi = 320, bg = "white", width = 15)
p <- ggplot(sigma_erc_effects, aes(x = linear, y = effect, group = region)) + 
  geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE), show.legend = FALSE) +
  facet_wrap(. ~ covar, scales = "free_x") + theme_minimal() + ggtitle("sigma_erc")
file_name <- paste0("full-model/figures/model-comp/g1_sigma_erc", ".png")
ggsave(file_name, p, dpi = 320, bg = "white")

stan_data_fwi <- readRDS("full-model/data/stan_data_fwi.RDS")
X <- stan_data_fwi$X_train
vars <- c('log_housing_density', 'fwi')
X_covar <- c()
X_cols <- vector("list", length(vars))
start <- 2
for(i in seq_along(vars)) {
  X_covar[i] <- paste0("X_", vars[i])
  X_cols[[i]] <- c(1, start:(start+5))
  start = start + 6
}

best_fwi <- `g1_xi-ri_fwi`$betas %>%
  as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>% 
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  separate_wider_delim(cols = "name", delim = ",", names = c("param", "coef", "region"))
best_fwi <- best_fwi %>% 
  mutate(param = as.numeric(gsub("beta\\[", "", param)),
         coef = as.numeric(coef),
         region = as.numeric(gsub("\\]", "", region)))
kappa_fwi <- best_fwi %>% filter(param == 1) %>% select(-param) %>% 
  group_by(region, coef) %>% summarize(med_val = median(value)) %>% ungroup() %>%
  pivot_wider(names_from = "region", values_from = "med_val") %>% select(-coef) %>% as.matrix()
sigma_fwi <- best_fwi %>% filter(param == 2) %>% select(-param) %>% 
  group_by(region, coef) %>% summarize(med_val = median(value)) %>% ungroup() %>%
  pivot_wider(names_from = "region", values_from = "med_val") %>% select(-coef) %>% as.matrix()

coef_df_list_kappa_fwi <- list()
coef_df_list_sigma_fwi <- list()
for(k in seq_along(vars)) {
  stored_df_kappa_fwi <- matrix(NA, t, r)
  stored_df_sigma_fwi <- matrix(NA, t, r)
  for(j in 1:r) {
    stored_df_kappa_fwi[, j] <- X[j, , X_cols[[k]]] %*% kappa_fwi[X_cols[[k]], j]
    stored_df_sigma_fwi[, j] <- X[j, , X_cols[[k]]] %*% sigma_fwi[X_cols[[k]], j]
  }
  coef_df_list_kappa_fwi[[k]] <- covar_effect(stored_df_kappa_fwi, vars[k], c(X[,,X_cols[[k]][2]]))
  coef_df_list_sigma_fwi[[k]] <- covar_effect(stored_df_sigma_fwi, vars[k], c(X[,,X_cols[[k]][2]]))
}
kappa_fwi_effects <- bind_rows(coef_df_list_kappa_fwi) %>% as_tibble() %>% left_join(., full_reg_key)
sigma_fwi_effects <- bind_rows(coef_df_list_sigma_fwi) %>% as_tibble() %>% left_join(., full_reg_key)
p <- ggplot(kappa_fwi_effects, aes(x = linear, y = effect, group = region)) + 
  geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE), show.legend = FALSE) +
  facet_wrap(. ~ covar, scales = "free_x") + theme_minimal() + ggtitle("kappa_fwi")
file_name <- paste0("full-model/figures/model-comp/g1_kappa_fwi", ".png")
ggsave(file_name, p, dpi = 320, bg = "white")
p <- ggplot(sigma_fwi_effects, aes(x = linear, y = effect, group = region)) + 
  geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE), show.legend = FALSE) +
  facet_wrap(. ~ covar, scales = "free_x") + theme_minimal() + ggtitle("sigma_fwi")
file_name <- paste0("full-model/figures/model-comp/g1_sigma_fwi", ".png")
ggsave(file_name, p, dpi = 320, bg = "white")

