library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(tidyverse)
library(stringr)
library(posterior)

# load in extracted gen quants done through Alpine
score_files <- paste0("full-model/fire-sims/model_comparison/extracted_values/",
                      list.files("full-model/fire-sims/model_comparison/extracted_values/", pattern = "scores"))
score_files <- score_files[!grepl("all-reg", score_files)]
model_names <- str_remove(str_remove(basename(score_files), "_scores.RDS"), "joint_")
new_models <- score_files[grepl("theta-time", score_files)]
new_names <- str_remove(str_remove(basename(new_models), "_scores.RDS"), "joint_")
for(i in seq_along(new_names)) {
  assign(new_names[i], readRDS(new_models[i]))
}

erc_fwi_files <- paste0("full-model/fire-sims/burns/g1/csv-fits/",
                        list.files("full-model/fire-sims/burns/g1/csv-fits/", pattern = "erc_fwi"))
erc_fwi_files <- erc_fwi_files[!grepl("all-reg", erc_fwi_files)]
nfits <- length(erc_fwi_files)/3
fit_groups <- vector(mode = "list", nfits)
for(i in 1:nfits) {
  fit_groups[[i]] <- erc_fwi_files[(3*i-2):(3*i)]
}

extraction <- function(file_group, burn_name) {
  object <- as_cmdstan_fit(file_group)
  
  reg <- object$draws(variables = "reg")
  ri_init <- object$draws(variables = "ri_init")
  pi_prob <- object$draws(variables = "pi_prob")
  lambda <- object$draws(variables = "lambda")
  delta <- object$draws(variables = "delta")
  theta <- object$draws(variables = "theta")
  gamma <- object$draws(variables = "gamma")
  phi <- object$draws(variables = "phi")
  temp <- list(reg, ri_init, pi_prob, lambda, delta, theta, gamma, phi)
  names(temp) <- c("reg", "ri_init", "pi_prob", "lambda", "delta", "theta", "gamma", "phi")
  assign(burn_name, temp, parent.frame())
  rm(object)
  gc()
}

extraction <- function(file_group, burn_name) {
  object <- as_cmdstan_fit(file_group)
  betas <- object$draws(variables = "beta")
  train_loglik <- object$draws(variables = "train_loglik")
  holdout_loglik <- object$draws(variables = "holdout_loglik")
  train_twcrps <- object$draws(variables = "train_twcrps")
  holdout_twcrps <- object$draws(variables = "holdout_twcrps")
  temp <- list(betas, train_loglik, holdout_loglik, train_twcrps, holdout_twcrps)
  names(temp) <- c("betas", "train_loglik", "holdout_loglik", "train_twcrps", "holdout_twcrps")
  assign(burn_name, temp, parent.frame())
  rm(object)
  gc()
}

burn_names <- lapply(fit_groups, function(x) str_remove(str_remove(basename(x[1]), "_\\d{2}\\w{3}2023_\\d{4}_\\d{1}.csv"), "g1_")) %>% unlist()

for(i in seq_along(burn_names)) {
  extraction(fit_groups[[i]], burn_names[i])
}


## log score calculations ---------
nfits <- length(burn_names)
# holdout_loglik_count <- vector("list", nfits)
# train_loglik_count <- vector("list", nfits)
holdout_loglik_burns <- vector("list", nfits)
train_loglik_burns <- vector("list", nfits)
train_twcrps <- vector("list", nfits)
holdout_twcrps <- vector("list", nfits)
for(i in seq_along(model_names)) {
  model_string <- str_split(model_names[i], pattern = "_")[[1]]
  model <- model_string[1]
  rand_effect <- model_string[2]
  link <- model_string[3]
  holdout_loglik_burns[[i]] <- get(model_names[i])[["holdout_loglik_burn"]] %>%
    as_draws_df() %>%
    select(-c(".iteration", ".chain")) %>% 
    pivot_longer(cols = !".draw") %>%
    rename(draw = ".draw") %>%
    group_by(draw) %>% 
    summarize(loglik = sum(value)) %>%
    mutate(model = model,
           rand_effect = rand_effect,
           link = link,
           train = FALSE)
  # holdout_loglik_count[[i]] <- get(model_names[i])[["holdout_loglik_count"]] %>%
  #   as_draws_df() %>%
  #   select(-c(".iteration", ".chain")) %>% 
  #   pivot_longer(cols = !".draw") %>%
  #   rename(draw = ".draw") %>%
  #   group_by(draw) %>% 
  #   summarize(loglik = sum(value)) %>%
  #   mutate(model = model,
  #          rand_effect = rand_effect,
  #          link = link,
  #          train = FALSE)
  train_loglik_burns[[i]] <- get(model_names[i])[["train_loglik_burn"]] %>%
    as_draws_df() %>%
    select(-c(".iteration", ".chain")) %>% 
    pivot_longer(cols = !".draw") %>%
    rename(draw = ".draw") %>%
    group_by(draw) %>% 
    summarize(loglik = sum(value)) %>%
    mutate(model = model,
           rand_effect = rand_effect,
           link = link,
           train = TRUE)
  # train_loglik_count[[i]] <- get(model_names[i])[["train_loglik_count"]] %>%
  #   as_draws_df() %>%
  #   select(-c(".iteration", ".chain")) %>% 
  #   pivot_longer(cols = !".draw") %>%
  #   rename(draw = ".draw") %>%
  #   group_by(draw) %>% 
  #   summarize(loglik = sum(value)) %>%
  #   mutate(model = model,
  #          rand_effect = rand_effect,
  #          link = link,
  #          train = TRUE)
  train_twcrps[[i]] <- get(model_names[i])[["train_twcrps"]] %>%
    as_draws_df() %>%
    select(-c(".iteration", ".chain")) %>% 
    pivot_longer(cols = !".draw") %>%
    rename(draw = ".draw") %>%
    group_by(draw) %>% 
    summarize(twcrps = mean(value)) %>%
    mutate(model = model,
           rand_effect = rand_effect,
           link = link,
           train = TRUE)
  holdout_twcrps[[i]] <- get(model_names[i])[["holdout_twcrps"]] %>%
    as_draws_df() %>%
    select(-c(".iteration", ".chain")) %>% 
    pivot_longer(cols = !".draw") %>%
    rename(draw = ".draw") %>%
    group_by(draw) %>% 
    summarize(twcrps = mean(value)) %>%
    mutate(model = model,
           rand_effect = rand_effect,
           link = link,
           train = FALSE)
}

# holdout_loglik_count <- bind_rows(holdout_loglik_count)
# train_loglik_count <- bind_rows(train_loglik_count)
holdout_loglik_burns <- bind_rows(holdout_loglik_burns)
train_loglik_burns <- bind_rows(train_loglik_burns)
holdout_twcrps <- bind_rows(holdout_twcrps)
train_twcrps <- bind_rows(train_twcrps)

ll_full_count <- holdout_loglik_count %>%
  full_join(train_loglik_count) %>% mutate(full_name = paste(model, rand_effect, link, sep = "_"))

# ll_boxplot_train <- ll_full %>% 
#   ggplot(aes(full_name, loglik, color = train)) + geom_boxplot() + theme_minimal()
# ggsave("full-model/figures/model-comp/logscores_burns_train_31jul2023.png", plot = ll_boxplot_train,
#        dpi = 320, bg = "white")

train_ll_count_ranked <- ll_full_count %>% filter(train == TRUE) %>% 
  group_by(model, rand_effect, link, full_name) %>% 
  summarize(med_train_ll = median(loglik)) %>% arrange(-med_train_ll)
top_mod_count_train <- as.character(train_ll_count_ranked$full_name[1])

train_ll_count_comp <- ll_full_count %>% filter(train == TRUE) %>% 
  select(c(draw, full_name, loglik)) %>% 
  pivot_wider(names_from = full_name, values_from = loglik, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_count_train))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(med_diff = median(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% arrange(-med_diff)

test_ll_count_ranked <- ll_full_count %>% filter(train == FALSE) %>% 
  group_by(model, rand_effect, link, full_name) %>% 
  summarize(med_train_ll = median(loglik)) %>% arrange(-med_train_ll)
top_mod_count_test <- as.character(test_ll_count_ranked$full_name[1])
test_ll_count_comp <- ll_full_count %>% filter(train == FALSE) %>% 
  select(c(draw, full_name, loglik)) %>% 
  pivot_wider(names_from = full_name, values_from = loglik, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_count_test))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(med_diff = median(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% arrange(-med_diff)
test_ll_count_comp


ll_full_burns <- holdout_loglik_burns %>%
  full_join(train_loglik_burns) %>% mutate(full_name = paste(model, rand_effect, link, sep = "_"))

train_ll_burns_ranked <- ll_full_burns %>% filter(train == TRUE) %>% 
  group_by(model, rand_effect, link, full_name) %>% 
  summarize(med_train_ll = median(loglik)) %>% arrange(-med_train_ll)
top_mod_burns_train <- as.character(train_ll_burns_ranked$full_name[1])

train_ll_burns_comp <- ll_full_burns %>% filter(train == TRUE) %>% 
  select(c(draw, full_name, loglik)) %>% 
  pivot_wider(names_from = full_name, values_from = loglik, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_burns_train))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(med_diff = median(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% arrange(-med_diff)

test_ll_burns_ranked <- ll_full_burns %>% filter(train == FALSE) %>% 
  group_by(model, rand_effect, link, full_name) %>% 
  summarize(med_train_ll = median(loglik)) %>% arrange(-med_train_ll)
top_mod_burns_test <- as.character(test_ll_burns_ranked$full_name[1])
test_ll_burns_comp <- ll_full_burns %>% filter(train == FALSE) %>% 
  select(c(draw, full_name, loglik)) %>% 
  pivot_wider(names_from = full_name, values_from = loglik, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_burns_test))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(med_diff = median(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% arrange(-med_diff)
test_ll_burns_comp


## twCRPS calculations ---------
twcrps_full <- holdout_twcrps %>%
  full_join(train_twcrps) %>% mutate(full_name = paste(model, rand_effect, link, sep = "_"))

train_twcrps_ranked <- twcrps_full %>% filter(train == TRUE) %>% 
  group_by(model, rand_effect, link, full_name) %>% 
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
  group_by(model, rand_effect, link, full_name) %>% 
  summarize(mean_twcrps = mean(twcrps, na.rm = TRUE)) %>% arrange(mean_twcrps)
top_mod_test_twcrps <- as.character(test_twcrps_ranked$full_name[1])
test_twcrps_comp <- twcrps_full %>% filter(train == FALSE) %>% 
  select(c(draw, full_name, twcrps)) %>% 
  pivot_wider(names_from = full_name, values_from = twcrps, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_test_twcrps))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(mean_diff = mean(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% arrange(mean_diff)
test_twcrps_comp


## effects plots ----------
stan_data_joint <- readRDS("full-model/data/stan_data_joint.RDS")
X_count <- stan_data_joint$X_train_count
X_burn <- stan_data_joint$X_train_burn
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
t <- stan_data_joint$T_train

best_joint_lambda <- `joint_sigma-ri_theta-ri_gamma-ri_betas`$beta_count %>%
  as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>% 
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  separate_wider_delim(cols = "name", delim = ",", names = c("coef", "region"))
best_joint_lambda <- best_joint_lambda %>% 
  mutate(coef = as.numeric(gsub("beta_count\\[", "", coef)),
         region = as.numeric(gsub("\\]", "", region)))
best_joint_kappa <- `joint_sigma-ri_theta-ri_gamma-ri_betas`$beta_burn %>%
  as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>% 
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  separate_wider_delim(cols = "name", delim = ",", names = c("coef", "region"))
best_joint_kappa <- best_joint_kappa %>% 
  mutate(coef = as.numeric(gsub("beta_burn\\[", "", coef)),
         region = as.numeric(gsub("\\]", "", region)))

best_joint_lambda <- best_joint_lambda %>% group_by(region, coef) %>% summarize(med_val = median(value)) %>% ungroup() %>%
  pivot_wider(names_from = "region", values_from = "med_val") %>% select(-coef) %>% as.matrix()
best_joint_kappa <- best_joint_kappa %>% group_by(region, coef) %>% summarize(med_val = median(value)) %>% ungroup() %>%
  pivot_wider(names_from = "region", values_from = "med_val") %>% select(-coef) %>% as.matrix()

covar_effect <- function(egpd_param_df, covar_term, linear_term) {
  return(
    egpd_param_df %>% as_tibble() %>% rename_with(., ~ as.character(reg_cols)) %>%
      mutate(time = c(1:t)) %>%
      pivot_longer(cols = c(1:all_of(r)), values_to = "effect", names_to = "region") %>%
      mutate(region = as.numeric(region), covar = covar_term, linear = linear_term)
  )
}

coef_df_list_lambda <- list()
for(k in seq_along(vars)) {
  stored_df_lambda <- matrix(NA, t, r)
  for(j in 1:r) {
    stored_df_lambda[, j] <- X_count[j, , X_cols[[k]]] %*% best_joint_lambda[X_cols[[k]], j]
  }
  coef_df_list_lambda[[k]] <- covar_effect(stored_df_lambda, vars[k], c(X_count[,,X_cols[[k]][2]]))
}
lambda_effects <- bind_rows(coef_df_list_lambda) %>% as_tibble() %>% left_join(., full_reg_key)
p <- ggplot(lambda_effects, aes(x = linear, y = effect, group = region)) + 
  geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE), show.legend = FALSE) +
  facet_wrap(. ~ covar, scales = "free_x") + theme_minimal() + ggtitle("lambda_effects")
file_name <- paste0("full-model/figures/joint/effects/lambda_effects", ".png")
ggsave(file_name, p, dpi = 320, bg = "white")


vars <- c('log_housing_density', 'erc')
X_covar <- c()
X_cols <- vector("list", length(vars))
start <- 2
for(i in seq_along(vars)) {
  X_covar[i] <- paste0("X_", vars[i])
  X_cols[[i]] <- c(1, start:(start+5))
  start = start + 6
}

coef_df_list_kappa <- list()
for(k in seq_along(vars)) {
  stored_df_kappa <- matrix(NA, t, r)
  for(j in 1:r) {
    stored_df_kappa[, j] <- X_burn[j, , X_cols[[k]]] %*% best_joint_kappa[X_cols[[k]], j]
  }
  coef_df_list_kappa[[k]] <- covar_effect(stored_df_kappa, vars[k], c(X_burn[,,X_cols[[k]][2]]))
}
kappa_effects <- bind_rows(coef_df_list_kappa) %>% as_tibble() %>% left_join(., full_reg_key)
p <- ggplot(kappa_effects, aes(x = linear, y = effect, group = region)) + 
  geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE), show.legend = FALSE) +
  facet_wrap(. ~ covar, scales = "free_x") + theme_minimal() + ggtitle("kappa_effects")
file_name <- paste0("full-model/figures/joint/effects/kappa_effects", ".png")
ggsave(file_name, p, dpi = 320, bg = "white")

