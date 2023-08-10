library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(tidyverse)
library(stringr)
library(posterior)

# load in extracted gen quants done through Alpine
load("~/Desktop/research/egpd-fires/full-model/fire-sims/model_comparison/gq_newdata_counts.RData")

## log score calculations ---------
nfits <- length(count_names)
holdout_loglik_counts <- vector("list", nfits)
train_loglik_counts <- vector("list", nfits)
for(i in seq_along(count_names)) {
  model_string <- str_split(count_names[i], pattern = "_")[[1]]
  if(length(model_string) > 3) {
    model <- "zinb_er"
    params <- model_string[3]
    dataset <- model_string[4]
  } else {
    model <- model_string[1]
    params <- model_string[2]
    dataset <- model_string[3]
  }
  holdout_loglik_counts[[i]] <- get(count_names[i])[["holdout_loglik"]] %>%
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
  train_loglik_counts[[i]] <- get(count_names[i])[["train_loglik"]] %>%
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
}

holdout_loglik_c <- bind_rows(holdout_loglik_counts)
train_loglik_c <- bind_rows(train_loglik_counts)
# 
ll_full <- holdout_loglik_c %>%
  full_join(train_loglik_c) %>% mutate(full_name = paste(model, params, dataset, sep = "_"))

# ll_boxplot_train <- ll_full %>% 
#   ggplot(aes(full_name, loglik, color = train)) + geom_boxplot() + theme_minimal()
# ggsave("full-model/figures/model-comp/logscores_counts_train_31jul2023.png", plot = ll_boxplot_train,
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
  summarize(med_train_ll = median(loglik)) %>% arrange(-med_train_ll)
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
  summarize(med_train_ll = median(loglik)) %>% arrange(-med_train_ll)
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
  summarize(med_train_ll = median(loglik)) %>% arrange(-med_train_ll)
top_mod_test_erc <- as.character(test_ll_ranked_erc$full_name[1])
test_ll_comp_erc <- ll_full %>% filter(train == FALSE, dataset == "erc") %>% 
  select(c(draw, full_name, loglik)) %>% 
  pivot_wider(names_from = full_name, values_from = loglik, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_test_erc))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(med_diff = median(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% arrange(-med_diff)
test_ll_comp_erc

test_ll_ranked_zinber <- ll_full %>% filter(train == FALSE, model == "zinb_er") %>% 
  group_by(dataset, params, full_name) %>% 
  summarize(med_train_ll = median(loglik)) %>% arrange(-med_train_ll)
top_mod_test_zinber <- as.character(test_ll_ranked_zinber$full_name[1])
test_ll_comp_zinber <- ll_full %>% filter(train == FALSE, model == "zinb_er") %>% 
  select(c(draw, full_name, loglik)) %>% 
  pivot_wider(names_from = full_name, values_from = loglik, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_test_zinber))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(med_diff = median(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% arrange(-med_diff)
test_ll_comp_zinber

saveRDS(ll_full, file = "full-model/figures/model-comp/ll_scores_counts_31jul2023.RDS")

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
lambda_counts <- vector("list", length(count_names))

best_climate <- `zinb_er_pi-ri_climate`$betas %>%
  as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>% 
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  separate_wider_delim(cols = "name", delim = ",", names = c("param", "coef", "region"))
best_climate <- best_climate %>% 
  mutate(param = as.numeric(gsub("beta\\[", "", param)),
         coef = as.numeric(coef),
         region = as.numeric(gsub("\\]", "", region)))
lambda_climate <- best_climate %>% select(-param) %>% 
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

coef_df_list_climate <- list()
for(k in seq_along(vars)) {
  stored_df_climate <- matrix(NA, t, r)
  for(j in 1:r) {
    stored_df_climate[, j] <- X[j, , X_cols[[k]]] %*% lambda_climate[X_cols[[k]], j]
  }
  coef_df_list_climate[[k]] <- covar_effect(stored_df_climate, vars[k], c(X[,,X_cols[[k]][2]]))
}
lambda_climate_effects <- bind_rows(coef_df_list_climate) %>% as_tibble() %>% left_join(., full_reg_key)
p <- ggplot(lambda_climate_effects, aes(x = linear, y = effect, group = region)) + 
  geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE), show.legend = FALSE) +
  facet_wrap(. ~ covar, scales = "free_x") + theme_minimal() + ggtitle("climate")
file_name <- paste0("full-model/figures/model-comp/lambda_climate", ".png")
ggsave(file_name, p, dpi = 320, bg = "white")

stan_data_erc <- readRDS("full-model/data/stan_data_erc.RDS")
X <- stan_data_erc$X_train
vars <- c('log_housing_density', 'erc')
X_covar <- c()
X_cols <- vector("list", length(vars))
start <- 2
for(i in seq_along(vars)) {
  X_covar[i] <- paste0("X_", vars[i])
  X_cols[[i]] <- c(1, start:(start+5))
  start = start + 6
}


best_erc <- `zinb_er_pi-ri_erc`$betas %>%
  as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>% 
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  separate_wider_delim(cols = "name", delim = ",", names = c("param", "coef", "region"))
best_erc <- best_erc %>% 
  mutate(param = as.numeric(gsub("beta\\[", "", param)),
         coef = as.numeric(coef),
         region = as.numeric(gsub("\\]", "", region)))
lambda_erc <- best_erc %>% select(-param) %>% 
  group_by(region, coef) %>% summarize(med_val = median(value)) %>% ungroup() %>%
  pivot_wider(names_from = "region", values_from = "med_val") %>% select(-coef) %>% as.matrix()

coef_df_list_erc <- list()
for(k in seq_along(vars)) {
  stored_df_erc <- matrix(NA, t, r)
  for(j in 1:r) {
    stored_df_erc[, j] <- X[j, , X_cols[[k]]] %*% lambda_erc[X_cols[[k]], j]
  }
  coef_df_list_erc[[k]] <- covar_effect(stored_df_erc, vars[k], c(X[,,X_cols[[k]][2]]))
}
lambda_erc_effects <- bind_rows(coef_df_list_erc) %>% as_tibble() %>% left_join(., full_reg_key)
p <- ggplot(lambda_erc_effects, aes(x = linear, y = effect, group = region)) + 
  geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE), show.legend = FALSE) +
  facet_wrap(. ~ covar, scales = "free_x") + theme_minimal() + ggtitle("erc")
file_name <- paste0("full-model/figures/model-comp/lambda_erc", ".png")
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

best_fwi <- `zinb_er_pi-ri_fwi`$betas %>%
  as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>% 
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  separate_wider_delim(cols = "name", delim = ",", names = c("param", "coef", "region"))
best_fwi <- best_fwi %>% 
  mutate(param = as.numeric(gsub("beta\\[", "", param)),
         coef = as.numeric(coef),
         region = as.numeric(gsub("\\]", "", region)))
lambda_fwi <- best_fwi %>% select(-param) %>% 
  group_by(region, coef) %>% summarize(med_val = median(value)) %>% ungroup() %>%
  pivot_wider(names_from = "region", values_from = "med_val") %>% select(-coef) %>% as.matrix()

coef_df_list_fwi <- list()
for(k in seq_along(vars)) {
  stored_df_fwi <- matrix(NA, t, r)
  for(j in 1:r) {
    stored_df_fwi[, j] <- X[j, , X_cols[[k]]] %*% lambda_fwi[X_cols[[k]], j]
  }
  coef_df_list_fwi[[k]] <- covar_effect(stored_df_fwi, vars[k], c(X[,,X_cols[[k]][2]]))
}
lambda_fwi_effects <- bind_rows(coef_df_list_fwi) %>% as_tibble() %>% left_join(., full_reg_key)
p <- ggplot(lambda_fwi_effects, aes(x = linear, y = effect, group = region)) + 
  geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE), show.legend = FALSE) +
  facet_wrap(. ~ covar, scales = "free_x") + theme_minimal() + ggtitle("fwi")
file_name <- paste0("full-model/figures/model-comp/lambda_fwi", ".png")
ggsave(file_name, p, dpi = 320, bg = "white")


for(i in seq_along(count_names)) {
  # count_beta <- get(count_names[[i]])[[1]]
  temp_df <- count_beta[[i]] %>%
    apply(., c(2,3), median) %>% 
    as_tibble() %>% pivot_longer(cols=everything(), names_to = "param_coef_reg", values_to = "value") %>%
    group_by(param_coef_reg) %>% summarize(value = mean(value)) %>%
    mutate(param_coef_reg = str_extract(param_coef_reg, "\\d{1},\\d{1,},\\d{1,}")) %>%
    separate(., param_coef_reg, into=c("param", "coef", "region"), ",") %>%
    mutate(param = case_when(
      grepl("1", param) ~ "lambda",
      grepl("2", param) ~ "pi",
      grepl("3", param) ~ "delta",
      TRUE ~ param),
      param = as.factor(param),
      coef = as.numeric(coef),
      region = as.numeric(region)) %>%
    filter(param == "lambda") %>% select(-param) %>% arrange(coef, region) %>%
    pivot_wider(names_from = region, values_from = value) %>% select(-coef) %>% as.matrix()
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