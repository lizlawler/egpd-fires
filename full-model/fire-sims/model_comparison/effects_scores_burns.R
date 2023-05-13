library(cmdstanr)
# set_cmdstan_path(path = "/projects/eslawler@colostate.edu/software/anaconda/envs/lawler/bin/cmdstan") # this is only relevant to Alpine
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(tidyverse)
# library(cowplot)
# library(patchwork)
# library(ggrepel)
library(stringr)
library(posterior)


# following code is for the counts -------
burn_fits <- paste0("full-model/fire-sims/burns/g1/", list.files(path = "full-model/fire-sims/burns/g1", pattern = "*.csv", recursive = TRUE))
nfits <- length(burn_fits)/3
fit_groups <- vector(mode = "list", nfits)
for(i in 1:nfits) {
  fit_groups[[i]] <- burn_fits[(3*i-2):(3*i)]
}

extraction <- function(file_group) {
  object <- as_cmdstan_fit(file_group)
  file <- basename(file_group[1])
  model <- str_remove(file, "_\\d{2}\\w{3}2023_\\d{4}_\\d{1}.csv")
  # train_loglik <- object$draws(variables = "train_loglik")
  # holdout_loglik <- object$draws(variables = "holdout_loglik")
  betas <- object$draws(variables = "beta")
  temp <- list(betas)
  assign(model, temp, parent.frame())
  rm(object)
  gc()
}

# extraction(fit_groups[[1]])
for(i in 1:3) {
  extraction(fit_groups[[i]])
}

save.image()
## log score calculations ---------
# burn_names <- lapply(fit_groups, function(x) str_remove(basename(x[1]), "_\\d{2}\\w{3}2023_\\d{4}_\\d{1}.csv")) %>% unlist()
# holdout_loglik_counts <- vector("list", nfits)
# train_loglik_counts <- vector("list", nfits)
# for(i in seq_along(count_names)) {
#   count_loglik <- get(count_names[i])[2:3]
#   holdout_loglik_counts[[i]] <- count_loglik[[2]] %>%
#     apply(., c(1,2), sum) %>% 
#     as_tibble() %>% 
#     rowid_to_column(var = "iter") %>%
#     mutate(mean_ax_cns = rowMeans(select(., !iter))) %>%
#     pivot_longer(cols = !iter, names_to = "chain") %>% 
#     mutate(chain = as.factor(chain)) %>%
#     mutate(model = count_names[i], train = FALSE)
#   train_loglik_counts[[i]] <- count_loglik[[1]] %>%
#     apply(., c(1,2), sum) %>% 
#     as_tibble() %>% 
#     rowid_to_column(var = "iter") %>%
#     mutate(mean_ax_cns = rowMeans(select(., !iter))) %>%
#     pivot_longer(cols = !iter, names_to = "chain") %>% 
#     mutate(chain = as.factor(chain)) %>%
#     mutate(model = count_names[i], train = TRUE)
# }
# 
# holdout_loglik_c <- bind_rows(holdout_loglik_counts)
# train_loglik_c <- bind_rows(train_loglik_counts)
# # 
# ll_full <- holdout_loglik_c %>%
#   full_join(train_loglik_c) %>%
#   mutate(train = ifelse(train == TRUE, 'train', 'test')) %>%
#   pivot_wider(names_from = train, values_from = value)
# 
# ll_full_0.81 <- ll_full %>% filter(grepl("0.81", model)) %>% 
#   mutate(model = gsub("_0.81", "", model),
#          model = gsub("_og", "", model))
# ll_full_0.90 <- ll_full %>% filter(!grepl("0.81", model)) %>% 
#   mutate(model = gsub("_0.9", "", model),
#          model = gsub("_og", "", model))
# 
# ll_boxplot_0.81 <- ll_full_0.81 %>% pivot_longer(cols = c("test", "train"), names_to = "dataset") %>% 
#   filter(chain == 'mean_ax_cns') %>% 
#   ggplot(aes(model, value, color = dataset)) + geom_boxplot() + theme_minimal()
# ggsave("full-model/figures/model-comp/logscores_counts_0.81_12may2023.png", plot = ll_boxplot_0.81,
#        dpi = 320, type = "cairo", bg = "white")
# 
# ll_boxplot_0.90 <- ll_full_0.90 %>% pivot_longer(cols = c("test", "train"), names_to = "dataset") %>% 
#   filter(chain == 'mean_ax_cns') %>% 
#   ggplot(aes(model, value, color = dataset)) + geom_boxplot() + theme_minimal()
# ggsave("full-model/figures/model-comp/logscores_counts_0.90_12may2023.png", plot = ll_boxplot_0.90,
#        dpi = 320, type = "cairo", bg = "white")
# 
# ll_ranked_test <- ll_full %>% filter(chain == 'mean_ax_cns') %>%
#   group_by(model) %>%
#   summarize(mean_test = mean(test),
#             sd_test = sd(test)) %>%
#   arrange(-mean_test)
# top_mod_test <- ll_ranked_test$model[1]
# 
# ll_comp_test <- ll_full %>% filter(chain == 'mean_ax_cns') %>% 
#   select(c(1, 3:4)) %>% pivot_wider(names_from = model, values_from = test) %>% 
#   mutate(across(.cols = -1, ~ .x - get(top_mod_test))) %>% 
#   pivot_longer(cols = -1, names_to = "model") %>%
#   group_by(model) %>%
#   summarize(mean_diff = mean(value),
#             sd_diff = sd(value)) %>%
#   arrange(-mean_diff)
# ll_comp_test
# 
# ll_ranked_test_0.81 <- ll_full_0.81 %>% filter(chain == 'mean_ax_cns') %>%
#   group_by(model) %>%
#   summarize(mean_test = mean(test),
#             sd_test = sd(test)) %>%
#   arrange(-mean_test)
# top_mod_test_0.81 <- ll_ranked_test_0.81$model[1]
# 
# ll_comp_test_0.81 <- ll_full_0.81 %>% filter(chain == 'mean_ax_cns') %>% 
#   select(c(1, 3:4)) %>% pivot_wider(names_from = model, values_from = test) %>% 
#   mutate(across(.cols = -1, ~ .x - get(top_mod_test_0.81))) %>% 
#   pivot_longer(cols = -1, names_to = "model") %>%
#   group_by(model) %>%
#   summarize(mean_diff = mean(value),
#             sd_diff = sd(value)) %>%
#   arrange(-mean_diff)
# ll_comp_test_0.81
# 
# ll_ranked_test_0.90 <- ll_full_0.90 %>% filter(chain == 'mean_ax_cns') %>%
#   group_by(model) %>%
#   summarize(mean_test = mean(test),
#             sd_test = sd(test)) %>%
#   arrange(-mean_test)
# top_mod_test_0.90 <- ll_ranked_test_0.90$model[1]
# 
# ll_comp_test_0.90 <- ll_full_0.90 %>% filter(chain == 'mean_ax_cns') %>% 
#   select(c(1, 3:4)) %>% pivot_wider(names_from = model, values_from = test) %>% 
#   mutate(across(.cols = -1, ~ .x - get(top_mod_test_0.90))) %>% 
#   pivot_longer(cols = -1, names_to = "model") %>%
#   group_by(model) %>%
#   summarize(mean_diff = mean(value),
#             sd_diff = sd(value)) %>%
#   arrange(-mean_diff)
# ll_comp_test_0.90
# 
# ll_ranked_train <- ll_full %>% filter(chain == 'mean_ax_cns') %>%
#   group_by(model) %>%
#   summarize(mean_train = mean(train),
#             sd_train = sd(train)) %>%
#   arrange(-mean_train)
# top_mod_train <- ll_ranked_train$model[1]
# 
# ll_comp_train <- ll_full %>% filter(chain == 'mean_ax_cns') %>% 
#   select(c(1,3,5)) %>% pivot_wider(names_from = model, values_from = train) %>% 
#   mutate(across(.cols = -1, ~ .x - get(top_mod_train))) %>% 
#   pivot_longer(cols = -1, names_to = "model") %>%
#   group_by(model) %>%
#   summarize(mean_diff = mean(value),
#             sd_diff = sd(value)) %>%
#   arrange(-mean_diff)
# ll_comp_train
# 
# ll_ranked_train_0.81 <- ll_full_0.81 %>% filter(chain == 'mean_ax_cns') %>%
#   group_by(model) %>%
#   summarize(mean_train = mean(train),
#             sd_train = sd(train)) %>%
#   arrange(-mean_train)
# top_mod_train_0.81 <- ll_ranked_train_0.81$model[1]
# 
# ll_comp_train_0.81 <- ll_full_0.81 %>% filter(chain == 'mean_ax_cns') %>% 
#   select(c(1,3,5)) %>% pivot_wider(names_from = model, values_from = train) %>% 
#   mutate(across(.cols = -1, ~ .x - get(top_mod_train_0.81))) %>% 
#   pivot_longer(cols = -1, names_to = "model") %>%
#   group_by(model) %>%
#   summarize(mean_diff = mean(value),
#             sd_diff = sd(value)) %>%
#   arrange(-mean_diff)
# ll_comp_train_0.81
# 
# ll_ranked_train_0.90 <- ll_full_0.90 %>% filter(chain == 'mean_ax_cns') %>%
#   group_by(model) %>%
#   summarize(mean_train = mean(train),
#             sd_train = sd(train)) %>%
#   arrange(-mean_train)
# top_mod_train_0.90 <- ll_ranked_train_0.90$model[1]
# 
# ll_comp_train_0.90 <- ll_full_0.90 %>% filter(chain == 'mean_ax_cns') %>% 
#   select(c(1,3,5)) %>% pivot_wider(names_from = model, values_from = train) %>% 
#   mutate(across(.cols = -1, ~ .x - get(top_mod_train_0.90))) %>% 
#   pivot_longer(cols = -1, names_to = "model") %>%
#   group_by(model) %>%
#   summarize(mean_diff = mean(value),
#             sd_diff = sd(value)) %>%
#   arrange(-mean_diff)
# ll_comp_train_0.90
# 
# saveRDS(ll_full, file = "full-model/figures/model-comp/ll_full_counts_12may2023.RDS")


stan_data <- readRDS("full-model/data/stan_data_og.rds")
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
t <- stan_data$T_train
lambda_counts <- vector("list", length(count_names))

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