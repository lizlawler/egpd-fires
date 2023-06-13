library(cmdstanr)
# set_cmdstan_path(path = "/projects/eslawler@colostate.edu/software/anaconda/envs/lawler/bin/cmdstan") # this is only relevant to Alpine
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(tidyverse)
library(stringr)
library(posterior)


# following code is for the counts -------
burn_fits <- paste0("full-model/fire-sims/burns/g1/", 
                    list.files(path = "full-model/fire-sims/burns/g1", 
                               pattern = "*.csv", recursive = TRUE))
gqidx <- which(grepl("gq_", burn_fits))
burn_fits <- burn_fits[-gqidx]
nfits <- length(burn_fits)/3
fit_groups <- vector(mode = "list", nfits)
for(i in 1:nfits) {
  fit_groups[[i]] <- burn_fits[(3*i-2):(3*i)]
}

burn_names <- lapply(fit_groups, function(x) str_remove(str_remove(basename(x[1]), "_\\d{2}\\w{3}2023_\\d{4}_\\d{1}.csv"), "cfcns_")) %>% unlist()
all_idx <- which(grepl("all-reg", burn_names))
sqrt_xi_idx <- which(grepl("sqrt_xi-ri", burn_names))
og_xi_idx <- which(grepl("og_xi-ri", burn_names))
kappa_xi_idx <- which(grepl("kappa-ri_xi-ri", burn_names))
nu_reg_idx <- c(all_idx, sqrt_xi_idx, og_xi_idx, kappa_xi_idx) %>% sort()
nu_reg_burns <- burn_names[nu_reg_idx]
sqrt_only <- grepl("sqrt", nu_reg_burns)
nu_reg_burns <- nu_reg_burns[sqrt_only]
ogdelta <- grepl("0.81", nu_reg_burns)
nu_reg_burns <- nu_reg_burns[ogdelta]

idx <- which(burn_names %in% nu_reg_burns)

extraction <- function(file_group) {
  object <- as_cmdstan_fit(file_group)
  file <- basename(file_group[1])
  model <- str_remove(str_remove(file, "_\\d{2}\\w{3}2023_\\d{4}_\\d{1}.csv"), "cfcns_")
  if (grepl("xi-ri", file)) {
    betas <- object$draws(variables = "beta")
    ri_matrix <- object$draws(variables = "ri_matrix")
    temp <- list(betas = betas, ri_matrix = ri_matrix)
  } else {
    betas <- object$draws(variables = "beta")
    temp <- list(betas = betas)
  }
  assign(model, temp, parent.frame())
  rm(object)
  gc()
}

for(i in idx) {
  extraction(fit_groups[[i]])
}

save(list=ls(pattern="g1_sqrt_xi-ri_0.81"), file = "full-model/fire-sims/model_comparison/nu_data/g1_sqrt_xi-ri_0.81.RData")
rm("g1_sqrt_xi-ri_0.81")
gc()
# extraction(fit_groups[[1]])
for(i in 1:20) {
  extraction(fit_groups[[i]])
}

burn_names <- lapply(fit_groups, function(x) str_remove(basename(x[1]), "_\\d{2}\\w{3}2023_\\d{4}_\\d{1}.csv")) %>% unlist()
rm("g1_sqrt_all-reg_0.81")
gc()
save.image()
## log score calculations ---------
# holdout_loglik_counts <- vector("list", nfits)
# train_loglik_counts <- vector("list", nfits)
# for(i in seq_along(burn_names)) {
#   count_loglik <- get(burn_names[i])[2:3]
#   holdout_loglik_counts[[i]] <- count_loglik[[2]] %>%
#     apply(., c(1,2), sum) %>% 
#     as_tibble() %>% 
#     rowid_to_column(var = "iter") %>%
#     mutate(mean_ax_cns = rowMeans(select(., !iter))) %>%
#     pivot_longer(cols = !iter, names_to = "chain") %>% 
#     mutate(chain = as.factor(chain)) %>%
#     mutate(model = burn_names[i], train = FALSE)
#   train_loglik_counts[[i]] <- count_loglik[[1]] %>%
#     apply(., c(1,2), sum) %>% 
#     as_tibble() %>% 
#     rowid_to_column(var = "iter") %>%
#     mutate(mean_ax_cns = rowMeans(select(., !iter))) %>%
#     pivot_longer(cols = !iter, names_to = "chain") %>% 
#     mutate(chain = as.factor(chain)) %>%
#     mutate(model = burn_names[i], train = TRUE)
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
## ---------

stan_data_og <- readRDS("full-model/data/stan_data_og.rds")
X_og <- stan_data_og$X_train
stan_data_sqrt <- readRDS("full-model/data/stan_data_sqrt.rds")
X_sqrt <- stan_data_sqrt$X_train
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
t <- stan_data_og$T_train

kappa_burns <- list()
nu_burns <- list()
xi_burns <- list()
egpd_param <- function(df, param_num) {
  return(
    df %>% select(-c("combo", "dataset", "delta")) %>%
      filter(param == param_num) %>% select(-c("param", "model_name")) %>%
      arrange(coef, region) %>% pivot_wider(names_from = region, values_from = value) %>%
      select(-coef) %>% as.matrix()
  ) 
}
covar_effect <- function(egpd_param_df, covar_term, linear_term) {
  return(
    egpd_param_df %>% as_tibble() %>% rename_with(., ~ reg_cols) %>% 
      mutate(time = c(1:t)) %>% 
      pivot_longer(cols = c(1:all_of(r)), values_to = "effect", names_to = "region") %>%
      mutate(region = as.numeric(region), covar = covar_term, linear = linear_term)
  )
}

for(i in seq_along(burn_beta_df_list)) {
  temp_df <- burn_beta_df_list[[i]] %>%
    mutate(dataset = 
             case_when(grepl("g1_og", model_name) ~ "og",
                       grepl("g1_sqrt", model_name) ~ "sqrt"),
           combo = 
             case_when(grepl("kappa-ri", model_name) ~ "kappa-ri",
                       grepl("sigma-ri", model_name) ~ "sigma-ri",
                       grepl("nu-ri", model_name) ~ "nu-ri",
                       grepl("xi-ri", model_name) ~ "xi-ri",
                       grepl("all-reg", model_name) ~ "all-reg",
                       TRUE ~ model_name),
           delta = 
             case_when(grepl("0.81", model_name) ~ "0.81",
                       grepl("0.9", model_name) ~ "0.9"),
           model_name = paste0(dataset, "_", combo, "_", delta),
           param = as.numeric(param),
           coef = as.numeric(coef),
           region = as.numeric(region))
  model_name <- unique(temp_df$model_name)
  if (unique(temp_df$dataset) == "og") {
    X <- X_og
  } else {
    X <- X_sqrt
  }
  if (unique(temp_df$combo) == "all-reg") {
    kappa_matrix <- egpd_param(temp_df, 1)
    nu_matrix <- egpd_param(temp_df, 2)
    xi_matrix <- egpd_param(temp_df, 3)
    coef_df_list_kappa <- list()
    coef_df_list_nu <- list()
    coef_df_list_xi <- list()
    for(k in seq_along(vars)) {
      stored_df_kappa <- matrix(NA, t, r)
      stored_df_nu <- matrix(NA, t, r)
      stored_df_xi <- matrix(NA, t, r)
      for(j in 1:r) {
        stored_df_kappa[, j] <- X[j, , X_cols[[k]]] %*% kappa_matrix[X_cols[[k]], j]
        stored_df_nu[, j] <- X[j, , X_cols[[k]]] %*% nu_matrix[X_cols[[k]], j]
        stored_df_xi[, j] <- X[j, , X_cols[[k]]] %*% xi_matrix[X_cols[[k]], j]
      }
      coef_df_list_kappa[[k]] <- covar_effect(stored_df_kappa, vars[k], c(X[,,X_cols[[k]][2]]))
      coef_df_list_nu[[k]] <- covar_effect(stored_df_nu, vars[k], c(X[,,X_cols[[k]][2]]))
      coef_df_list_xi[[k]] <- covar_effect(stored_df_xi, vars[k], c(X[,,X_cols[[k]][2]]))
    }
    kappa_burns[[i]] <- bind_rows(coef_df_list_kappa) %>% as_tibble() %>% mutate(model = model_name) %>% left_join(., full_reg_key)
    nu_burns[[i]] <- bind_rows(coef_df_list_nu) %>% as_tibble() %>% mutate(model = model_name) %>% left_join(., full_reg_key)
    xi_burns[[i]] <- bind_rows(coef_df_list_xi) %>% as_tibble() %>% mutate(model = model_name) %>% left_join(., full_reg_key)
  }  else if (unique(temp_df$combo) == "xi-ri") {
    kappa_matrix <- egpd_param(temp_df, 1)
    nu_matrix <- egpd_param(temp_df, 2)
    coef_df_list_kappa <- list()
    coef_df_list_nu <- list()
    for(k in seq_along(vars)) {
      stored_df_kappa <- matrix(NA, t, r)
      stored_df_nu <- matrix(NA, t, r)
      for(j in 1:r) {
        stored_df_kappa[, j] <- X[j, , X_cols[[k]]] %*% kappa_matrix[X_cols[[k]], j]
        stored_df_nu[, j] <- X[j, , X_cols[[k]]] %*% nu_matrix[X_cols[[k]], j]
      }
      coef_df_list_kappa[[k]] <- covar_effect(stored_df_kappa, vars[k], c(X[,,X_cols[[k]][2]]))
      coef_df_list_nu[[k]] <- covar_effect(stored_df_nu, vars[k], c(X[,,X_cols[[k]][2]]))
    }
    kappa_burns[[i]] <- bind_rows(coef_df_list_kappa) %>% as_tibble() %>% mutate(model = model_name) %>% left_join(., full_reg_key)
    nu_burns[[i]] <- bind_rows(coef_df_list_nu) %>% as_tibble() %>% mutate(model = model_name) %>% left_join(., full_reg_key)
  } else if (unique(temp_df$combo) == "kappa-ri") {
    nu_matrix <- egpd_param(temp_df, 1)
    coef_df_list_nu <- list()
    for(k in seq_along(vars)) {
      stored_df_nu <- matrix(NA, t, r)
      for(j in 1:r) {
        stored_df_nu[, j] <- X[j, , X_cols[[k]]] %*% nu_matrix[X_cols[[k]], j]
      }
      coef_df_list_nu[[k]] <- covar_effect(stored_df_nu, vars[k], c(X[,,X_cols[[k]][2]]))
    }
    nu_burns[[i]] <- bind_rows(coef_df_list_nu) %>% as_tibble() %>% mutate(model = model_name) %>% left_join(., full_reg_key)    
  } else {
    kappa_matrix <- egpd_param(temp_df, 1)
    coef_df_list_kappa <- list()
    for(k in seq_along(vars)) {
      stored_df_kappa <- matrix(NA, t, r)
      for(j in 1:r) {
        stored_df_kappa[, j] <- X[j, , X_cols[[k]]] %*% kappa_matrix[X_cols[[k]], j]
      }
      coef_df_list_kappa[[k]] <- covar_effect(stored_df_kappa, vars[k], c(X[,,X_cols[[k]][2]]))
    }
    kappa_burns[[i]] <- bind_rows(coef_df_list_kappa) %>% as_tibble() %>% mutate(model = model_name) %>% left_join(., full_reg_key)
  }
}


for(i in seq_along(burn_names)) {
  # count_beta <- get(burn_names[[i]])[[1]]
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
  lambda_counts[[i]] <- bind_rows(coef_df_list) %>% as_tibble() %>% mutate(model = burn_names[i]) %>% left_join(., full_reg_key)
}


ggplot(nu_burns[[3]], aes(x = linear, y = effect, group = region)) + 
  geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE), show.legend = FALSE) +
  facet_wrap(. ~ covar, scales = "free_x") + theme_minimal() + ggtitle(burn_names[3])
file_name <- paste0("full-model/figures/model-comp/", burn_names[i], ".png")


for(i in 1:length(kappa_burns)) {
  if (!is.null(kappa_burns[[i]])) {
    p <- ggplot(kappa_burns[[i]], aes(x = linear, y = effect, group = region)) + 
      geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE)) +
      facet_wrap(. ~ covar, scales = "free_x") + theme_minimal() + ggtitle(paste0("kappa_", unique(kappa_burns[[i]]$model)))
    file_name <- paste0("full-model/figures/model-comp/g1_kappa_", unique(kappa_burns[[i]]$model), ".png")
    ggsave(file_name, p, dpi = 320, type = "cairo", bg = "white")
  }
  else {
    print("element of list is empty")
  }
}

# for(i in 1:length(kappa_burns)) {
#   if (!is.null(kappa_burns[[i]])) {
#     print(paste0(i, "_", unique(kappa_burns[[i]]$model)))
#   }
#   else {
#     print("element of list is empty")
#   }
# }

for(i in 1:length(nu_burns)) {
  if (!is.null(nu_burns[[i]])) {
    p <- ggplot(nu_burns[[i]], aes(x = linear, y = effect, group = region)) + 
      geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE), show.legend = FALSE) +
      facet_wrap(. ~ covar, scales = "free_x") + theme_minimal() + ggtitle(paste0("nu_", unique(nu_burns[[i]]$model)))
    file_name <- paste0("full-model/figures/model-comp/g1_nu_", unique(nu_burns[[i]]$model), ".png")
    ggsave(file_name, p, dpi = 320, type = "cairo", bg = "white")
  }
  else {
    print("element of list is empty")
  }
}

for(i in 1:length(xi_burns)) {
  if (!is.null(xi_burns[[i]])) {
    p <- ggplot(xi_burns[[i]], aes(x = linear, y = effect, group = region)) + 
      geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE), show.legend = FALSE) +
      facet_wrap(. ~ covar, scales = "free_x") + theme_minimal() + ggtitle(paste0("xi_", unique(xi_burns[[i]]$model)))
    file_name <- paste0("full-model/figures/model-comp/g1_xi_", unique(xi_burns[[i]]$model), ".png")
    ggsave(file_name, p, dpi = 320, type = "cairo", bg = "white")
  }
  else {
    print("element of list is empty")
  }
}


poster_kappa <- kappa_burns[[17]] %>% filter(covar != "log_housing_density", covar != "pr", covar != "rmin", covar != "tmmx") %>%
  mutate(covar = case_when(
    covar == 'vs' ~ 'Wind speed',
    covar == 'prev_12mo_precip' ~ 'Precipitation: 12 month'),
    covar = factor(covar, levels = c("Wind speed", "Precipitation: 12 month")))

p <- ggplot(poster_kappa, aes(x=linear, y=effect, group = region)) + 
  geom_line(aes(linetype=NA_L1CODE, color = NA_L2CODE), show.legend = FALSE) + 
  facet_grid(. ~ covar) + xlab("Linear term") + ylab(expression("Partial effect on"~kappa))+
  theme_classic(base_size = 25)

ggsave("~/Desktop/research/posters-presentations/prelim_kappa.png", dpi = 700, type = "cairo", bg="white",
       height = 11, width = 17)



