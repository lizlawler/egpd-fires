library(cmdstanr)
# set_cmdstan_path(path = "/projects/eslawler@colostate.edu/software/anaconda/envs/lawler/bin/cmdstan") # this is only relevant to Alpine
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(tidyverse)
library(stringr)
library(posterior)

# following code is for extracting from the actual model fit -------
burn_fits <- paste0("full-model/fire-sims/burns/g1/csv-fits/",
                    list.files(path = "full-model/fire-sims/burns/g1/csv-fits",
                               pattern = "sigma-reg", recursive = TRUE))
allregidx <- which(grepl("all-reg", burn_fits))
burn_fits <- burn_fits[-allregidx]
nfits <- length(burn_fits)/3
fit_groups <- vector(mode = "list", nfits)
for(i in 1:nfits) {
  fit_groups[[i]] <- burn_fits[(3*i-2):(3*i)]
}

burn_names <- lapply(fit_groups, function(x) str_remove(str_remove(basename(x[1]), "_\\d{2}\\w{3}2023_\\d{4}_\\d{1}.csv"), "cfcns_")) %>% unlist()
# 
extraction <- function(file_group, burn_name) {
  object <- as_cmdstan_fit(file_group)
  betas <- object$draws(variables = "beta")
  ri_matrix <- object$draws(variables = "ri_matrix")
  temp <- list(betas, ri_matrix)
  names(temp) <- c("betas", "ri_matrix")
  assign(burn_name, temp, parent.frame())
  rm(object)
  gc()
}

for(i in seq_along(burn_names)) {
  extraction(fit_groups[[i]], burn_names[i])
}
# 
save(list=c(ls(pattern="g1"), ls(pattern = "burn_names")), file = "full-model/figures/g1/sigma_reg_draws.RData")
# rm("g1_sqrt_xi-ri_0.81")
# gc()
# # extraction(fit_groups[[1]])
# for(i in 1:20) {
#   extraction(fit_groups[[i]])
# }
# 
# burn_names <- lapply(fit_groups, function(x) str_remove(basename(x[1]), "_\\d{2}\\w{3}2023_\\d{4}_\\d{1}.csv")) %>% unlist()
# rm("g1_sqrt_all-reg_0.81")
# gc()
# save.image()
# 

# read in best model - G1 with regression on sigma and kappa
files <- paste0("full-model/fire-sims/burns/g1/csv-fits/", 
                list.files(path = "full-model/fire-sims/burns/g1/csv-fits/",
                           pattern = "og_xi-ri_sigma-reg", recursive = TRUE))
g1_og_xi_ri <- read_cmdstan_csv(files, variables = "ri_init")

stan_data_og <- readRDS("full-model/data/stan_data_og.rds")
X <- stan_data_og$X_train
# stan_data_sqrt <- readRDS("full-model/data/stan_data_sqrt.rds")
# X_sqrt <- stan_data_sqrt$X_train
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
# 
kappa_burns <- list()
nu_burns <- list()
xi_burns <- list()
# kappa_idx <- which(grepl("beta\\[1,", variables(`g1_sqrt_all-reg_0.81`[[1]])))
# nu_idx <- which(grepl("beta\\[2,", variables(`g1_sqrt_all-reg_0.81`[[1]])))
# xi_idx <- which(grepl("beta\\[3,", variables(`g1_sqrt_all-reg_0.81`[[1]])))
# betas_kappa <- `g1_sqrt_all-reg_0.81`[[1]][,,kappa_idx]
# betas_nu <- `g1_sqrt_all-reg_0.81`[[1]][,,nu_idx]
# betas_xi <- `g1_sqrt_all-reg_0.81`[[1]][,,xi_idx]

idx <- which(grepl("beta", variables(g1_og_xi_ri$post_warmup_draws)))
all_betas <- g1_og_xi_ri$post_warmup_draws[,,idx] %>% as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>% 
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>% 
  separate_wider_delim(cols = "name", delim = ",", names = c("param", "coef", "region"))
all_betas <- all_betas %>% 
  mutate(param = as.numeric(gsub("beta\\[", "", param)),
         coef = as.numeric(coef),
         region = as.numeric(gsub("\\]", "", region)))
kappa_betas <- all_betas %>% filter(param == 1) %>% select(-param)
sigma_betas <- all_betas %>% filter(param == 2) %>% select(-param)
kappa_medians <- kappa_betas %>% group_by(region, coef) %>% summarize(med_val = median(value)) %>% ungroup() %>%
  pivot_wider(names_from = "region", values_from = "med_val") %>% select(-coef) %>% as.matrix()

sigma_medians <- sigma_betas %>% group_by(region, coef) %>% summarize(med_val = median(value)) %>% ungroup() %>%
  pivot_wider(names_from = "region", values_from = "med_val") %>% select(-coef) %>% as.matrix()

covar_effect <- function(egpd_param_df, covar_term, linear_term) {
  return(
    egpd_param_df %>% as_tibble() %>% rename_with(., ~ as.character(reg_cols)) %>%
      mutate(time = c(1:t)) %>%
      pivot_longer(cols = c(1:all_of(r)), values_to = "effect", names_to = "region") %>%
      mutate(region = as.numeric(region), covar = covar_term, linear = linear_term)
  )
}

coef_df_list_kappa <- list()
coef_df_list_sigma <- list()
for(k in seq_along(vars)) {
  stored_df_kappa <- matrix(NA, t, r)
  stored_df_sigma <- matrix(NA, t, r)
  for(j in 1:r) {
    stored_df_kappa[, j] <- X[j, , X_cols[[k]]] %*% kappa_medians[X_cols[[k]], j]
    stored_df_sigma[, j] <- X[j, , X_cols[[k]]] %*% sigma_medians[X_cols[[k]], j]
  }
  coef_df_list_kappa[[k]] <- covar_effect(stored_df_kappa, vars[k], c(X[,,X_cols[[k]][2]]))
  coef_df_list_sigma[[k]] <- covar_effect(stored_df_sigma, vars[k], c(X[,,X_cols[[k]][2]]))
}

kappa_burns <- bind_rows(coef_df_list_kappa) %>% as_tibble() %>% left_join(., full_reg_key)
sigma_burns <- bind_rows(coef_df_list_sigma) %>% as_tibble() %>% left_join(., full_reg_key)

p <- ggplot(sigma_burns, aes(x = linear, y = effect, group = region)) + 
        geom_line(aes(color = NA_L2CODE)) +
        facet_wrap(. ~ covar, scales = "free_x") + theme_minimal() + ggtitle("sigma_effects")
file_name <- "full-model/figures/g1/effects/sigma_effects_bestmod_nolevel1.png"
ggsave(file_name, p, dpi = 320, type = "cairo", bg = "white")

p <- ggplot(kappa_burns, aes(x = linear, y = effect, group = region)) + 
  geom_line(aes(color = NA_L2CODE)) +
  facet_wrap(. ~ covar, scales = "free_x") + theme_minimal() + ggtitle("kappa_effects")
file_name <- "full-model/figures/g1/effects/kappa_effects_bestmod_nolevel1.png"
ggsave(file_name, p, dpi = 320, type = "cairo", bg = "white")

xi_idx <- which(grepl("ri_init", variables(g1_og_xi_ri$post_warmup_draws)))
xi_vals <- g1_og_xi_ri$post_warmup_draws[,,xi_idx] %>% as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>% 
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  mutate(region = readr::parse_number(name)) %>%
  select(-name) %>%
  mutate(exp_val = exp(value))

xi_vals_med <- xi_vals %>% group_by(region) %>% summarize(med_val = median(exp_val)) %>% ungroup()


# kappa_df_list <- lapply(split(kappa_df, kappa_df$iter), function(x) {select(x, -c("iter", "covar")) %>% as.matrix()})
# 
# iter_kappa_effects_list <- lapply(kappa_df_list, function(x) {
#   coef_df_list_kappa <- list()
#   for(k in seq_along(vars)) {
#     stored_df_kappa <- matrix(NA, t, r)
#     for(j in 1:r) {
#       stored_df_kappa[,j] <- X_sqrt[j,,X_cols[[k]]] %*% x[X_cols[[k]], j]
#     }
#     coef_df_list_kappa[[k]] <- covar_effect(stored_df_kappa, vars[k], c(X_sqrt[,,X_cols[[k]][2]]))
#   }
#   return(coef_df_list_kappa)
# })
# 
# 
# egpd_param <- function(df, param_num) {
#   return(
#     df %>% select(-c("combo", "dataset", "delta")) %>%
#       filter(param == param_num) %>% select(-c("param", "model_name")) %>%
#       arrange(coef, region) %>% pivot_wider(names_from = region, values_from = value) %>%
#       select(-coef) %>% as.matrix()
#   ) 
# }
# covar_effect <- function(egpd_param_df, covar_term, linear_term) {
#   return(
#     egpd_param_df %>% as_tibble() %>% rename_with(., ~ as.character(reg_cols)) %>% 
#       mutate(time = c(1:t)) %>% 
#       pivot_longer(cols = c(1:all_of(r)), values_to = "effect", names_to = "region") %>%
#       mutate(region = as.numeric(region), covar = covar_term, linear = linear_term)
#   )
# }
# 
# for(i in seq_along(burn_beta_df_list)) {
#   temp_df <- burn_beta_df_list[[i]] %>%
#     mutate(dataset = 
#              case_when(grepl("g1_og", model_name) ~ "og",
#                        grepl("g1_sqrt", model_name) ~ "sqrt"),
#            combo = 
#              case_when(grepl("kappa-ri", model_name) ~ "kappa-ri",
#                        grepl("sigma-ri", model_name) ~ "sigma-ri",
#                        grepl("nu-ri", model_name) ~ "nu-ri",
#                        grepl("xi-ri", model_name) ~ "xi-ri",
#                        grepl("all-reg", model_name) ~ "all-reg",
#                        TRUE ~ model_name),
#            delta = 
#              case_when(grepl("0.81", model_name) ~ "0.81",
#                        grepl("0.9", model_name) ~ "0.9"),
#            model_name = paste0(dataset, "_", combo, "_", delta),
#            param = as.numeric(param),
#            coef = as.numeric(coef),
#            region = as.numeric(region))
#   model_name <- unique(temp_df$model_name)
#   if (unique(temp_df$dataset) == "og") {
#     X <- X_og
#   } else {
#     X <- X_sqrt
#   }
#   if (unique(temp_df$combo) == "all-reg") {
#     kappa_matrix <- egpd_param(temp_df, 1)
#     nu_matrix <- egpd_param(temp_df, 2)
#     xi_matrix <- egpd_param(temp_df, 3)
#     coef_df_list_kappa <- list()
#     coef_df_list_nu <- list()
#     coef_df_list_xi <- list()
#     for(k in seq_along(vars)) {
#       stored_df_kappa <- matrix(NA, t, r)
#       stored_df_nu <- matrix(NA, t, r)
#       stored_df_xi <- matrix(NA, t, r)
#       for(j in 1:r) {
#         stored_df_kappa[, j] <- X[j, , X_cols[[k]]] %*% kappa_matrix[X_cols[[k]], j]
#         stored_df_nu[, j] <- X[j, , X_cols[[k]]] %*% nu_matrix[X_cols[[k]], j]
#         stored_df_xi[, j] <- X[j, , X_cols[[k]]] %*% xi_matrix[X_cols[[k]], j]
#       }
#       coef_df_list_kappa[[k]] <- covar_effect(stored_df_kappa, vars[k], c(X[,,X_cols[[k]][2]]))
#       coef_df_list_nu[[k]] <- covar_effect(stored_df_nu, vars[k], c(X[,,X_cols[[k]][2]]))
#       coef_df_list_xi[[k]] <- covar_effect(stored_df_xi, vars[k], c(X[,,X_cols[[k]][2]]))
#     }
#     kappa_burns[[i]] <- bind_rows(coef_df_list_kappa) %>% as_tibble() %>% mutate(model = model_name) %>% left_join(., full_reg_key)
#     nu_burns[[i]] <- bind_rows(coef_df_list_nu) %>% as_tibble() %>% mutate(model = model_name) %>% left_join(., full_reg_key)
#     xi_burns[[i]] <- bind_rows(coef_df_list_xi) %>% as_tibble() %>% mutate(model = model_name) %>% left_join(., full_reg_key)
#   }  else if (unique(temp_df$combo) == "xi-ri") {
#     kappa_matrix <- egpd_param(temp_df, 1)
#     nu_matrix <- egpd_param(temp_df, 2)
#     coef_df_list_kappa <- list()
#     coef_df_list_nu <- list()
#     for(k in seq_along(vars)) {
#       stored_df_kappa <- matrix(NA, t, r)
#       stored_df_nu <- matrix(NA, t, r)
#       for(j in 1:r) {
#         stored_df_kappa[, j] <- X[j, , X_cols[[k]]] %*% kappa_matrix[X_cols[[k]], j]
#         stored_df_nu[, j] <- X[j, , X_cols[[k]]] %*% nu_matrix[X_cols[[k]], j]
#       }
#       coef_df_list_kappa[[k]] <- covar_effect(stored_df_kappa, vars[k], c(X[,,X_cols[[k]][2]]))
#       coef_df_list_nu[[k]] <- covar_effect(stored_df_nu, vars[k], c(X[,,X_cols[[k]][2]]))
#     }
#     kappa_burns[[i]] <- bind_rows(coef_df_list_kappa) %>% as_tibble() %>% mutate(model = model_name) %>% left_join(., full_reg_key)
#     nu_burns[[i]] <- bind_rows(coef_df_list_nu) %>% as_tibble() %>% mutate(model = model_name) %>% left_join(., full_reg_key)
#   } else if (unique(temp_df$combo) == "kappa-ri") {
#     nu_matrix <- egpd_param(temp_df, 1)
#     coef_df_list_nu <- list()
#     for(k in seq_along(vars)) {
#       stored_df_nu <- matrix(NA, t, r)
#       for(j in 1:r) {
#         stored_df_nu[, j] <- X[j, , X_cols[[k]]] %*% nu_matrix[X_cols[[k]], j]
#       }
#       coef_df_list_nu[[k]] <- covar_effect(stored_df_nu, vars[k], c(X[,,X_cols[[k]][2]]))
#     }
#     nu_burns[[i]] <- bind_rows(coef_df_list_nu) %>% as_tibble() %>% mutate(model = model_name) %>% left_join(., full_reg_key)    
#   } else {
#     kappa_matrix <- egpd_param(temp_df, 1)
#     coef_df_list_kappa <- list()
#     for(k in seq_along(vars)) {
#       stored_df_kappa <- matrix(NA, t, r)
#       for(j in 1:r) {
#         stored_df_kappa[, j] <- X[j, , X_cols[[k]]] %*% kappa_matrix[X_cols[[k]], j]
#       }
#       coef_df_list_kappa[[k]] <- covar_effect(stored_df_kappa, vars[k], c(X[,,X_cols[[k]][2]]))
#     }
#     kappa_burns[[i]] <- bind_rows(coef_df_list_kappa) %>% as_tibble() %>% mutate(model = model_name) %>% left_join(., full_reg_key)
#   }
# }
# 
# 
# for(i in seq_along(burn_names)) {
#   # count_beta <- get(burn_names[[i]])[[1]]
#   temp_df <- count_beta[[i]] %>%
#     apply(., c(2,3), median) %>% 
#     as_tibble() %>% pivot_longer(cols=everything(), names_to = "param_coef_reg", values_to = "value") %>%
#     group_by(param_coef_reg) %>% summarize(value = mean(value)) %>%
#     mutate(param_coef_reg = str_extract(param_coef_reg, "\\d{1},\\d{1,},\\d{1,}")) %>%
#     separate(., param_coef_reg, into=c("param", "coef", "region"), ",") %>%
#     mutate(param = case_when(
#       grepl("1", param) ~ "lambda",
#       grepl("2", param) ~ "pi",
#       grepl("3", param) ~ "delta",
#       TRUE ~ param),
#       param = as.factor(param),
#       coef = as.numeric(coef),
#       region = as.numeric(region)) %>%
#     filter(param == "lambda") %>% select(-param) %>% arrange(coef, region) %>%
#     pivot_wider(names_from = region, values_from = value) %>% select(-coef) %>% as.matrix()
#   coef_df_list <- list()
#   for(k in seq_along(vars)) {
#     stored_df <- matrix(NA, t, r)
#     for(j in 1:r) {
#       stored_df[, j] <- X[j, , X_cols[[k]]] %*% temp_df[X_cols[[k]], j]
#     }
#     coef_df_list[[k]] <- stored_df %>% 
#       as_tibble() %>% 
#       rename_with(., ~ reg_cols) %>% 
#       mutate(time = c(1:t)) %>% 
#       pivot_longer(cols = c(1:all_of(r)), values_to = "effect", names_to = "region") %>%
#       mutate(region = as.numeric(region), covar = vars[k], linear = c(X[,,X_cols[[k]][2]]))
#   }
#   lambda_counts[[i]] <- bind_rows(coef_df_list) %>% as_tibble() %>% mutate(model = burn_names[i]) %>% left_join(., full_reg_key)
# }
# 
# 
# ggplot(nu_burns[[3]], aes(x = linear, y = effect, group = region)) + 
#   geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE), show.legend = FALSE) +
#   facet_wrap(. ~ covar, scales = "free_x") + theme_minimal() + ggtitle(burn_names[3])
# file_name <- paste0("full-model/figures/model-comp/", burn_names[i], ".png")
# 
# 
# for(i in 1:length(kappa_burns)) {
#   if (!is.null(kappa_burns[[i]])) {
#     p <- ggplot(kappa_burns[[i]], aes(x = linear, y = effect, group = region)) + 
#       geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE)) +
#       facet_wrap(. ~ covar, scales = "free_x") + theme_minimal() + ggtitle(paste0("kappa_", unique(kappa_burns[[i]]$model)))
#     file_name <- paste0("full-model/figures/model-comp/g1_kappa_", unique(kappa_burns[[i]]$model), ".png")
#     ggsave(file_name, p, dpi = 320, type = "cairo", bg = "white")
#   }
#   else {
#     print("element of list is empty")
#   }
# }
# 
# # for(i in 1:length(kappa_burns)) {
# #   if (!is.null(kappa_burns[[i]])) {
# #     print(paste0(i, "_", unique(kappa_burns[[i]]$model)))
# #   }
# #   else {
# #     print("element of list is empty")
# #   }
# # }
# 
# for(i in 1:length(nu_burns)) {
#   if (!is.null(nu_burns[[i]])) {
#     p <- ggplot(nu_burns[[i]], aes(x = linear, y = effect, group = region)) + 
#       geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE), show.legend = FALSE) +
#       facet_wrap(. ~ covar, scales = "free_x") + theme_minimal() + ggtitle(paste0("nu_", unique(nu_burns[[i]]$model)))
#     file_name <- paste0("full-model/figures/model-comp/g1_nu_", unique(nu_burns[[i]]$model), ".png")
#     ggsave(file_name, p, dpi = 320, type = "cairo", bg = "white")
#   }
#   else {
#     print("element of list is empty")
#   }
# }
# 
# for(i in 1:length(xi_burns)) {
#   if (!is.null(xi_burns[[i]])) {
#     p <- ggplot(xi_burns[[i]], aes(x = linear, y = effect, group = region)) + 
#       geom_line(aes(linetype = NA_L1CODE, color = NA_L2CODE), show.legend = FALSE) +
#       facet_wrap(. ~ covar, scales = "free_x") + theme_minimal() + ggtitle(paste0("xi_", unique(xi_burns[[i]]$model)))
#     file_name <- paste0("full-model/figures/model-comp/g1_xi_", unique(xi_burns[[i]]$model), ".png")
#     ggsave(file_name, p, dpi = 320, type = "cairo", bg = "white")
#   }
#   else {
#     print("element of list is empty")
#   }
# }
# 
# 
# poster_kappa <- kappa_burns[[17]] %>% filter(covar != "log_housing_density", covar != "pr", covar != "rmin", covar != "tmmx") %>%
#   mutate(covar = case_when(
#     covar == 'vs' ~ 'Wind speed',
#     covar == 'prev_12mo_precip' ~ 'Precipitation: 12 month'),
#     covar = factor(covar, levels = c("Wind speed", "Precipitation: 12 month")))
# 
# p <- ggplot(poster_kappa, aes(x=linear, y=effect, group = region)) + 
#   geom_line(aes(linetype=NA_L1CODE, color = NA_L2CODE), show.legend = FALSE) + 
#   facet_grid(. ~ covar) + xlab("Linear term") + ylab(expression("Partial effect on"~kappa))+
#   theme_classic(base_size = 25)
# 
# ggsave("~/Desktop/research/posters-presentations/prelim_kappa.png", dpi = 700, type = "cairo", bg="white",
#        height = 11, width = 17)

## following code is for model evaluation (scoring) ---------
files <- paste0("full-model/fire-sims/burns/", 
                    list.files(path = "full-model/fire-sims/burns/", 
                               pattern = ".csv", recursive = TRUE))
gq_files <- c(files[which(grepl("gq_", files))], files[which(grepl("sigma-reg", files))])
gq_files <- gq_files[-which(grepl("nu-ri", gq_files))]
gq_files <- gq_files[-which(grepl("0.9", gq_files))]

g1_files <- gq_files[which(grepl("g1", gq_files))]
g2_files <- gq_files[which(grepl("g2", gq_files))]
lognorm_files <- gq_files[which(grepl("lognorm", gq_files))]
g1_files_sigma <- c(g1_files[which(grepl("sigma-reg", g1_files))], g1_files[which(grepl("sigma-ri", g1_files))])
g1_files_sigma <- g1_files_sigma[-which(grepl("14Jun2023", g1_files_sigma))]
# g2_files_03jun <- g2_files[which(grepl("03Jun", g2_files))]
# g2_files_03jun <- g2_files_03jun[which(grepl("sigma-ri", g2_files_03jun))]
g2_files_07may <- g2_files[which(grepl("07May", g2_files))]
g2_files_07may <- g2_files_07may[c(which(grepl("sqrt_xi-ri", g2_files_07may)), which(grepl("sigma-ri", g2_files_07may)), which(grepl("sqrt_kappa-ri", g2_files_07may)))]
# g2_fit_mixed <- c(g2_files_03jun, g2_files_07may)
lognorm_files <- lognorm_files[-which(grepl("sqrt_all-reg", lognorm_files))]
all_gq_burns <- c(g1_files_sigma, g2_files_07may, lognorm_files)

# gq_fits <- c(gq_fits, sigma_reg_fits)
nfits <- length(all_gq_burns)/3
gq_fit_groups <- vector(mode = "list", nfits)
for(i in 1:nfits) {
  gq_fit_groups[[i]] <- all_gq_burns[(3*i-2):(3*i)]
}

gq_mod_names <- lapply(gq_fit_groups, function(x) str_remove(str_remove(str_remove(basename(x[1]), "_\\d{2}\\w{3}2023_\\d{4}_\\d{1}.csv"), 
                                                                     "cfcns_"), 
                                                          "gq_\\d{2}\\w{3}2023_\\d{4}_")) %>% unlist()

# ensure every model name is unique
sum(duplicated(gq_mod_names))

# pull indices of individual sets of generated quantities to separate within function (applicable to "gq" csv files only)
one_fit <- read_cmdstan_csv(gq_fit_groups[[19]])
train_ll_idx <- which(grepl("train_loglik", variables(one_fit$generated_quantities)))
holdout_ll_idx <- which(grepl("holdout_loglik", variables(one_fit$generated_quantities)))
train_twcrps_idx <- which(grepl("train_twcrps", variables(one_fit$generated_quantities)))
holdout_twcrps_idx <- which(grepl("holdout_twcrps", variables(one_fit$generated_quantities)))

extraction <- function(file_group, model_name) {
  if(file.info(file_group[1])$size/1e6 > 300){
    model_object <- as_cmdstan_fit(file_group)
    train_loglik <- model_object$draws(variables = "train_loglik")
    holdout_loglik <- model_object$draws(variables = "holdout_loglik")
    train_twcrps <- model_object$draws(variables = "train_twcrps")
    holdout_twcrps <- model_object$draws(variables = "holdout_twcrps")
    rm(model_object)
  } else {
    gen_quants <- read_cmdstan_csv(file_group)
    train_loglik <- gen_quants$generated_quantities[,,train_ll_idx]
    holdout_loglik <- gen_quants$generated_quantities[,,holdout_ll_idx]
    train_twcrps <- gen_quants$generated_quantities[,,train_twcrps_idx]
    holdout_twcrps <- gen_quants$generated_quantities[,,holdout_twcrps_idx]
    rm(gen_quants)
  }
  temp <- list(train_loglik, holdout_loglik, train_twcrps, holdout_twcrps)
  names(temp) <- c("train_loglik", "holdout_loglik", "train_twcrps", "holdout_twcrps")
  assign(model_name, temp, parent.frame())
  gc()
}

rm(one_fit)
gc()

# remove gen-quant file that didn't complete (will investigate later)
gq_fit_groups <- c(gq_fit_groups[1:8], list(gq_fit_groups[[9]][2:3]), gq_fit_groups[10:nfits])

for(i in 1:nfits) {
  extraction(gq_fit_groups[[i]], gq_mod_names[i])
}

# also excluding lognorm with mu-ri
gq_mod_names <- gq_mod_names[-which(grepl("mu-ri", gq_mod_names))]

nfits <- length(gq_mod_names)
train_loglik_list <- vector("list", nfits)
holdout_loglik_list <- vector("list", nfits)
train_twcrps_list <- vector("list", nfits)
holdout_twcrps_list <- vector("list", nfits)

for(i in seq_along(gq_mod_names)) {
  model_string <- str_split(gq_mod_names[i], pattern = "_")[[1]]
  model_string <- if(length(model_string) > 4) model_string[-4] else model_string
  train_loglik_list[[i]] <- get(gq_mod_names[i])[["train_loglik"]] %>%
    as_draws_df() %>%
    select(-c(".iteration", ".chain")) %>% 
    pivot_longer(cols = !".draw") %>%
    rename(draw = ".draw") %>%
    group_by(draw) %>% 
    summarize(loglik = sum(value)) %>%
    mutate(model = model_string[1],
           dataset = model_string[2],
           params = model_string[3],
           stepsize = model_string[4],
           train = TRUE)
  holdout_loglik_list[[i]] <- get(gq_mod_names[i])[["holdout_loglik"]] %>%
    as_draws_df() %>%
    select(-c(".iteration", ".chain")) %>% 
    pivot_longer(cols = !".draw") %>%
    rename(draw = ".draw") %>%
    group_by(draw) %>% 
    summarize(loglik = sum(value)) %>%
    mutate(model = model_string[1],
           dataset = model_string[2],
           params = model_string[3],
           stepsize = model_string[4],
           train = FALSE)
  train_twcrps_list[[i]] <- get(gq_mod_names[i])[["train_twcrps"]] %>%
    as_draws_df() %>%
    select(-c(".iteration", ".chain")) %>% 
    pivot_longer(cols = !".draw") %>%
    rename(draw = ".draw") %>%
    group_by(draw) %>% 
    summarize(twcrps = mean(value)) %>%
    mutate(model = model_string[1],
           dataset = model_string[2],
           params = model_string[3],
           stepsize = model_string[4],
           train = TRUE)
  holdout_twcrps_list[[i]] <- get(gq_mod_names[i])[["holdout_twcrps"]] %>%
    as_draws_df() %>%
    select(-c(".iteration", ".chain")) %>% 
    pivot_longer(cols = !".draw") %>%
    rename(draw = ".draw") %>%
    group_by(draw) %>% 
    summarize(twcrps = mean(value)) %>%
    mutate(model = model_string[1],
           dataset = model_string[2],
           params = model_string[3],
           stepsize = model_string[4],
           train = FALSE)
}

train_loglik <- bind_rows(train_loglik_list)
holdout_loglik <- bind_rows(holdout_loglik_list)
train_twcrps <- bind_rows(train_twcrps_list)
holdout_twcrps <- bind_rows(holdout_twcrps_list)

## log-likelihood aggregation and comparisons --------
ll_full <- train_loglik %>% 
  full_join(holdout_loglik) %>% 
  mutate(params = case_when(stepsize == "sigma-reg" ~ paste0(params, "_", stepsize),
                            TRUE ~ params),
         fullname = paste0(model, "_", dataset, "_", params)) %>%
  select(-stepsize)

train_ll_sort <- ll_full %>% filter(train == TRUE) %>% 
  group_by(fullname, model, dataset, params) %>% 
  summarize(med_train_ll = median(loglik)) %>% arrange(-med_train_ll)
top_mod_train <- as.character(train_ll_sort$fullname[1])

test_ll_sort <- ll_full %>% filter(train == FALSE) %>% 
  group_by(fullname, model, dataset, params) %>% 
  summarize(med_test_ll = median(loglik)) %>% arrange(-med_test_ll)
top_mod_test <- as.character(test_ll_sort$fullname[1])

ll_comp_train_full <- ll_full %>% filter(train == TRUE) %>% 
  select(c(draw, fullname, loglik)) %>% 
  pivot_wider(names_from = fullname, values_from = loglik, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_train))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(med_diff = median(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% arrange(-med_diff)

ll_comp_test_full <- ll_full %>% filter(train == FALSE) %>% 
  select(c(draw, fullname, loglik)) %>% 
  pivot_wider(names_from = fullname, values_from = loglik, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_test))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(med_diff = median(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% arrange(-med_diff)
# -----------

## twCRPS aggregation and comparions ---------
twcrps_full <- train_twcrps %>%
  full_join(holdout_twcrps) %>%
  mutate(params = case_when(stepsize == "sigma-reg" ~ paste0(params, "_", stepsize),
                            TRUE ~ params),
         fullname = paste0(model, "_", dataset, "_", params)) %>%
  select(-stepsize)

limits_twcrps_full <- twcrps_full %>% reframe(limits = quantile(twcrps, c(0.05,0.95), na.rm = TRUE))
twcrps_boxplot_full <- twcrps_full %>%
  ggplot(aes(fullname, twcrps, color = train)) +
  geom_boxplot(outlier.shape = NA) + scale_y_continuous(limits = limits_twcrps_og$limits) +
  theme_minimal()
ggsave("full-model/figures/model-comp/twcrps_g1_full_18jun2023.png", plot = twcrps_boxplot_full,
       dpi = 320, bg = "white")

train_twcrps_sort <- twcrps_full %>% filter(train == TRUE) %>% 
  group_by(fullname, model, dataset, params) %>% 
  summarize(mean_train_twcrps = mean(twcrps, na.rm = TRUE)) %>% arrange(mean_train_twcrps)
top_mod_train_tw <- as.character(train_twcrps_sort$fullname[1])

test_twcrps_sort <- twcrps_full %>% filter(train == FALSE) %>% 
  group_by(fullname, model, dataset, params) %>% 
  summarize(mean_test_twcrps = mean(twcrps, na.rm = TRUE)) %>% arrange(mean_test_twcrps)
top_mod_test_tw <- as.character(test_twcrps_sort$fullname[1])

twcrps_comp_train_full <- twcrps_full %>% filter(train == TRUE) %>% 
  select(c(draw, fullname, twcrps)) %>% 
  pivot_wider(names_from = fullname, values_from = twcrps, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_train_tw))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(mean_diff = mean(value, na.rm = TRUE), sd_diff = sd(value, na.rm = TRUE)) %>% arrange(mean_diff)

twcrps_comp_test_full <- twcrps_full %>% filter(train == FALSE) %>% 
  select(c(draw, fullname, twcrps)) %>% 
  pivot_wider(names_from = fullname, values_from = twcrps, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_test_tw))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(mean_diff = mean(value, na.rm = TRUE), sd_diff = sd(value, na.rm = TRUE)) %>% arrange(mean_diff)

