##############################################################################
##  This script produces two partial effects plots:                         ##
##  (1) climate and housing covariates on rate parameter of count submodel  ##
##  (2) fire indices and housing covariates on kappa parameter of           ##
##      sizes submodel                                                      ##
##############################################################################

library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(tidyverse)
library(stringr)
library(lubridate)
library(sf)
library(RColorBrewer)
library(patchwork)
library(colorspace)
source("./data/load_eco.R")

## load in ecoregion shapes ----------------
ecoregion_shp <- load_ecoregions()
ecoregion_df <- ecoregion_shp %>% st_as_sf()
levels_l1 <- levels(as.factor(as.numeric(ecoregion_df$NA_L1CODE)))
levels_l2 <- levels(as.factor(as.numeric(ecoregion_df$NA_L2CODE)))
ecoregion_df <- ecoregion_df %>% 
  mutate(NA_L2CODE = factor(NA_L2CODE, levels = levels_l2),
         NA_L1CODE = factor(NA_L1CODE, levels = levels_l1),
         NA_L3CODE = as.factor(NA_L3CODE),
         NA_L1NAME = as.factor(str_to_title(NA_L1NAME)))

## read in L3 ecoregion key
region_key <- readRDS(file = "./data/processed/region_key.rds") %>%
  mutate(NA_L2CODE = factor(NA_L2CODE, levels = levels_l2),
         NA_L1CODE = factor(NA_L1CODE, levels = levels_l1),
         NA_L3CODE = as.factor(NA_L3CODE),
         NA_L1NAME = as.factor(str_to_title(NA_L1NAME)))

## read in posterior beta values --------------------------------------------
beta_count <- readRDS("./figures/mcmc_draws/beta_count.RDS")
beta_size <- readRDS("./figures/mcmc_draws/beta_size.RDS")
r <- 84
t <- 12*21 

# partial effects on lambda
stan_data <- readRDS("./data/stan_lists/data_joint.RDS")
X <- stan_data$X_train_count
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

lambda_betas <- beta_count %>%
  separate_wider_delim(cols = "name", delim = ",", names = c("coef", "region")) %>%
  mutate(coef = as.numeric(gsub("beta_count\\[", "", coef)),
         region = as.numeric(gsub("\\]", "", region))) %>%
  group_by(region, coef) %>% summarize(med_val = median(value)) %>% ungroup() %>%
  pivot_wider(names_from = "region", values_from = "med_val") %>% select(-coef) %>% as.matrix()

# rescale data
un_std_data <- readRDS("./data/processed/un_std_all.RDS")
X_unstd <- X
for(k in seq_along(vars)) {
  var_mean <- un_std_data %>% filter(variable == vars[[k]]) %>% select(mean) %>% as.numeric()
  var_sd <- un_std_data %>% filter(variable == vars[[k]]) %>% select(sd) %>% as.numeric()
  X_unstd[,,X_cols[[k]][2]] <- X[,,X_cols[[k]][2]] * var_sd + var_mean
}

reg_cols <- region_key$region
covar_effect <- function(param_df, covar_term, linear_term) {
  return(
    param_df %>% as_tibble() %>% rename_with(., ~ as.character(reg_cols)) %>%
      mutate(time = c(1:t)) %>%
      pivot_longer(cols = c(1:all_of(r)), values_to = "effect", names_to = "region") %>%
      mutate(region = as.numeric(region), covar = covar_term, linear = linear_term)
  )
}

coef_df_list_lambda <- list()
for(k in seq_along(vars)) {
  var_mean <- un_std_data %>% filter(variable == vars[[k]]) %>% select(mean) %>% as.numeric()
  var_sd <- un_std_data %>% filter(variable == vars[[k]]) %>% select(sd) %>% as.numeric()
  stored_df_lambda <- matrix(NA, t, r)
  for(j in 1:r) {
    stored_df_lambda[, j] <- X[j, , X_cols[[k]][3:7]] %*% lambda_betas[X_cols[[k]][3:7], j] +
      X_unstd[j, , X_cols[[k]][2]] * lambda_betas[X_cols[[k]][2], j]/var_sd +
      lambda_betas[1, j] - (lambda_betas[X_cols[[k]][2], j] * var_mean)/var_sd
  }
  coef_df_list_lambda[[k]] <- covar_effect(stored_df_lambda, vars[k], c(X_unstd[,,X_cols[[k]][2]]))
}

lambda_effects <- bind_rows(coef_df_list_lambda) %>% as_tibble() %>% 
  left_join(., region_key) %>% 
  mutate(NA_L2CODE = factor(NA_L2CODE, levels = levels_l2))
lambda_plot <- lambda_effects %>% 
  mutate(linear = case_when(covar == 'tmmx' ~ linear - 273.15,
                            TRUE ~ linear),
         covar = case_when(covar == 'log_housing_density' ~ "log(housing density (units/sq. km))",
                           covar == 'pr' ~ 'Precipitation (mm), same month',
                           covar == 'prev_12mo_precip' ~ 'Precipitation (mm), past 12months',
                           covar == 'rmin' ~ 'Min. relative humidity (%)',
                           covar == 'tmmx' ~ 'Max. air temperature (C)',
                           covar == 'vs' ~ 'Wind speed (m/s)',
                           TRUE ~ covar)) %>% 
  ggplot(aes(x = linear, y = effect, group = region)) + 
  geom_line(aes(color = NA_L2CODE)) +
  facet_wrap(. ~ covar, scales = "free_x") + 
  theme_classic() + 
  theme(legend.position = "none") +
  ylab("Partial effect") + xlab("")
ggsave("./figures/paper_figures/lambda_partial_effects.pdf", 
       lambda_plot, 
       width = 7, 
       height = 4,
       bg = 'transparent')
knitr::plot_crop("./figures/paper_figures/lambda_partial_effects.pdf")

# determine regions with particularly prominent effects
lambda_effects %>% filter(covar == 'vs', effect > 0.5, linear > 5) %>% group_by(NA_L2CODE, NA_L1NAME, NA_L3NAME) %>% count() %>% arrange(-n)
lambda_effects %>% filter(covar == 'pr', effect < -2.5, linear > 290) %>% group_by(NA_L2CODE, NA_L1NAME, NA_L3NAME) %>% count() %>% arrange(-n)
lambda_effects %>% filter(covar == 'pr', effect > 0, linear > 200) %>% group_by(NA_L2CODE, NA_L1NAME, NA_L3NAME) %>% count() %>% arrange(-n)

# partial effects on kappa
X <- stan_data$X_train_size
vars <- c('log_housing_density', 'erc', 'fwi')
X_covar <- c()
X_cols <- vector("list", length(vars))
start <- 2
for(i in seq_along(vars)) {
  X_covar[i] <- paste0("X_", vars[i])
  X_cols[[i]] <- c(1, start:(start+5))
  start = start + 6
}

kappa_betas <- beta_size %>%
  separate_wider_delim(cols = "name", delim = ",", names = c("coef", "region")) %>%
  mutate(coef = as.numeric(gsub("beta_size\\[", "", coef)),
         region = as.numeric(gsub("\\]", "", region))) %>%
  group_by(region, coef) %>% summarize(med_val = median(value)) %>% ungroup() %>%
  pivot_wider(names_from = "region", values_from = "med_val") %>% select(-coef) %>% as.matrix()

# rescale data
X_unstd <- X
for(k in seq_along(vars)) {
  var_mean <- un_std_data %>% filter(variable == vars[[k]]) %>% select(mean) %>% as.numeric()
  var_sd <- un_std_data %>% filter(variable == vars[[k]]) %>% select(sd) %>% as.numeric()
  X_unstd[,,X_cols[[k]][2]] <- X[,,X_cols[[k]][2]] * var_sd + var_mean
}

coef_df_list_kappa <- list()
for(k in seq_along(vars)) {
  var_mean <- un_std_data %>% filter(variable == vars[[k]]) %>% select(mean) %>% as.numeric()
  var_sd <- un_std_data %>% filter(variable == vars[[k]]) %>% select(sd) %>% as.numeric()
  stored_df_kappa <- matrix(NA, t, r)
  for(j in 1:r) {
    stored_df_kappa[, j] <- X[j, , X_cols[[k]][3:7]] %*% kappa_betas[X_cols[[k]][3:7], j] +
      X_unstd[j, , X_cols[[k]][2]] * kappa_betas[X_cols[[k]][2], j]/var_sd +
      kappa_betas[1, j] - (kappa_betas[X_cols[[k]][2], j] * var_mean)/var_sd
  }
  coef_df_list_kappa[[k]] <- covar_effect(stored_df_kappa, vars[k], c(X_unstd[,,X_cols[[k]][2]]))
}

kappa_effects <- bind_rows(coef_df_list_kappa) %>% as_tibble() %>% 
  left_join(., region_key) %>% 
  mutate(NA_L2CODE = factor(NA_L2CODE, levels = levels_l2))
kappa_plot <- kappa_effects %>% 
  mutate(covar = case_when(covar == 'log_housing_density' ~ "log(housing density (units/sq. km))",
                           covar == 'erc' ~ 'Energy Release Component',
                           covar == 'fwi' ~ 'Fire Weather Index',
                           TRUE ~ covar)) %>% 
  ggplot(aes(x = linear, y = effect, group = region)) + 
  geom_line(aes(color = NA_L2CODE)) +
  facet_wrap(. ~ covar, scales = "free_x") + theme_classic() + 
  theme(legend.position = "none") +
  ylab("Partial effect") + xlab("")
ggsave("./figures/paper_figures/kappa_partial_effects.pdf", 
       kappa_plot, 
       width = 7, 
       height = 2,
       bg = 'transparent')
knitr::plot_crop("./figures/paper_figures/kappa_partial_effects.pdf")

# determine regions with particularly prominent effects
kappa_effects %>% filter(covar == 'erc', effect > 0, linear < 10) %>% 
  group_by(NA_L2CODE, NA_L1NAME, NA_L3NAME) %>% 
  count() %>% arrange(-n)
kappa_effects %>% filter(covar == 'erc', effect > -0.5, effect < 0, linear < 10) %>% 
  group_by(NA_L2CODE, NA_L1NAME, NA_L3NAME) %>% 
  count() %>% arrange(-n)
kappa_effects %>% filter(covar == 'erc', effect > 0.25, linear > 60) %>% 
  group_by(NA_L2CODE) %>% 
  count() %>% arrange(-n) 
kappa_effects %>% filter(covar == 'fwi', linear > 40, effect > 0) %>% 
  group_by(NA_L2CODE) %>% count() %>% arrange(-n)

kappa_effects |> filter(covar == 'fwi', NA_L2CODE == 12.1) |> arrange(-effect)
kappa_effects |> filter(covar == 'erc', NA_L2CODE == 12.1) |> arrange(-effect)
