library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(tidyverse)
library(stringr)
library(posterior)
library(lubridate)
library(sf)
library(classInt)
library(RColorBrewer)

## read in region key 
region_key <- readRDS(file = "./full-model/data/processed/region_key.rds")
full_reg_key <- as_tibble(region_key) %>% 
  mutate(region = c(1:84),
         NA_L2CODE = as.factor(NA_L2CODE),
         NA_L1CODE = as.factor(NA_L1CODE),
         NA_L3CODE = as.factor(NA_L3CODE),
         NA_L1NAME = as.factor(str_to_title(NA_L1NAME)))

## read in ecoregion shapes
ecoregions <- read_rds(file = "ecoregions.RDS")
ecoregions_geom <- ecoregions %>% filter(!NA_L2NAME == "UPPER GILA MOUNTAINS (?)") %>% 
  mutate(NA_L2CODE = as.factor(NA_L2CODE),
         NA_L1CODE = as.factor(NA_L1CODE),
         NA_L3CODE = as.factor(NA_L3CODE),
         NA_L1NAME = as.factor(str_to_title(NA_L1NAME)))

# egpd functions - G1
pegpd <- function(y, kappa, sigma, xi) (1 - (1 + xi * (y/sigma))^(-1/xi))^kappa
qegpd <- function(p, kappa, sigma, xi) (sigma/xi) * ( (1 - p^(1/kappa) )^-xi - 1)
regpd <- function(n, kappa, sigma, xi) { # truncated version
  u <- runif(n)
  u_adj <- u * (1 - pegpd(1.001, kappa, sigma, xi)) + pegpd(1.001, kappa, sigma, xi)
  return(qegpd(u_adj, kappa, sigma, xi))
}
exp_egpd <- function(n, kappa, sigma, xi) {
  mcmc_rng <- regpd(n, kappa, sigma, xi) 
  return(mean(mcmc_rng[is.finite(mcmc_rng)]))
}

date_seq <- seq(as.Date("1990-01-01"), by = "1 month", length.out = 372) %>% as_tibble() %>% rename(date = value)
time_df <- date_seq %>% mutate(time = 1:372)

all_years <- 1990:2020
first_five <- 1990:1994
last_five <- 2020:2016
test_years <- sort(c(first_five, last_five))
train_years <- setdiff(all_years, test_years)

time_df <- time_df %>% mutate(train = case_when(year(date) %in% train_years ~ TRUE,
                                                year(date) %in% test_years ~ FALSE))
test_tmpts <- time_df %>% 
  filter(train == FALSE) %>% 
  select(time) %>% 
  rowid_to_column(var = "model_tmpt") %>% 
  rename(true_tmpt = time)

train_tmpts <- time_df %>% 
  filter(train == TRUE) %>% 
  select(time) %>% 
  rowid_to_column(var = "model_tmpt") %>% 
  rename(true_tmpt = time)

# extract sigma and xi 
g1_files <- paste0("full-model/fire-sims/burns/g1/csv-fits/",
                   list.files(path = "full-model/fire-sims/burns/g1/csv-fits", pattern = "erc_12Sep"))
rand_int_draws <- read_cmdstan_csv(g1_files, variables = "ri_init")$post_warmup_draws
rand_int <- rand_int_draws %>%
  as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>%
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  separate_wider_delim(cols = "name", delim = ",", names = c("param", "region")) %>%
  mutate(param = as.character(gsub("ri_init\\[", "", param)),
         region = as.numeric(gsub("\\]", "", region)),
         param = as.character(param),
         param = case_when(param == "1" ~ "sigma",
                           param == "2" ~ "xi",
                           TRUE ~ param),
         value = exp(value)) %>%
  pivot_wider(names_from = param, values_from = value)
saveRDS(rand_int, file = paste0("full-model/fire-sims/burns/g1/best_fit_", "rand_int", ".RDS"))
rm(rand_int_draws)

# extract kappa values from gen quant block (so it's all timepoints) of best G1 model
g1_gq_files <- paste0("full-model/fire-sims/model_comparison/extracted_values/",
                      list.files(path = "full-model/fire-sims/model_comparison/extracted_values/", 
                                 pattern = "g1_sigma-ri_xi-ri_erc_12"))
g1_gq_files <- g1_gq_files[grepl("long", g1_gq_files)]
g1_gq_files <- g1_gq_files[grepl(".csv", g1_gq_files)]
kappa_draws <- read_cmdstan_csv(g1_gq_files, variables = "reg_full")$generated_quantities
kappa <- kappa_draws %>%
  as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>%
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  separate_wider_delim(cols = "name", delim = ",", names = c("param", "time", "region")) %>%
  mutate(param = as.numeric(gsub("reg_full\\[", "", param)),
         time = as.numeric(time),
         region = as.numeric(gsub("\\]", "", region)),
         kappa = exp(value)) %>% select(-value)
kappa <- kappa %>% select(-param)
saveRDS(kappa, file = paste0("full-model/fire-sims/burns/g1/best_fit_", "kappa", ".RDS"))
rm(kappa_draws)
gc()

# extract count model params
zinb_files <- paste0("full-model/fire-sims/counts/zinb_er/csv-fits/",
                   list.files(path = "full-model/fire-sims/counts/zinb_er/csv-fits", pattern = "pi-ri_climate"))
lambda_train_draws <- read_cmdstan_csv(zinb_files, variables = "lambda")$post_warmup_draws
lambda_train <- lambda_train_draws %>%
  as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>%
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  separate_wider_delim(cols = "name", delim = ",", names = c("time", "region")) %>%
  mutate(time = as.numeric(gsub("lambda\\[", "", time)),
         region = as.numeric(gsub("\\]", "", region)),
         param = "lambda") %>%
  pivot_wider(names_from = param, values_from = value)
lambda_train <- lambda_train %>% left_join(train_tmpts, join_by(time == model_tmpt))
lambda_train <- lambda_train %>% select(-time) %>% rename(time = true_tmpt)

lambda_hold_draws <- read_cmdstan_csv(zinb_files, variables = "lambda_hold")$post_warmup_draws
lambda_hold <- lambda_hold_draws %>%
  as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>%
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  separate_wider_delim(cols = "name", delim = ",", names = c("time", "region")) %>%
  mutate(time = as.numeric(gsub("lambda_hold\\[", "", time)),
         region = as.numeric(gsub("\\]", "", region)),
         param = "lambda") %>%
  pivot_wider(names_from = param, values_from = value) %>%
  left_join(test_tmpts, join_by(time == model_tmpt)) %>%
  select(-time) %>% 
  rename(time = true_tmpt)
lambda <- rbind(lambda_hold, lambda_train)
saveRDS(lambda, file = paste0("full-model/fire-sims/counts/zinb_er/best_fit_", "lambda", ".RDS"))

pi_prob_draws <- read_cmdstan_csv(zinb_files, variables = "pi_prob")$post_warmup_draws
pi_prob <- pi_prob_draws %>%
  as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>%
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  mutate(region = as.character(gsub("pi_prob\\[", "", name)),
         region = as.numeric(gsub("\\]", "", region)),
         name = "pi") %>%
  pivot_wider(names_from = name, values_from = value)
saveRDS(pi_prob, file = paste0("full-model/fire-sims/counts/zinb_er/best_fit_", "pi_prob", ".RDS"))

delta_draws <- read_cmdstan_csv(zinb_files, variables = "delta")$post_warmup_draws
delta <- delta_draws %>%
  as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>%
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  mutate(region = as.character(gsub("delta\\[", "", name)),
         region = as.numeric(gsub("\\]", "", region)),
         name = "delta") %>%
  pivot_wider(names_from = name, values_from = value)
saveRDS(delta, file = paste0("full-model/fire-sims/counts/zinb_er/best_fit_", "delta", ".RDS"))

burn_preds_gen <- function(pi_prob, delta, lambda, kappa, sigma, xi) {
  burn_draws = rep(NA, 500)
  pi <- exp(pi_prob)/(1+exp(pi_prob))
  for(i in 1:500) {
    zero <- rbinom(1, 1, pi)
    count_draw <- (1-zero) * rnbinom(1, size = delta, mu = exp(lambda))
    if(count_draw == 0) {
      burn_draws[i] <- 0
    } else {
      temp <- regpd(count_draw, kappa, sigma, xi)
      burn_draws[i] = sum(temp[is.finite(temp)])
    }
  }
  return(mean(burn_draws, na.rm = TRUE))
}

kappa <- readRDS(paste0("full-model/fire-sims/burns/g1/best_fit_", "kappa", ".RDS"))
rand_int <- readRDS(paste0("full-model/fire-sims/burns/g1/best_fit_", "rand_int", ".RDS"))
lambda <- readRDS(paste0("full-model/fire-sims/counts/zinb_er/best_fit_", "lambda", ".RDS"))
delta <- readRDS(paste0("full-model/fire-sims/counts/zinb_er/best_fit_", "delta", ".RDS"))
pi_prob <- readRDS(paste0("full-model/fire-sims/counts/zinb_er/best_fit_", "pi_prob", ".RDS"))

burn_preds <- kappa %>% left_join(rand_int) %>% left_join(lambda) %>% left_join(pi_prob) %>% left_join(delta)
burn_preds <- burn_preds %>%
  mutate(preds = burn_preds_gen(pi, delta, lambda, kappa, sigma, xi))

# // generate burn area predictions, based on how many fires occurred (draw from ZINB)
# for (r in 1:R) {
#   for (t in 1:T_all) {
#     vector[500] burn_draws;
#     real sigma = exp(ri_init[1][r]);
#     real xi = exp(ri_init[2][r]);
#     real kappa = exp(reg_full[t, r]);
#     for (i in 1:500) {
#       int zero = bernoulli_logit_rng(pi_prob[r]);
#       int count_draw = (1 - zero) * neg_binomial_2_log_rng(lambda_full[t, r], delta[r]);
#       if (count_draw == 0) {
#         burn_draws[i] = 0; 
#       } else {
#         burn_draws[i] = sum(egpd_rng(count_draw, y_min, sigma, xi, kappa));
#       }
#     }
#     burn_pred[t, r] = mean(burn_draws);
#   }
# }


## create time series plots of phi values for kappa and lambda -------
date_seq <- seq(as.Date("1990-01-01"), by = "1 month", length.out = 372) %>% as_tibble() %>% rename(date = value)
time_df <- date_seq %>% mutate(time = 1:372)


## create boxplot of predicted burn area for entire US annually for all timepoints (holdout and training) and overlay truth ------
# burn_params <- kappa %>% 
#   left_join(rand_int)
# 
# burn_preds <- readRDS("full-model/figures/paper/burn_preds_df.RDS")
# burn_preds <- burn_preds %>% 
#   left_join(delta) %>% 
#   left_join(lambda) %>% 
#   left_join(pi_prob) %>% 
#   mutate(condl_preds = preds * fire_occur(pi, delta, lambda))
# burn_preds <- burn_preds %>% select(-c("delta", "lambda", "pi"))
# burn_preds <- burn_params %>%
#   mutate(preds = purrr::pmap_dbl(list(500, kappa, sigma, xi), med_egpd)) %>%
#   left_join(full_reg_key) %>%
#   left_join(time_df) %>%
#   mutate(year = year(date))
# saveRDS(burn_preds, file = "full-model/figures/paper/burn_preds_df.RDS")