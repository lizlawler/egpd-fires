args <- commandArgs(trailingOnly=TRUE)
chain <- args[1]

library(cmdstanr)
set_cmdstan_path(path = "/projects/eslawler@colostate.edu/software/anaconda/envs/lawler/bin/cmdstan") # this is only relevant to Alpine
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(tidyr)
library(dplyr)
library(tibble)
library(stringr)
library(posterior)
library(lubridate)
library(extraDistr)

# egpd functions - G1
egpd_cdf <- function(y, sigma, xi, kappa) {
  lcdf <- kappa * pgpd(y, mu = 0, sigma = sigma, xi = xi, log.p = TRUE)
  return(exp(lcdf))
}
egpd_icdf <- function(u, sigma, xi, kappa) {
  p <- exp( exp(-log(kappa)) * log(u))
  return(qgpd(p, mu = 0, sigma = sigma, xi = xi))
}
egpd_rng <- function(n, sigma, xi, kappa) {
  cst <- egpd_cdf(1.001, sigma, xi, kappa)
  u <- runif(n)
  rng_val <- rep(NA, n)
  for(i in 1:n) {
    u_adj = u[i] * (1-cst) + cst
    rng_val[i] = egpd_icdf(u_adj, sigma, xi, kappa)
  }
  return(rng_val)
}

burn_preds_gen <- function(n = 500, pi_prob, delta, lambda, sigma, xi, kappa) {
  pi <- exp(pi_prob)/(1+exp(pi_prob))
  burns_rng <- function(pi, delta, lambda, sigma, xi, kappa) {
    zero <- rbinom(1, 1, pi)
    count_draw <- (1-zero) * rnbinom(1, size = delta, mu = exp(lambda))
    if(count_draw == 0) {
      return(0)
    } else {
      temp <- egpd_rng(count_draw, sigma, xi, kappa)
      return(sum(temp[is.finite(temp)]))
    }
  }
  burn_draws <- replicate(n, burns_rng(pi, delta, lambda, sigma, xi, kappa))
  return(mean(burn_draws, na.rm = TRUE))
}
  
#   for(i in 1:500) {
#     zero <- rbinom(1, 1, pi)
#     count_draw <- (1-zero) * rnbinom(1, size = delta, mu = exp(lambda))
#     if(count_draw == 0) {
#       burn_draws[i] <- 0
#     } else {
#       temp <- egpd_rng(count_draw, sigma, xi, kappa)
#       burn_draws[i] = sum(temp[is.finite(temp)])
#     }
#   }
#   return(mean(burn_draws, na.rm = TRUE))
# }

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

# extract sigma and xi from ONE CHAIN
g1_files <- paste0("full-model/fire-sims/burns/g1/csv-fits/",
                   list.files(path = "full-model/fire-sims/burns/g1/csv-fits", pattern = "sigma-ri_xi-ri_erc_12Sep"))
g1_files <- g1_files[grepl("long", g1_files)]
g1_files <- g1_files[grepl(paste0(chain, ".csv"), g1_files)]

print("Extracting sigma and xi...")
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
print("Sigma and xi extracted, moving on to kappa extraction...")
rm(rand_int_draws)
gc()

# extract kappa values from gen quant block (so it's all timepoints) of best G1 model
g1_gq_files <- paste0("full-model/fire-sims/model_comparison/extracted_values/",
                      list.files(path = "full-model/fire-sims/model_comparison/extracted_values/", 
                                 pattern = "g1_sigma-ri_xi-ri_erc_12"))
g1_gq_files <- g1_gq_files[grepl("long", g1_gq_files)]
g1_gq_files <- g1_gq_files[grepl(paste0(chain, ".csv"), g1_gq_files)]
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
print("Kappa has been extracted, moving on to ZINB model...")
rm(kappa_draws)
gc()

# extract count model params
zinb_files <- paste0("full-model/fire-sims/counts/zinb_er/csv-fits/",
                   list.files(path = "full-model/fire-sims/counts/zinb_er/csv-fits", pattern = "pi-ri_climate"))
zinb_files <- zinb_files[grepl("long", zinb_files)]
zinb_files <- zinb_files[grepl(paste0(chain, ".csv"), zinb_files)]

print("Extracting lambda train values...")
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
  pivot_wider(names_from = param, values_from = value) %>% 
  left_join(train_tmpts, join_by(time == model_tmpt)) %>% 
  select(-time) %>% 
  rename(time = true_tmpt)

print("Extracting lambda hold values...")
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
print("Full tibble of lambda values has been created, moving on to extraction of pi...")
rm(list = ls(pattern = "lambda_train"))
rm(list = ls(pattern = "lambda_hold"))
gc()

pi_prob_draws <- read_cmdstan_csv(zinb_files, variables = "pi_prob")$post_warmup_draws
pi_prob <- pi_prob_draws %>%
  as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>%
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  mutate(region = as.character(gsub("pi_prob\\[", "", name)),
         region = as.numeric(gsub("\\]", "", region)),
         name = "pi_prob") %>%
  pivot_wider(names_from = name, values_from = value)
print("Pi prob values have been extracted, moving on to extraction of delta...")
rm(pi_prob_draws)
gc()

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
print("Delta values have been extracted, moving on to joining all count and burn parameters...")
rm(delta_draws)
gc()

burn_preds <- kappa %>% left_join(rand_int) %>% left_join(lambda) %>% left_join(pi_prob) %>% left_join(delta)
print("Burn and count parameters have been joined for one chain, moving on to predictions. May take a while...")
burn_preds <- burn_preds %>%
  mutate(preds = purrr::pmap_dbl(list(pi_prob, delta, lambda, sigma, xi, kappa), burn_preds_gen, .progress = TRUE)) %>%
  select(c("draw", "time", "region", "preds"))
saveRDS(burn_preds, file = paste0("full-model/fire-sims/model_comparison/g1_burn_preds_", chain, ".RDS"))


