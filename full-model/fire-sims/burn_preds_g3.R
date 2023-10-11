args <- commandArgs(trailingOnly=TRUE)
chain <- args[1]
ineq <- args[2]

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

# egpd functions - G3, incorporating truncation directly
egpd_ccdf <- function(y, sigma, xi, delta_burn) {
  alpha <- 1/delta_burn
  beta <- 2
  inner <- exp(delta_burn * pgpd(y, mu = 0, sigma = sigma, xi = xi, lower.tail = FALSE, log.p = TRUE))
  return(exp(pbeta(inner, shape1 = alpha, shape2 = beta, log.p = TRUE)))
}

egpd_icdf <- function(u, sigma, xi, delta_burn) {
  alpha <- 1/delta_burn
  beta <- 2
  inner <- 1 - qbeta(u, shape1 = alpha, shape2 = beta)^alpha
  return(qgpd(inner, mu = 0, sigma = sigma, xi = xi))
}

egpd_rng <- function(n, sigma, xi, delta_burn) {
  cst <- egpd_ccdf(1.001, sigma, xi, delta_burn)
  u <- runif(n)
  rng_val <- rep(NA, n)
  for(i in 1:n) {
    u_adj = cst * (1 - u[i])
    rng_val[i] = egpd_icdf(u_adj, sigma, xi, delta_burn)
  }
  return(rng_val)
}

burn_preds_gen <- function(n = 500, pi_prob, delta_count, lambda, sigma, xi, delta_burn) {
  pi <- exp(pi_prob)/(1+exp(pi_prob))
  burns_rng <- function(pi, delta_count, lambda, sigma, xi, delta_burn) {
    zero <- rbinom(1, 1, pi)
    count_draw <- (1-zero) * rnbinom(1, size = delta_count, mu = exp(lambda))
    if(count_draw == 0) {
      return(0)
    } else {
      temp <- egpd_rng(count_draw, sigma, xi, delta_burn)
      return(sum(temp[is.finite(temp)]))
    }
  }
  burn_draws <- replicate(n, burns_rng(pi, delta_count, lambda, sigma, xi, delta_burn))
  return(mean(burn_draws, na.rm = TRUE))
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

# extract xi and delta_burn from ONE CHAIN
g3_files <- paste0("full-model/fire-sims/burns/g3/csv-fits/",
                   list.files(path = "full-model/fire-sims/burns/g3/csv-fits", pattern = "erc_xi-ri_nu"))
g3_files <- g3_files[grepl(paste0(chain, ".csv"), g3_files)]

print("Extracting xi and delta...")
rand_int_draws <- read_cmdstan_csv(g3_files, variables = "ri_init")$post_warmup_draws
rand_int <- rand_int_draws %>%
  as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>%
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  separate_wider_delim(cols = "name", delim = ",", names = c("param", "region")) %>%
  mutate(param = as.character(gsub("ri_init\\[", "", param)),
         region = as.numeric(gsub("\\]", "", region)),
         param = as.character(param),
         param = case_when(param == "1" ~ "xi",
                           param == "2" ~ "delta_burn",
                           TRUE ~ param),
         value = exp(value)) %>%
  pivot_wider(names_from = param, values_from = value)
print("Xi and delta extracted, moving on to nu extraction...")
rm(rand_int_draws)
gc()

# extract nu values of best g3 model
nu_draws <- read_cmdstan_csv(g3_files, variables = "reg_full")$generated_quantities
nu <- nu_draws %>%
  as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>%
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  separate_wider_delim(cols = "name", delim = ",", names = c("param", "time", "region")) %>%
  mutate(param = as.numeric(gsub("reg_full\\[", "", param)),
         time = as.numeric(time),
         region = as.numeric(gsub("\\]", "", region)),
         nu = exp(value)) %>% select(-value)
nu <- nu %>% select(-param)
print("Nu has been extracted, moving on to ZINB model...")
rm(nu_draws)
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
delta_count <- delta_draws %>%
  as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>%
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  mutate(region = as.character(gsub("delta\\[", "", name)),
         region = as.numeric(gsub("\\]", "", region)),
         name = "delta_count") %>%
  pivot_wider(names_from = name, values_from = value)
print("Delta values have been extracted, moving on to joining all count and burn parameters...")
rm(delta_draws)
gc()

burn_preds <- kappa %>% 
  left_join(rand_int) %>% 
  left_join(lambda) %>% 
  left_join(pi_prob) %>% 
  left_join(delta_count) %>% 
  mutate(sigma = nu / (1+xi)) %>% 
  select(-nu)

# attempts to speed up this process by cutting it in half and distributing on more nodes
if(ineq == "less") {
  burn_preds <- burn_preds %>% filter(draw <= 500)
} else {
  burn_preds <- burn_preds %>% filter(draw > 500)
}

print("Burn and count parameters have been joined for one chain, moving on to predictions. May take a while...")
burn_preds <- burn_preds %>%
  mutate(preds = purrr::pmap_dbl(list(500, pi_prob, delta_count, lambda, sigma, xi, delta_burn), 
                                 burn_preds_gen, 
                                 .progress = TRUE)) %>%
  select(c("draw", "time", "region", "preds"))
saveRDS(burn_preds, file = paste0("full-model/fire-sims/model_comparison/g3_burn_preds_", chain, ".RDS"))


