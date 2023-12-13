args <- commandArgs(trailingOnly=TRUE)
type <- args[1]
model <- args[2]
params <- args[3]
dataset <- args[4]
sttime <- args[5]

library(cmdstanr)
set_cmdstan_path(path = "/projects/eslawler@colostate.edu/software/anaconda/envs/lawler/bin/cmdstan") # this is only relevant to Alpine
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(MCMCvis)
library(posterior)
library(tidyr)
library(dplyr)

csvbase <- paste0("./full-model/fire-sims/", type, "/csv-fits/")
csvpattern <- paste0(type, "_", model, "_", params, "_", dataset, "_", sttime)
files <- paste0(csvbase, list.files(path = csvbase, pattern = csvpattern))
print("Filenames being used are:")
files

# create directory to save mcmc draws 
draws_path <- paste0("full-model/figures/paper/mcmc_draws/joint_", model, "_", params, "_", dataset, "/")
dir.create(path = draws_path, recursive = TRUE)

burn_pred_draws <- read_cmdstan_csv(files, variables = "burn_pred")$post_warmup_draws
burn_pred <- burn_pred_draws %>%
  as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>%
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  separate_wider_delim(cols = "name", delim = ",", names = c("time", "region")) %>%
  mutate(time = as.numeric(gsub("burn_pred\\[", "", time)),
         region = as.numeric(gsub("\\]", "", region)))
saveRDS(burn_pred, file = paste0(draws_path, "burn_pred", ".RDS"))
rm(list=ls(pattern = "burn_pred"))
gc()

kappa_draws <- read_cmdstan_csv(files, variables = "reg_full")$post_warmup_draws
kappa <- kappa_draws %>%
  as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>%
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  separate_wider_delim(cols = "name", delim = ",", names = c("time", "region")) %>%
  mutate(time = as.numeric(gsub("reg_full\\[", "", time)),
         region = as.numeric(gsub("\\]", "", region)),
         kappa = exp(value)) %>% select(-value)
saveRDS(kappa, file = paste0(draws_path, "kappa", ".RDS"))
rm(list=ls(pattern = "kappa"))
gc()

# inv_logit <- function(value) exp(value) / (exp(value) + 1)
rand_int_draws <- read_cmdstan_csv(files, variables = "ri_init")$post_warmup_draws
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
         # value = case_when(param == "sigma" ~ exp(value),
         #                   param == "xi" ~ inv_logit(value) * 1.5),
         value = exp(value)) %>%
  pivot_wider(names_from = param, values_from = value)
saveRDS(rand_int, file = paste0(draws_path, "rand_int", ".RDS"))
rm(list=ls(pattern = "rand_int"))
gc()

lambda_draws <- read_cmdstan_csv(files, variables = "lambda_full")$post_warmup_draws
lambda <- lambda_draws %>%
  as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>%
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  separate_wider_delim(cols = "name", delim = ",", names = c("time", "region")) %>%
  mutate(time = as.numeric(gsub("lambda_full\\[", "", time)),
         region = as.numeric(gsub("\\]", "", region)),
         param = "lambda") %>%
  pivot_wider(names_from = param, values_from = value)
saveRDS(lambda, file = paste0(draws_path, "lambda", ".RDS"))
rm(list=ls(pattern = "lambda"))
gc()

phi_draws <- read_cmdstan_csv(files, variables = "phi")$post_warmup_draws
phi <- phi_draws %>%
  as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>%
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  separate_wider_delim(cols = "name", delim = ",", names = c("param", "time", "region")) %>%
  mutate(param = as.character(gsub("phi\\[", "", param)),
         time = as.numeric(time),
         region = as.numeric(gsub("\\]", "", region)),
         param = case_when(param == "1" ~ "lambda",
                           param == "2" ~ "kappa",
                           TRUE ~ param)) %>%
  pivot_wider(names_from = param, values_from = value)
saveRDS(phi, file = paste0(draws_path, "phi", ".RDS"))
rm(list=ls(pattern = "phi"))
gc()

pi_prob_draws <- read_cmdstan_csv(files, variables = "pi_prob")$post_warmup_draws
pi_prob <- pi_prob_draws %>%
  as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>%
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  mutate(region = as.character(gsub("pi_prob\\[", "", name)),
         region = as.numeric(gsub("\\]", "", region)),
         name = "pi") %>%
  pivot_wider(names_from = name, values_from = value)
saveRDS(pi_prob, file = paste0(draws_path, "pi_prob", ".RDS"))
rm(list=ls(pattern = "pi_prob"))
gc()

gamma_draws <- read_cmdstan_csv(files, variables = "gamma")$post_warmup_draws
gamma <- gamma_draws %>%
  as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>%
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  mutate(region = as.character(gsub("gamma\\[", "", name)),
         region = as.numeric(gsub("\\]", "", region)),
         name = "gamma") %>%
  pivot_wider(names_from = name, values_from = value)
saveRDS(gamma, file = paste0(draws_path, "gamma", ".RDS"))
rm(list=ls(pattern = "gamma"))
gc()

theta_draws <- read_cmdstan_csv(files, variables = "theta")$post_warmup_draws
theta <- theta_draws %>% as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>%
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  mutate(timepoint = as.character(gsub("theta\\[", "", name)),
         timepoint = as.numeric(gsub("\\]", "", timepoint)),
         name = "theta") %>%
  pivot_wider(names_from = name, values_from = value)
saveRDS(theta, file = paste0(draws_path, "theta", ".RDS"))
rm(list=ls(pattern = "theta"))
gc()

delta_draws <- read_cmdstan_csv(files, variables = "delta")$post_warmup_draws
delta <- delta_draws %>% as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>%
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  mutate(name = gsub("delta\\[", "", name),
         name = as.integer(gsub("\\]", "", name))) %>%
  rename(region = name, delta = value)
saveRDS(delta, file = paste0(draws_path, "delta", ".RDS"))
rm(list=ls(pattern = "delta"))
gc()

beta_count_draws <- read_cmdstan_csv(files, variables = "beta_count")$post_warmup_draws
beta_count <- beta_count_draws %>% as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>%
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw")
saveRDS(beta_count, file = paste0(draws_path, "beta_count", ".RDS"))
rm(list=ls(pattern = "beta_count"))
gc()

beta_burn_draws <- read_cmdstan_csv(files, variables = "beta_burn")$post_warmup_draws
beta_burn <- beta_burn_draws %>% as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>%
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw")
saveRDS(beta_burn, file = paste0(draws_path, "beta_burn", ".RDS"))
rm(list=ls(pattern = "beta_burn"))
gc()

tau_draws <- read_cmdstan_csv(files, variables = "tau")$post_warmup_draws
tau <- tau_draws %>% as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>%
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw")
saveRDS(tau, file = paste0(draws_path, "tau", ".RDS"))
rm(list=ls(pattern = "tau"))
gc()

eta_draws <- read_cmdstan_csv(files, variables = "eta")$post_warmup_draws
eta <- eta_draws %>% as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>%
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw")
saveRDS(eta, file = paste0(draws_path, "eta", ".RDS"))
rm(list=ls(pattern = "eta"))
gc()

bp_draws <- read_cmdstan_csv(files, variables = "bp")$post_warmup_draws
bp <- bp_draws %>% as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>%
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw")
saveRDS(bp, file = paste0(draws_path, "bp", ".RDS"))
rm(list=ls(pattern = "bp"))
gc()

rho1_draws <- read_cmdstan_csv(files, variables = "rho1")$post_warmup_draws
rho1 <- rho1_draws %>% as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>%
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw")
saveRDS(rho1, file = paste0(draws_path, "rho1", ".RDS"))
rm(list=ls(pattern = "rho1"))
gc()

rho2_draws <- read_cmdstan_csv(files, variables = "rho2")$post_warmup_draws
rho2 <- rho2_draws %>% as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>%
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw")
saveRDS(rho2, file = paste0(draws_path, "rho2", ".RDS"))

print("All MCMC draws have been successfully written to disk")
