library(readr)
library(rstan)
library(MCMCvis)
library(tidyverse)
library(splines)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# base version -----
toy_data <- read_rds('manuscript/scripts/toy_sim/data/toy_data_simple.rds')
egpd_init <- stan_model('manuscript/scripts/toy_sim/toy_egpd_g1_base.stan')
egpd_fit <- sampling(egpd_init,
                    data = toy_data,
                    iter = 1000,
                    chains = 3)

MCMCtrace(egpd_fit, params = c("beta_kappa", "beta_nu", "beta_xi"), 
          ind = TRUE, 
          gvals = c(toy_data$truth$betas_kappa, toy_data$truth$betas_nu, toy_data$truth$betas_xi))

# write_rds(egpd_fit_base, file = 'manuscript/scripts/toy_sim/egpd_fit_base.rds')

# addition of correlation matrix ----
toy_data <- read_rds('manuscript/scripts/toy_sim/data/toy_data_corr.rds')

egpd_init <- stan_model('manuscript/scripts/toy_sim/toy_egpd_g1_corr.stan')
egpd_fit <- sampling(egpd_init,
                     data = toy_data,
                     iter = 1000,
                     chains = 3)

MCMCtrace(egpd_fit, params = c("beta_kappa", "beta_nu", "beta_xi", "rho1", "rho2"), 
          ind = TRUE, 
          gvals = c(toy_data$truth$betas_kappa, toy_data$truth$betas_nu, toy_data$truth$betas_xi, 0.1, 0.3))

# addition of AR(1) covariance matrix ----
toy_data <- read_rds('manuscript/scripts/toy_sim/g1/data/toy_data_ar1_xi.rds')

egpd_init_orig <- stan_model('manuscript/scripts/toy_sim/g1/stan/toy_egpd_g1_ar1.stan')
egpd_fit_orig <- sampling(egpd_init_orig,
                     data = toy_data,
                     iter = 1000,
                     chains = 3)

MCMCtrace(egpd_fit_orig, params = c("beta_kappa", "beta_nu", "beta_xi", "rho1", "rho2", "bp"), 
          ind = TRUE, 
          gvals = c(toy_data$truth$betas_kappa, toy_data$truth$betas_nu, toy_data$truth$betas_xi, 
                    toy_data$truth$rho1, toy_data$truth$rho2, toy_data$truth$bp))

# AR(1) and correlation, with all 3 regression parameters; same corr and cov matrices------
toy_data <- read_rds('manuscript/scripts/toy_sim/g1/data/toy_data_ar1_all3.rds')

egpd_init_orig <- stan_model('manuscript/scripts/toy_sim/g1/stan/g1_ar1_corr_all3.stan')
egpd_fit_orig <- sampling(egpd_init_orig,
                          data = toy_data,
                          iter = 1000,
                          chains = 3)

MCMCtrace(egpd_fit_orig, params = c("beta_kappa", "beta_nu", "beta_xi", "rho1", "rho2", "bp"), 
          ind = TRUE, 
          gvals = c(toy_data$truth$betas_kappa, toy_data$truth$betas_nu, toy_data$truth$betas_xi, 
                    toy_data$truth$rho1, toy_data$truth$rho2, toy_data$truth$bp))
# -----
# AR(1) and correlation, with all 3 regression parameters; same corr but 3 diff cov matrices------
toy_data <- read_rds('manuscript/scripts/toy_sim/g1/data/toy_data_diffar1_all3.rds')

egpd_init <- stan_model('manuscript/scripts/toy_sim/g1/stan/g1_diffar1_corr_all3.stan')
egpd_fit <- sampling(egpd_init, 
                     data = toy_data, 
                     iter = 1000,
                     chains = 3)

MCMCtrace(egpd_fit, params = c("beta_kappa", "beta_nu", "beta_xi", "rho1", "rho2", "bp_kappa", "bp_nu", "bp_xi"), 
          ind = TRUE, 
          gvals = c(toy_data$truth$betas_kappa, toy_data$truth$betas_nu, toy_data$truth$betas_xi, 
                    toy_data$truth$rho1, toy_data$truth$rho2, toy_data$truth$bp_kappa, toy_data$truth$bp_nu, toy_data$truth$bp_xi))
# -------

# AR(1) and correlation, with all 3 regression parameters; 3 diff corr and 3 diff cov matrices------
toy_data <- read_rds('manuscript/scripts/toy_sim/g1/data/toy_data_diffar1_all3.rds')

egpd_init <- stan_model('manuscript/scripts/toy_sim/g1/stan/g1_diffar1_diffcorr_all3.stan')
egpd_fit <- sampling(egpd_init, 
                     data = toy_data, 
                     iter = 1000,
                     chains = 3)

MCMCtrace(egpd_fit, params = c("beta_kappa", "beta_nu", "beta_xi", 
                               "rho1_kappa", "rho2_kappa", "rho1_nu", "rho2_nu",  "rho1_xi", "rho2_xi",
                               "bp_kappa", "bp_nu", "bp_xi"), 
          ind = TRUE, 
          gvals = c(toy_data$truth$betas_kappa, toy_data$truth$betas_nu, toy_data$truth$betas_xi, 
                    toy_data$truth$rho1, toy_data$truth$rho2, toy_data$truth$rho1, toy_data$truth$rho2, toy_data$truth$rho1, toy_data$truth$rho2, 
                    toy_data$truth$bp_kappa, toy_data$truth$bp_nu, toy_data$truth$bp_xi))

# AR(1) and correlation, with regression on kappa only; added in ICAR for spatial effect (haven't added in time AR piece yet) ------
toy_data <- read_rds('manuscript/scripts/toy_sim/g1/data/toy_data_ar1_icarphi_spatial.rds')

egpd_init <- stan_model('manuscript/scripts/toy_sim/g1/stan/g1_ar1_icarphi_kappa.stan')
egpd_fit <- sampling(egpd_init, 
                     data = toy_data, 
                     iter = 1000,
                     chains = 3)
egpd_fit_vb <- vb(egpd_init, data = toy_data)

MCMCtrace(egpd_fit, params = c("beta_kappa", "beta_nu", "beta_xi", 
                               "rho1", "rho2", "bp", "phi"), 
          ind = TRUE, 
          gvals = c(toy_data$truth$betas_kappa, toy_data$truth$betas_nu, toy_data$truth$betas_xi, 
                    toy_data$truth$rho1, toy_data$truth$rho2, 
                    toy_data$truth$bp_kappa, toy_data$truth$phi),
          open_pdf = FALSE,
          filename = "test.pdf")


# make effects plots
X_full <- toy_data$X
X <- X_full[, 2]
load(file = "region_key.RData")
mod_reg_key <- as_tibble(region_key) %>% 
  mutate(region = sprintf("reg%d", 1:84),
         NA_L2CODE = as.factor(NA_L2CODE),
         NA_L1CODE = as.factor(NA_L1CODE)) %>% 
  select(3:5)

# effect plots
post <- rstan::extract(egpd_fit, pars = 'beta_kappa')

# kappa
# post effects
median_kappa <- apply(post$beta_kappa, c(2,3), median)
kappa_effects_df <- X_full %*% median_kappa
post_kappa <- kappa_effects_df %>% as_tibble() %>%
  rename_with(., ~ gsub("V", "reg", .x, fixed = TRUE)) %>% 
  mutate(design = as.vector(X)) %>% 
  pivot_longer(cols = c(1:15), values_to = "effect", names_to = "region") %>%
  left_join(., mod_reg_key) %>% mutate(type = "sim")

# regenerate truth
df_kappa <- X_full %*% toy_data$truth$betas_kappa
kappa_effects <- df_kappa %>% as_tibble() %>%
  rename_with(., ~ gsub("V", "reg", .x, fixed = TRUE)) %>%
  mutate(design = as.vector(X)) %>%
  # change cols to appropriate number of regions; if 15 regions, make 1:15, etc
  pivot_longer(cols = c(1:15), values_to = "effect", names_to = "region") %>%
  left_join(., mod_reg_key) %>%
  mutate(type = "truth")

kappa_full <- rbind(kappa_effects, post_kappa) %>% mutate(type = factor(type, levels = c("truth", "sim")))

ggplot(kappa_full, aes(x=design, y=effect, group = region)) + 
  geom_line(aes(linetype=NA_L1CODE, color = NA_L2CODE)) + facet_grid(type ~ .)


# nu
# post effects
median_nu <- apply(post$beta_nu, c(2,3), median)
nu_effects_df <- X_full %*% median_nu
post_nu <- nu_effects_df %>% as_tibble() %>%
  rename_with(., ~ gsub("V", "reg", .x, fixed = TRUE)) %>% 
  mutate(design = as.vector(X)) %>% 
  pivot_longer(cols = c(1:15), values_to = "effect", names_to = "region") %>%
  left_join(., mod_reg_key) %>% mutate(type = "sim")

# regenerate truth
df_nu <- X_full %*% toy_data$truth$betas_nu
nu_effects <- df_nu %>% as_tibble() %>%
  rename_with(., ~ gsub("V", "reg", .x, fixed = TRUE)) %>%
  mutate(design = as.vector(X)) %>%
  # change cols to appropriate number of regions; if 15 regions, make 1:15, etc
  pivot_longer(cols = c(1:15), values_to = "effect", names_to = "region") %>%
  left_join(., mod_reg_key) %>%
  mutate(type = "truth")

nu_full <- rbind(nu_effects, post_nu) %>% mutate(type = factor(type, levels = c("truth", "sim")))

ggplot(nu_full, aes(x=design, y=effect, group = region)) + 
  geom_line(aes(linetype=NA_L1CODE, color = NA_L2CODE)) + facet_grid(type ~ .)

# xi
# post effects
median_xi <- apply(post$beta_xi, c(2,3), median)
xi_effects_df <- X_full %*% median_xi
post_xi <- xi_effects_df %>% as_tibble() %>%
  rename_with(., ~ gsub("V", "reg", .x, fixed = TRUE)) %>% 
  mutate(design = as.vector(X)) %>% 
  pivot_longer(cols = c(1:15), values_to = "effect", names_to = "region") %>%
  left_join(., mod_reg_key) %>% mutate(type = "sim")

# regenerate truth
df_xi <- X_full %*% toy_data$truth$betas_xi
xi_effects <- df_xi %>% as_tibble() %>%
  rename_with(., ~ gsub("V", "reg", .x, fixed = TRUE)) %>%
  mutate(design = as.vector(X)) %>%
  # change cols to appropriate ximber of regions; if 15 regions, make 1:15, etc
  pivot_longer(cols = c(1:15), values_to = "effect", names_to = "region") %>%
  left_join(., mod_reg_key) %>%
  mutate(type = "truth")

xi_full <- rbind(xi_effects, post_xi) %>% mutate(type = factor(type, levels = c("truth", "sim")))

ggplot(xi_full, aes(x=design, y=effect, group = region)) + 
  geom_line(aes(linetype=NA_L1CODE, color = NA_L2CODE)) + facet_grid(type ~ .)

