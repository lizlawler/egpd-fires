library(readr)
library(rstan)
library(MCMCvis)
options(mc.cores = parallel::detectCores())

## baseline ------
toy_data <- read_rds('manuscript/scripts/toy_sim/g3/data/toy_data_g3_base.rds')
rstan_options(auto_write = TRUE)

egpd_init <- stan_model('manuscript/scripts/toy_sim/g3/stan/toy_egpd_g3_base.stan')
egpd_fit <- sampling(egpd_init,
                     data = toy_data,
                     iter = 1000,
                     chains = 3)

MCMCtrace(egpd_fit, params = c("beta_delta", "beta_nu", "beta_xi"), 
          ind = TRUE, 
          gvals = c(toy_data$truth$betas_delta, toy_data$truth$betas_nu, toy_data$truth$betas_xi))
# write_rds(egpd_fit_base, file = 'manuscript/scripts/toy_sim/egpd_fit_base.rds')
# ------

## small AR(1) ------
toy_data <- read_rds('manuscript/scripts/toy_sim/g3/data/toy_data_g3_ar1.rds')

egpd_init <- stan_model('manuscript/scripts/toy_sim/g3/stan/toy_egpd_g3_ar1.stan')
egpd_fit_g3 <- sampling(egpd_init,
                        data = toy_data,
                        iter = 1500,
                        chains = 3)
MCMCtrace(egpd_fit_g3, params = c("beta_delta", "beta_nu", "beta_xi", "rho1", "rho2", "bp_init"), 
          ind = TRUE, 
          gvals = c(toy_data$truth$betas_delta, toy_data$truth$betas_nu, toy_data$truth$betas_xi, 
                    toy_data$truth$rho1, toy_data$truth$rho2, toy_data$truth$bp*2))
# -------

# AR(1) with 10 regs ------
toy_data_g3 <- read_rds('manuscript/scripts/toy_sim/g3/data/toy_data_g3_corr_ar1_10reg.rds')

egpd_init_g3 <- stan_model('manuscript/scripts/toy_sim/g3/stan/toy_egpd_g3_ar1.stan')
egpd_fit_g3 <- sampling(egpd_init_g3,
                        data = toy_data_g3,
                        iter = 1000,
                        chains = 3)

MCMCtrace(egpd_fit_g3, params = c("beta_delta", "beta_nu", "beta_xi", "rho1", "rho2", "bp"), 
          ind = TRUE, 
          gvals = c(toy_data_g3$truth$betas_delta, toy_data_g3$truth$betas_nu, toy_data_g3$truth$betas_xi, 
                    toy_data_g3$truth$rho1, toy_data_g3$truth$rho2, toy_data_g3$truth$bp))

# diff priors ------
# egpd_init_g3_rhos <- stan_model('manuscript/scripts/toy_sim/g3/toy_egpd_g3_diffrhoprior.stan')
# egpd_fit_g3_rhos <- sampling(egpd_init_g3_rhos,
#                         data = toy_data_g3,
#                         iter = 1000,
#                         chains = 3)
# 
# egpd_init_g3_Zdelta <- stan_model('manuscript/scripts/toy_sim/g3/toy_egpd_g3_diffZdeltaprior.stan')
# egpd_fit_g3_Zdelta <- sampling(egpd_init_g3_Zdelta,
#                              data = toy_data_g3,
#                              iter = 1000,
#                              chains = 3)
# 
# MCMCtrace(egpd_fit_g3_Zdelta, params = c("beta_delta", "beta_nu", "beta_xi", "rho1", "rho2", "bp"), 
#           ind = TRUE, 
#           gvals = c(toy_data_g3$truth$betas_delta, toy_data_g3$truth$betas_nu, toy_data_g3$truth$betas_xi, 
#                     toy_data_g3$truth$rho1, toy_data_g3$truth$rho2, toy_data_g3$truth$bp))
### ----------

# create post effect plots
genSpline <- function(x, n, df = 5, degree, theta) {
  basis <- bs(x = x, df = df, degree = degree, 
              Boundary.knots = range(x), intercept = FALSE)
  full_design <- cbind(rep(1, n), x, basis)
  return(list(full_design, full_design %*% theta))
}

X_full <- toy_data_g3$X
X <- X_full[,2]
n <- toy_data_g3$n

# compare with truth
load(file = "region_key.RData")
library(tidyverse)
library(ggplot2)
library(splines)
mod_reg_key <- as_tibble(region_key) %>%
  mutate(region = sprintf("reg%d", 1:84),
         NA_L2CODE = as.factor(NA_L2CODE),
         NA_L1CODE = as.factor(NA_L1CODE)) %>%
  select(3:5)

post <- rstan::extract(egpd_fit_g3, pars = 'beta_delta')

median_delta <- apply(post$beta_delta, c(2,3), median)

delta_effects_df <- X_full %*% median_delta
post_delta <- delta_effects_df %>% as_tibble() %>%
  rename_with(., ~ gsub("V", "reg", .x, fixed = TRUE)) %>% 
  mutate(design = as.vector(X)) %>% 
  pivot_longer(cols = c(1:10), values_to = "effect", names_to = "region") %>%
  left_join(., mod_reg_key) %>% mutate(type = "sim")

# regenerate truth
df_delta <- X_full %*% toy_data_g3$truth$betas_delta
delta_effects <- df_delta %>% as_tibble() %>%
  rename_with(., ~ gsub("V", "reg", .x, fixed = TRUE)) %>%
  mutate(design = as.vector(X)) %>%
  # change cols to appropriate number of regions; if 15 regions, make 1:15, etc
  pivot_longer(cols = c(1:10), values_to = "effect", names_to = "region") %>%
  left_join(., mod_reg_key) %>%
  mutate(type = "truth")

delta_full <- rbind(delta_effects, post_delta) %>% mutate(type = factor(type, levels = c("truth", "sim")))

ggplot(delta_full, aes(x=design, y=effect, group = region)) + 
  geom_line(aes(linetype=NA_L1CODE, color = NA_L2CODE)) + facet_grid(type ~ .)

full_effects <- basis_df[[2]] %>% as_tibble() %>%
  rename_with(., ~ gsub("V", "reg", .x, fixed = TRUE)) %>% 
  mutate(design = as.vector(X)) %>% 
  # change cols to appropriate number of regions; if 15 regions, make 1:15, etc
  pivot_longer(cols = c(1:15), values_to = "effect", names_to = "region") %>%
  left_join(., mod_reg_key)

truth_plot <- ggplot(full_effects, aes(x=design, y=effect, group = region)) + 
  geom_line(aes(linetype=NA_L1CODE, color = NA_L2CODE))
truth_plot

# effect plot for original
post <- rstan::extract(egpd_fit_g3, pars = 'beta_delta')
median_all_regs <- apply(post$beta_delta, c(2,3), median)

post_effects_df <- genSpline(X, n, df = 5, degree = 3, median_all_regs)

post_effects <- post_effects_df[[2]] %>% as_tibble() %>%
  rename_with(., ~ gsub("V", "reg", .x, fixed = TRUE)) %>% 
  mutate(design = as.vector(X)) %>% 
  # change cols to appropriate number of regions; if 15 regions, make 1:15, etc
  pivot_longer(cols = c(1:10), values_to = "effect", names_to = "region") %>%
  left_join(., mod_reg_key)

post_plot_orig <- ggplot(post_effects, aes(x=design, y=effect, group = region)) + 
  geom_line(aes(linetype=NA_L1CODE, color = NA_L2CODE))
post_plot_orig

# effect plot for different rho priors
post <- rstan::extract(egpd_fit_g3_rhos, pars = 'beta_delta')
median_all_regs <- apply(post$beta_delta, c(2,3), median)

post_effects_df <- genSpline(X, n, df = 5, degree = 3, median_all_regs)

post_effects <- post_effects_df[[2]] %>% as_tibble() %>%
  rename_with(., ~ gsub("V", "reg", .x, fixed = TRUE)) %>% 
  mutate(design = as.vector(X)) %>% 
  # change cols to appropriate number of regions; if 15 regions, make 1:15, etc
  pivot_longer(cols = c(1:15), values_to = "effect", names_to = "region") %>%
  left_join(., mod_reg_key)

post_plot_rhos <- ggplot(post_effects, aes(x=design, y=effect, group = region)) + 
  geom_line(aes(linetype=NA_L1CODE, color = NA_L2CODE))
post_plot_rhos

# effect plot for different Zdelta prior
post <- rstan::extract(egpd_fit_g3_Zdelta, pars = 'beta_delta')
median_all_regs <- apply(post$beta_delta, c(2,3), median)

post_effects_df <- genSpline(X, n, df = 5, degree = 3, median_all_regs)

post_effects <- post_effects_df[[2]] %>% as_tibble() %>%
  rename_with(., ~ gsub("V", "reg", .x, fixed = TRUE)) %>% 
  mutate(design = as.vector(X)) %>% 
  # change cols to appropriate number of regions; if 15 regions, make 1:15, etc
  pivot_longer(cols = c(1:15), values_to = "effect", names_to = "region") %>%
  left_join(., mod_reg_key)

post_plot_Zdelta <- ggplot(post_effects, aes(x=design, y=effect, group = region)) + 
  geom_line(aes(linetype=NA_L1CODE, color = NA_L2CODE))
post_plot_Zdelta

