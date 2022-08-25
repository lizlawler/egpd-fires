library(readr)
library(rstan)
library(MCMCvis)
library(tidyverse)
library(splines)
library(spdep)
library(sf)
library(RColorBrewer)
library(patchwork)
library(classInt)
library(spatialreg)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# generate toy data
# source("manuscript/scripts/toy_sim/g1/rep_summit/g1_icar_full_data.R")

# run sampling
# toy_data <- toy_data_icarphi_st
toy_data <- stan_d
# toy_data <- read_rds("manuscript/scripts/toy_sim/zip/data/toy_data_zip_lambdaonly_betasonly.rds")
egpd_init <- stan_model('./toy-sim/zip/stan/zip_icar.stan')
egpd_fit <- sampling(egpd_init, 
                     data = toy_data, 
                     iter = 1000,
                     chains = 3,
                     refresh = 50)

saveRDS(egpd_fit, file = "manuscript/scripts/toy_sim/zip/zip_49reg_t1000.rds")

# MCMCtrace(egpd_fit, params = c("beta_lambda", "bp", "rho1", "rho2"), 
#           ind = TRUE,
#           gvals = c(betas_lambda, bp_lambda, rho1, rho2))
MCMCtrace(egpd_fit, params = c("beta_lambda", "bp", "rho1", "rho2", "phi"), 
          ind = TRUE)

# save traceplot
MCMCtrace(egpd_fit, params = c("beta_kappa", "beta_nu", "beta_xi", "bp", 
                               "eta", "phi_post", "rho1", "rho2", "tau"), 
          ind = TRUE, 
          gvals = c(betas_kappa, betas_nu, betas_xi, bp_kappa, eta, phi_mat, rho1, rho2, tau), 
          open_pdf = FALSE, 
          filename = paste0('manuscript/scripts/toy_sim/g1/rep_summit/trace_g1_full_icar_', 
                            format(as.POSIXlt(Sys.time(), "America/Denver"), "%d%b%Y_%H%M"), ".pdf"))



# save effects plot of comparison 
# need to regenerate truth first
# r <- toy_data$r
# t <- toy_data$t
# p <- toy_data$p
# 
# load(file = "region_key.RData")
# mod_reg_key <- as_tibble(region_key) %>% 
#   mutate(region = sprintf("reg%d", 1:84),
#          NA_L2CODE = as.factor(NA_L2CODE),
#          NA_L1CODE = as.factor(NA_L1CODE)) %>% 
#   select(3:5)
# 
# X_full <- toy_data$X
# betas_kappa <- toy_data$truth$betas_kappa
# df_kappa <- matrix(NA, r, t)
# for(i in 1:r) {
#   df_kappa[i,] <- X_full[i, , ] %*% betas_kappa[, i]
# }
# 
# X <- t(X_full[1:15, ,2])
# 
# X_long <- X %>% as_tibble() %>% mutate(time = c(1:t)) %>%
#   rename_with(., ~ gsub("V", "reg", .x, fixed = TRUE)) %>%
#   pivot_longer(cols = c(1:all_of(r)), values_to = "linear", names_to = "region")
# 
# kappa_effects <- t(df_kappa) %>% as_tibble() %>% mutate(time = c(1:t)) %>%
#   rename_with(., ~ gsub("V", "reg", .x, fixed = TRUE)) %>%
#   pivot_longer(cols = c(1:all_of(r)), values_to = "effect", names_to = "region") %>%
#   left_join(., X_long) %>%
#   left_join(., mod_reg_key) %>% mutate(type = "truth")


post <- rstan::extract(egpd_fit, pars = 'beta_lambda')
median_lambda <- apply(post$beta_lambda, c(2,3), median)

post_lambda_effects_df <- matrix(NA, r, t)
for(i in 1:r) {
  post_lambda_effects_df[i,] <- X_full[i, , ] %*% median_lambda[, i]
}

post_lambda <- t(post_lambda_effects_df) %>% as_tibble() %>% 
  rename_with(., ~ reg_cols) %>% 
  mutate(time = c(1:t)) %>%
  pivot_longer(cols = c(1:all_of(r)), values_to = "effect", names_to = "region") %>%
  left_join(., X_long) %>%
  left_join(., mod_reg_key) %>% mutate(type = "sim")

lambda_full <- rbind(lambda_effects, post_lambda) %>% mutate(type = factor(type, levels = c("truth", "sim")))

post_effects <- ggplot(lambda_full, aes(x=linear, y=effect, group = region)) + 
  geom_line(aes(linetype=NA_L1CODE, color = NA_L2CODE)) + 
  facet_grid(type ~ .)
post_effects
ggsave(paste0('~/Desktop/sim_plots/effects/zip/icar_lambda_betasonly', 
              format(as.POSIXlt(Sys.time(), "America/Denver"), "%d%b%Y_%H%M"),
              ".pdf"), 
       plot=post_effects, 
       device = "pdf")


## maps of phi -------
nb <- read_rds('~/Desktop/csu/research/josephs_paper/data/processed/nb.rds')
ecoregions <- read_rds(file = "ecoregions.RDS")

ecoregions_geom <- ecoregions %>% filter(!NA_L2NAME == "UPPER GILA MOUNTAINS (?)")

### maps of all regions labeled ------
# ecoregions_geom %>%
#   ggplot() +
#   geom_sf(size = .1, fill = 'white') +
#   geom_sf(data = ecoregions_geom,
#           aes(fill=NA_L1CODE), alpha = 0.6, lwd = 0, inherit.aes = FALSE) +
#   geom_sf_label(aes(label = NA_L2CODE)) +
#   theme_minimal() +
#   theme(panel.grid.major = element_line(colour = "lightgrey"))
# -------
pre_phi_times <- phi_mat %>% as_tibble() %>% 
  rename_with(., ~ reg_cols) %>%
  mutate(timepoint = 1:all_of(t)) %>%
  pivot_longer(cols = c(1:all_of(r)), values_to = "phi_val", names_to = "region") %>%
  left_join(., mod_reg_key) %>% mutate(type = "truth")

post_phi <- rstan::extract(egpd_fit, pars = "phi")
median_phi <- apply(post_phi$phi, c(2,3), median)
post_phi_times <- median_phi %>% as_tibble() %>%
  rename_with(., ~ reg_cols) %>%
  mutate(timepoint = 1:all_of(t)) %>%
  pivot_longer(cols = c(1:all_of(r)), values_to = "phi_val", names_to = "region") %>%
  left_join(., mod_reg_key) %>% mutate(type = "sim")

phi_full <- rbind(pre_phi_times, post_phi_times) %>% mutate(type = factor(type, levels = c("truth", "sim")))
joined_l3_phi_full <- right_join(ecoregions_geom, phi_full)
breaks <- classIntervals(c(min(joined_l3_phi_full$phi_val) - .00001, joined_l3_phi_full$phi_val), n=5, style = "quantile")
joined_l3_phi_full <- joined_l3_phi_full %>% mutate(phi_cat = cut(phi_val, unique(breaks$brks)))

full_phi_l3_onetime <- ecoregions_geom %>%
  ggplot() +
  geom_sf(size = .1, fill = 'white') + 
  geom_sf(data = joined_l3_phi_full[joined_l3_phi_full$timepoint == 99, ], 
          aes(fill=phi_cat), alpha = 0.6, lwd = 0, inherit.aes = FALSE) + 
  facet_grid(type ~ .) +
  theme_minimal() + 
  theme(panel.grid.major = element_line(colour = "lightgrey"))
full_phi_l3_onetime

