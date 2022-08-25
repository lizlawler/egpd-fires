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
toy_data <- stan_d
# toy_data <- read_rds("manuscript/scripts/toy_sim/g1/data/toy_data_ar1_icarphi_spatial_time_50regs.rds")
egpd_init <- stan_model('toy-sim/g1/stan/g1_icarphi_nu-only.stan')
egpd_fit <- sampling(egpd_init, 
                     data = toy_data, 
                     iter = 1000,
                     chains = 3,
                     refresh = 50)

# saveRDS(egpd_fit, file = "manuscript/scripts/toy_sim/g1/49reg_t1000.rds")

MCMCtrace(egpd_fit, params = c("rho1_kappa","rho1_nu", "rho1_xi", "rho2_kappa","rho2_nu", "rho2_xi",
                               "bp_kappa", "bp_nu", "bp_xi"), 
          ind = TRUE)

MCMCtrace(egpd_fit, params = c("beta_kappa", "beta_nu", "beta_xi", "phi_kappa", "phi_nu", "phi_xi"), 
          ind = TRUE)

# MCMCtrace(egpd_fit, params = c("beta_kappa", "bp","phi", "rho1", "rho2"), ind = TRUE,
#           gvals = c(betas_kappa, bp_kappa, phi_mat, rho1, rho2))
# 
# # save traceplot
# MCMCtrace(egpd_fit, params = c("beta_kappa", "bp", "phi", "rho1", "rho2"), 
#           ind = TRUE, 
#           gvals = c(betas_kappa, bp_kappa, phi_mat, rho1, rho2), 
#           open_pdf = FALSE, 
#           filename = paste0('manuscript/scripts/toy_sim/g1/49reg_t1000', 
#                             format(as.POSIXlt(Sys.time(), "America/Denver"), "%d%b%Y_%H%M"), ".pdf"))

# effects plot of comparison 
# need to regenerate truth first ------
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
## ------
post <- rstan::extract(egpd_fit, pars = c('beta_kappa', "beta_nu", "beta_xi"))
median_kappa <- apply(post$beta_kappa, c(2,3), median)
median_nu <- apply(post$beta_nu, c(2,3), median)
median_xi <- apply(post$beta_xi, c(2,3), median)

post_kappa_effects_df <- matrix(NA, r, t)
post_nu_effects_df <- matrix(NA, r, t)
post_xi_effects_df <- matrix(NA, r, t)
for(i in 1:r) {
  post_kappa_effects_df[i,] <- X_full[i, , ] %*% median_kappa[, i]
  post_nu_effects_df[i,] <- X_full[i, , ] %*% median_nu[, i]
  post_xi_effects_df[i,] <- X_full[i, , ] %*% median_xi[, i]
}

post_kappa <- t(post_kappa_effects_df) %>% as_tibble() %>% 
  rename_with(., ~ reg_cols) %>% 
  mutate(time = c(1:t)) %>%
  pivot_longer(cols = c(1:all_of(r)), values_to = "effect", names_to = "region") %>%
  left_join(., X_long) %>%
  left_join(., mod_reg_key) %>% mutate(type = "sim")
kappa_full <- rbind(kappa_effects, post_kappa) %>% mutate(type = factor(type, levels = c("truth", "sim")))
post_effects_kappa <- ggplot(kappa_full, aes(x=linear, y=effect, group = region)) + 
  geom_line(aes(linetype=NA_L1CODE, color = NA_L2CODE)) + 
  facet_grid(type ~ .)
post_effects_kappa

post_nu <- t(post_nu_effects_df) %>% as_tibble() %>% 
  rename_with(., ~ reg_cols) %>% 
  mutate(time = c(1:t)) %>%
  pivot_longer(cols = c(1:all_of(r)), values_to = "effect", names_to = "region") %>%
  left_join(., X_long) %>%
  left_join(., mod_reg_key) %>% mutate(type = "sim")
nu_full <- rbind(nu_effects, post_nu) %>% mutate(type = factor(type, levels = c("truth", "sim")))
post_effects_nu <- ggplot(nu_full, aes(x=linear, y=effect, group = region)) + 
  geom_line(aes(linetype=NA_L1CODE, color = NA_L2CODE)) + 
  facet_grid(type ~ .)
post_effects_nu

post_xi <- t(post_xi_effects_df) %>% as_tibble() %>% 
  rename_with(., ~ reg_cols) %>% 
  mutate(time = c(1:t)) %>%
  pivot_longer(cols = c(1:all_of(r)), values_to = "effect", names_to = "region") %>%
  left_join(., X_long) %>%
  left_join(., mod_reg_key) %>% mutate(type = "sim")
xi_full <- rbind(xi_effects, post_xi) %>% mutate(type = factor(type, levels = c("truth", "sim")))
post_effects_xi <- ggplot(xi_full, aes(x=linear, y=effect, group = region)) + 
  geom_line(aes(linetype=NA_L1CODE, color = NA_L2CODE)) + 
  facet_grid(type ~ .)
post_effects_xi

## maps of phi -----
ecoregions <- read_rds(file = "ecoregions.RDS")
ecoregions_geom <- ecoregions %>% filter(!NA_L2NAME == "UPPER GILA MOUNTAINS (?)")

pre_phi_kappa <- phi_mat_kappa %>% as_tibble() %>% 
  rename_with(., ~ reg_cols) %>%
  mutate(timepoint = 1:all_of(t)) %>%
  pivot_longer(cols = c(1:all_of(r)), values_to = "phi_val", names_to = "region") %>%
  left_join(., mod_reg_key) %>% mutate(type = "truth")

post_phi <- rstan::extract(egpd_fit, pars = c("phi_kappa", "phi_nu", "phi_xi"))
median_phi_kappa <- apply(post_phi$phi_kappa, c(2,3), median)
post_phi_kappa <- median_phi_kappa %>% as_tibble() %>%
  rename_with(., ~ reg_cols) %>%
  mutate(timepoint = 1:all_of(t)) %>%
  pivot_longer(cols = c(1:all_of(r)), values_to = "phi_val", names_to = "region") %>%
  left_join(., mod_reg_key) %>% mutate(type = "sim")

phi_full_kappa <- rbind(pre_phi_kappa, post_phi_kappa) %>% mutate(type = factor(type, levels = c("truth", "sim")))
joined_l3_phi_full_kappa <- right_join(ecoregions_geom, phi_full_kappa)
breaks <- classIntervals(c(min(joined_l3_phi_full_kappa$phi_val) - .00001, joined_l3_phi_full_kappa$phi_val), n=5, style = "quantile")
joined_l3_phi_full_kappa <- joined_l3_phi_full_kappa %>% mutate(phi_cat = cut(phi_val, unique(breaks$brks)))

full_phi_onetime_kappa <- ecoregions_geom %>%
  ggplot() +
  geom_sf(size = .1, fill = 'white') + 
  geom_sf(data = joined_l3_phi_full_kappa[joined_l3_phi_full_kappa$timepoint == 99, ], 
          aes(fill=phi_cat), alpha = 0.6, lwd = 0, inherit.aes = FALSE) + 
  facet_grid(type ~ .) +
  theme_minimal() + 
  theme(panel.grid.major = element_line(colour = "lightgrey"))
full_phi_onetime_kappa

pre_phi_nu <- phi_mat_nu %>% as_tibble() %>% 
  rename_with(., ~ reg_cols) %>%
  mutate(timepoint = 1:all_of(t)) %>%
  pivot_longer(cols = c(1:all_of(r)), values_to = "phi_val", names_to = "region") %>%
  left_join(., mod_reg_key) %>% mutate(type = "truth")

median_phi_nu <- apply(post_phi$phi_nu, c(2,3), median)
post_phi_nu <- median_phi_nu %>% as_tibble() %>%
  rename_with(., ~ reg_cols) %>%
  mutate(timepoint = 1:all_of(t)) %>%
  pivot_longer(cols = c(1:all_of(r)), values_to = "phi_val", names_to = "region") %>%
  left_join(., mod_reg_key) %>% mutate(type = "sim")

phi_full_nu <- rbind(pre_phi_nu, post_phi_nu) %>% mutate(type = factor(type, levels = c("truth", "sim")))
joined_l3_phi_full_nu <- right_join(ecoregions_geom, phi_full_nu)
breaks <- classIntervals(c(min(joined_l3_phi_full_nu$phi_val) - .00001, joined_l3_phi_full_nu$phi_val), n=5, style = "quantile")
joined_l3_phi_full_nu <- joined_l3_phi_full_nu %>% mutate(phi_cat = cut(phi_val, unique(breaks$brks)))

full_phi_onetime_nu <- ecoregions_geom %>%
  ggplot() +
  geom_sf(size = .1, fill = 'white') + 
  geom_sf(data = joined_l3_phi_full_nu[joined_l3_phi_full_nu$timepoint == 99, ], 
          aes(fill=phi_cat), alpha = 0.6, lwd = 0, inherit.aes = FALSE) + 
  facet_grid(type ~ .) +
  theme_minimal() + 
  theme(panel.grid.major = element_line(colour = "lightgrey"))
full_phi_onetime_nu

pre_phi_xi <- phi_mat_xi %>% as_tibble() %>% 
  rename_with(., ~ reg_cols) %>%
  mutate(timepoint = 1:all_of(t)) %>%
  pivot_longer(cols = c(1:all_of(r)), values_to = "phi_val", names_to = "region") %>%
  left_join(., mod_reg_key) %>% mutate(type = "truth")

median_phi_xi <- apply(post_phi$phi_xi, c(2,3), median)
post_phi_xi <- median_phi_xi %>% as_tibble() %>%
  rename_with(., ~ reg_cols) %>%
  mutate(timepoint = 1:all_of(t)) %>%
  pivot_longer(cols = c(1:all_of(r)), values_to = "phi_val", names_to = "region") %>%
  left_join(., mod_reg_key) %>% mutate(type = "sim")

phi_full_xi <- rbind(pre_phi_xi, post_phi_xi) %>% mutate(type = factor(type, levels = c("truth", "sim")))
joined_l3_phi_full_xi <- right_join(ecoregions_geom, phi_full_xi)
breaks <- classIntervals(c(min(joined_l3_phi_full_xi$phi_val) - .00001, joined_l3_phi_full_xi$phi_val), n=5, style = "quantile")
joined_l3_phi_full_xi <- joined_l3_phi_full_xi %>% mutate(phi_cat = cut(phi_val, unique(breaks$brks)))

full_phi_onetime_xi <- ecoregions_geom %>%
  ggplot() +
  geom_sf(size = .1, fill = 'white') + 
  geom_sf(data = joined_l3_phi_full_xi[joined_l3_phi_full_xi$timepoint == 99, ], 
          aes(fill=phi_cat), alpha = 0.6, lwd = 0, inherit.aes = FALSE) + 
  facet_grid(type ~ .) +
  theme_minimal() + 
  theme(panel.grid.major = element_line(colour = "lightgrey"))
full_phi_onetime_xi


# ggsave(paste0('~/Desktop/sim_plots/effects/g1/post_icarphi_spatial_time', 
#               format(as.POSIXlt(Sys.time(), "America/Denver"), "%d%b%Y_%H%M"),
#               ".pdf"), 
#        plot=post_effects, 
#        device = "pdf")

