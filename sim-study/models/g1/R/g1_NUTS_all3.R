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
options(mc.cores = parallel::detectCores())

# generate toy data
source("./sim-study/models/g1/R/g1_data_all3.R")

# run sampling
egpd_init <- stan_model('./sim-study/models/g1/stan/g1_all3.stan')
egpd_fit <- sampling(egpd_init, 
                     data = toy_data, 
                     iter = 1000,
                     chains = 3,
                     refresh = 50)

MCMCtrace(egpd_fit, params = c("beta_kappa", "beta_nu", "beta_xi",
                               "phi_kappa", "phi_nu", "phi_xi",
                               "rho1_kappa", "rho2_kappa", "rho1_nu", "rho2_nu", "rho1_xi", "rho2_xi",
                               "bp_kappa", "bp_nu", "bp_xi"),
          ind = TRUE,
          open_pdf = FALSE,
          filename = paste0('./sim-study/figures/g1/trace/g1_trace_all3_matnorm',
                            format(as.POSIXlt(Sys.time(), "America/Denver"), "%d%b%Y_%H%M"), ".pdf"))

# save MCMC object in case below dx plots don't save properly
saveRDS(egpd_fit, file = "./sim-study/models/g1/stan-fits/g1_all3_with-matnorm.RDS")

# effects plot of comparison --------
post <- rstan::extract(egpd_fit, pars = c("beta_kappa", "beta_nu", "beta_xi", "phi_kappa", "phi_nu", "phi_xi"))
median_kappa <- apply(post$beta_kappa, c(2,3), median)
median_nu <- apply(post$beta_nu, c(2,3), median)
median_xi <- apply(post$beta_nu, c(2,3), median)

post_kappa_effects_df <- matrix(NA, r, t)
for(i in 1:r) {
  post_kappa_effects_df[i,] <- X_full[i, , ] %*% median_kappa[, i]
}

post_kappa <- t(post_kappa_effects_df) %>% as_tibble() %>% 
  rename_with(., ~ reg_cols) %>% 
  mutate(time = c(1:t)) %>%
  pivot_longer(cols = c(1:all_of(r)), values_to = "effect", names_to = "region") %>%
  left_join(., X_long) %>%
  left_join(., mod_reg_key) %>% mutate(type = "sim")
kappa_full <- rbind(kappa_effects, post_kappa) %>% 
  mutate(type = case_when(type == 'truth' ~ 'Truth',
                          type == 'sim' ~ 'Simulated'),
         type = factor(type, levels = c("Truth", "Simulated")))
post_effects_kappa <- ggplot(kappa_full, aes(x=linear, y=effect, group = region)) + 
  geom_line(aes(linetype=NA_L1CODE, color = NA_L2CODE)) + 
  facet_grid(. ~ type) + xlab("Linear term") + ylab(expression("Partial effect on"~kappa)) +
  guides(color = guide_legend("Level 2"),
         linetype = guide_legend("Level 1")) +
  theme_classic()
ggsave("effects_plot_classic_v3.png", dpi = 700, type = "cairo")
ggsave(paste0('./sim-study/figures/g1/effects/g1_post-effects_kappa-nu-sim_kappa_',
              format(as.POSIXlt(Sys.time(), "America/Denver"), "%d%b%Y_%H%M"),
              ".pdf"),
       plot=post_effects_kappa,
       device = "pdf")

post_nu_effects_df <- matrix(NA, r, t)
for(i in 1:r) {
  post_nu_effects_df[i,] <- X_full[i, , ] %*% median_nu[,i]
}

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
ggsave(paste0('./sim-study/figures/g1/effects/g1_post-effects_kappa-nu-sim_nu_',
              format(as.POSIXlt(Sys.time(), "America/Denver"), "%d%b%Y_%H%M"),
              ".pdf"),
       plot=post_effects_nu,
       device = "pdf")

post_xi_effects_df <- matrix(NA, r, t)
for(i in 1:r) {
  post_xi_effects_df[i,] <- X_full[i, , ] %*% median_xi[i, ]
}

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
ggsave(paste0('./sim-study/figures/g1/effects/g1_post-effects_kappa-nu-sim_xi_',
              format(as.POSIXlt(Sys.time(), "America/Denver"), "%d%b%Y_%H%M"),
              ".pdf"),
       plot=post_effects_xi,
       device = "pdf")

## maps of phi -----
# kappa maps
ecoregions_geom <- ecoregions %>% filter(!NA_L2NAME == "UPPER GILA MOUNTAINS (?)")
pre_phi_kappa <- phi_mat_kappa %>% as_tibble() %>% 
  rename_with(., ~ reg_cols) %>%
  mutate(timepoint = 1:all_of(t)) %>%
  pivot_longer(cols = c(1:all_of(r)), values_to = "phi_val", names_to = "region") %>%
  left_join(., mod_reg_key) %>% mutate(type = "truth")

median_phi_kappa <- apply(post$phi_kappa, c(2,3), median)
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

save_kappa <- full_phi_onetime_kappa + scale_fill_brewer(palette = "YlOrRd")

ggsave(paste0('./sim-study/figures/g1/maps/g1_phi-map_kappa-nu-sim_kappa_',
              format(as.POSIXlt(Sys.time(), "America/Denver"), "%d%b%Y_%H%M"),
              ".pdf"),
       plot=save_kappa,
       device = "pdf")

# nu maps -----
pre_phi_nu <- phi_mat_nu %>% as_tibble() %>% 
  rename_with(., ~ reg_cols) %>%
  mutate(timepoint = 1:all_of(t)) %>%
  pivot_longer(cols = c(1:all_of(r)), values_to = "phi_val", names_to = "region") %>%
  left_join(., mod_reg_key) %>% mutate(type = "truth")

median_phi_nu <- apply(post$phi_nu, c(2,3), median)
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

save_nu <- full_phi_onetime_nu + scale_fill_brewer(palette = "YlOrRd")

ggsave(paste0('./sim-study/figures/g1/maps/g1_phi-map_kappa-nu-sim_nu_',
              format(as.POSIXlt(Sys.time(), "America/Denver"), "%d%b%Y_%H%M"),
              ".pdf"),
       plot=save_nu,
       device = "pdf")

# xi maps ------
pre_phi_xi <- phi_mat_xi %>% as_tibble() %>% 
  rename_with(., ~ reg_cols) %>%
  mutate(timepoint = 1:all_of(t)) %>%
  pivot_longer(cols = c(1:all_of(r)), values_to = "phi_val", names_to = "region") %>%
  left_join(., mod_reg_key) %>% mutate(type = "truth")

median_phi_xi <- apply(post$phi_xi, c(2,3), median)
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

save_xi <- full_phi_onetime_xi + scale_fill_brewer(palette = "YlOrRd")

ggsave(paste0('./sim-study/figures/g1/maps/g1_phi-map_kappa-nu-sim_xi_',
              format(as.POSIXlt(Sys.time(), "America/Denver"), "%d%b%Y_%H%M"),
              ".pdf"),
       plot=save_xi,
       device = "pdf")

# trace plot ------
MCMCtrace(egpd_fit, params = c("rho1_kappa", "rho2_kappa", "rho1_nu", "rho2_nu", "rho1_xi", "rho2_xi"),
          ind = TRUE,
          gvals = c(rho1_kappa, rho2_kappa, rho1_nu, rho2_nu, rho1_xi, rho2_xi),
          open_pdf = FALSE,
          filename = paste0('./sim-study/figures/g1/trace/g1_trace_all3_rhos',
                            format(as.POSIXlt(Sys.time(), "America/Denver"), "%d%b%Y_%H%M"), ".pdf"))
