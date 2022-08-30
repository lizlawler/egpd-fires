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
source("./toy-sim/zip/R/zip_data.R")

# run sampling
egpd_init <- stan_model('./toy-sim/zip/stan/zip_icar.stan')
egpd_fit <- sampling(egpd_init, 
                     data = toy_data, 
                     iter = 1000,
                     chains = 3,
                     refresh = 50)

# save traceplot
MCMCtrace(egpd_fit, params = c("beta_lambda", "phi", "rho1", "rho2"),
          ind = TRUE,
          gvals = c(betas_lambda, phi_mat, rho1, rho2),
          open_pdf = FALSE,
          filename = paste0('./toy-sim/figures/zip_trace-lambda_',
                            format(as.POSIXlt(Sys.time(), "America/Denver"), "%d%b%Y_%H%M"), ".pdf"))

# pre and post effects plot
post <- rstan::extract(egpd_fit, pars = c('beta_lambda', 'phi'))
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
ggsave(paste0('./toy-sim/figures/zip_post-effects_',
              format(as.POSIXlt(Sys.time(), "America/Denver"), "%d%b%Y_%H%M"),
              ".pdf"),
       plot=post_effects,
       device = "pdf")


## maps of phi -------
ecoregions_geom <- ecoregions %>% filter(!NA_L2NAME == "UPPER GILA MOUNTAINS (?)")
pre_phi <- phi_mat %>% as_tibble() %>% 
  rename_with(., ~ reg_cols) %>%
  mutate(timepoint = 1:all_of(t)) %>%
  pivot_longer(cols = c(1:all_of(r)), values_to = "phi_val", names_to = "region") %>%
  left_join(., mod_reg_key) %>% mutate(type = "truth")

median_phi <- apply(post$phi, c(2,3), median)
post_phi <- median_phi %>% as_tibble() %>%
  rename_with(., ~ reg_cols) %>%
  mutate(timepoint = 1:all_of(t)) %>%
  pivot_longer(cols = c(1:all_of(r)), values_to = "phi_val", names_to = "region") %>%
  left_join(., mod_reg_key) %>% mutate(type = "sim")

phi_full <- rbind(pre_phi, post_phi) %>% mutate(type = factor(type, levels = c("truth", "sim")))
joined_l3_phi_full <- right_join(ecoregions_geom, phi_full)
breaks <- classIntervals(c(min(joined_l3_phi_full$phi_val) - .00001, joined_l3_phi_full$phi_val), n=5, style = "quantile")
joined_l3_phi_full <- joined_l3_phi_full %>% mutate(phi_cat = cut(phi_val, unique(breaks$brks)))

full_phi_onetime <- ecoregions_geom %>%
  ggplot() +
  geom_sf(size = .1, fill = 'white') + 
  geom_sf(data = joined_l3_phi_full[joined_l3_phi_full$timepoint == 99, ], 
          aes(fill=phi_cat), alpha = 0.6, lwd = 0, inherit.aes = FALSE) + 
  facet_grid(type ~ .) +
  theme_minimal() + 
  theme(panel.grid.major = element_line(colour = "lightgrey"))
ggsave(paste0('./toy-sim/figures/zip_phi-map_',
              format(as.POSIXlt(Sys.time(), "America/Denver"), "%d%b%Y_%H%M"),
              ".pdf"),
       plot=full_phi_onetime,
       device = "pdf")
