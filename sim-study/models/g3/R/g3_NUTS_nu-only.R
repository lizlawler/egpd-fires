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
source("./sim-study/models/g3/R/g3_data_nu-only.R")

# run sampling
egpd_init <- stan_model('./sim-study/models/g3/stan/g3_nu-only.stan')
egpd_fit <- sampling(egpd_init, 
                     data = toy_data, 
                     iter = 500,
                     chains = 1,
                     refresh = 50)

saveRDS(egpd_fit, file = "./sim-study/models/g3/stan-fits/g3_nu-only.RDS")

## ------
post <- rstan::extract(egpd_fit, pars = c("beta_nu", "phi_nu"))
median_nu <- apply(post$beta_nu, c(2,3), median)

post_nu_effects_df <- matrix(NA, r, t)
for(i in 1:r) {
  post_nu_effects_df[i,] <- X_full[i, , ] %*% median_nu[, i]
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
ggsave(paste0('./sim-study/figures/g3/effects/g3_post-effects-nu_',
              format(as.POSIXlt(Sys.time(), "America/Denver"), "%d%b%Y_%H%M"),
              ".pdf"),
       plot=post_effects_nu,
       device = "pdf")


## maps of phi -----
ecoregions_geom <- ecoregions %>% filter(!NA_L2NAME == "UPPER GILA MOUNTAINS (?)")
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

ggsave(paste0('./sim-study/figures/g3/maps/g3_phi-map-nu_',
              format(as.POSIXlt(Sys.time(), "America/Denver"), "%d%b%Y_%H%M"),
              ".pdf"),
       plot=full_phi_onetime_nu,
       device = "pdf")

# save traceplot-------
MCMCtrace(egpd_fit, params = c("beta_nu", "phi_nu", "rho1_nu", "rho2_nu"),
          ind = TRUE,
          gvals = c(betas_nu, phi_mat_nu, rho1, rho2),
          open_pdf = FALSE,
          filename = paste0('./sim-study/figures/g3/trace/g3_trace-nu_',
                            format(as.POSIXlt(Sys.time(), "America/Denver"), "%d%b%Y_%H%M"), ".pdf"))
