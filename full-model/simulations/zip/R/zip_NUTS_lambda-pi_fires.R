rm(list = setdiff(ls(), "stan_data"))
gc()

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


# start time, for identification purposes later
st_time <- format(as.POSIXlt(Sys.time(), "America/Denver"), "%d-%b-%Y_%H%M")
stan_data <- readRDS(file = "full-model/simulations/zip/data/stan_data.RDS")

# run sampling
egpd_init <- stan_model('./full-model/simulations/zip/stan/zip_lambda-pi_fires.stan')
beta_rho_inits <- list(beta_lambda = matrix(0, 37, 84), beta_pi = matrix(0, 37, 84), 
                       rho1_lambda = 0.54, rho2_lambda = 0.45, rho1_pi = 0.54, rho2_pi = 0.45)
inits_list3 <- rep(list(beta_rho_inits), 3)
egpd_fit <- sampling(egpd_init, 
                     data = stan_data, 
                     iter = 1000,
                     chains = 3,
                     init = inits_list3,
                     refresh = 50)

end_time <- format(as.POSIXlt(Sys.time(), "America/Denver"), "%H%M")

saveRDS(egpd_fit, file = paste0("./full-model/simulations/zip/stan-fits/zip_lambda-pi_fires", st_time, "_", end_time, ".RDS"))

# saveRDS(count_data, file = paste0("./sim-study/models/zip/data/zip_lambda-pi_", st_time, "_", end_time, ".RDS"))

MCMCtrace(egpd_fit, params = c("rho1_lambda", "rho2_lambda", "rho1_pi", "rho2_pi"),
          ind = TRUE,
          open_pdf = FALSE,
          filename = paste0('./full-model/figures/zip/trace/zip_trace-lambda-pi_fires', st_time, "_", end_time, ".pdf"))

# run for 2000 iterations to compare the two
rm(list = setdiff(ls(), c("stan_data", "egpd_init", "inits_list3")))
gc()
st_time <- format(as.POSIXlt(Sys.time(), "America/Denver"), "%d%b%Y_%H%M")
egpd_fit_2000 <- sampling(egpd_init, 
                     data = stan_data, 
                     iter = 2000,
                     chains = 3,
                     init = inits_list3,
                     refresh = 50)
end_time <- format(as.POSIXlt(Sys.time(), "America/Denver"), "%H%M")
saveRDS(egpd_fit_2000, file = paste0("./full-model/simulations/zip/stan-fits/zip_lambda-pi_fires_2000iter", st_time, "_", end_time, ".RDS"))
MCMCtrace(egpd_fit_2000, params = c("rho1_lambda", "rho2_lambda", "rho1_pi", "rho2_pi"),
          ind = TRUE,
          open_pdf = FALSE,
          filename = paste0('./full-model/figures/zip/trace/zip_trace-lambda-pi_fires_2000iter', st_time, "_", end_time, ".pdf"))

quit(save = "no")


# MCMCtrace(egpd_fit, params = c("rho1_lambda", "rho2_lambda", "rho1_pi", "rho2_pi"),
#           ind = TRUE)
# 
# MCMCtrace(egpd_fit, params = c("beta_lambda", "beta_pi", "phi_lambda", "phi_pi"),
#           ind = TRUE)

# saveRDS(egpd_fit, file = "./sim-study/models/zip/stan-fits/zip_lambda-pi.RDS")

# pre and post effects plots ---------
post <- rstan::extract(egpd_fit, pars = c('beta_lambda', 'beta_pi', 'phi_lambda', 'phi_pi'))

## lambda ----
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
ggsave(paste0('./sim-study/figures/zip/effects/zip_lambda-pi_lambda-effects',
              format(as.POSIXlt(Sys.time(), "America/Denver"), "%d%b%Y_%H%M"),
              ".pdf"),
       plot=post_effects,
       device = "pdf")

## pi ----
median_pi <- apply(post$beta_pi, c(2,3), median)
post_pi_effects_df <- matrix(NA, r, t)
for(i in 1:r) {
  post_pi_effects_df[i,] <- X_full[i, , ] %*% median_pi[, i]
}

post_pi <- t(post_pi_effects_df) %>% as_tibble() %>% 
  rename_with(., ~ reg_cols) %>% 
  mutate(time = c(1:t)) %>%
  pivot_longer(cols = c(1:all_of(r)), values_to = "effect", names_to = "region") %>%
  left_join(., X_long) %>%
  left_join(., mod_reg_key) %>% mutate(type = "sim")

pi_full <- rbind(pi_effects, post_pi) %>% mutate(type = factor(type, levels = c("truth", "sim")))

post_effects <- ggplot(pi_full, aes(x=linear, y=effect, group = region)) + 
  geom_line(aes(linetype=NA_L1CODE, color = NA_L2CODE)) + 
  facet_grid(type ~ .)
ggsave(paste0('./sim-study/figures/zip/effects/zip_lambda-pi_pi-effects',
              format(as.POSIXlt(Sys.time(), "America/Denver"), "%d%b%Y_%H%M"),
              ".pdf"),
       plot=post_effects,
       device = "pdf")



## maps of phi -------
## lambda ----
ecoregions_geom <- ecoregions %>% filter(!NA_L2NAME == "UPPER GILA MOUNTAINS (?)")

pre_phi_lambda <- phi_mat_lambda %>% as_tibble() %>% 
  rename_with(., ~ reg_cols) %>%
  mutate(timepoint = 1:all_of(t)) %>%
  pivot_longer(cols = c(1:all_of(r)), values_to = "phi_val", names_to = "region") %>%
  left_join(., mod_reg_key) %>% mutate(type = "truth")

median_phi_lambda <- apply(post$phi_lambda, c(2,3), median)
post_phi_lambda <- median_phi_lambda %>% as_tibble() %>%
  rename_with(., ~ reg_cols) %>%
  mutate(timepoint = 1:all_of(t)) %>%
  pivot_longer(cols = c(1:all_of(r)), values_to = "phi_val", names_to = "region") %>%
  left_join(., mod_reg_key) %>% mutate(type = "sim")

phi_lambda_full <- rbind(pre_phi_lambda, post_phi_lambda) %>% mutate(type = factor(type, levels = c("truth", "sim")))
joined_l3_phi_lambda_full <- right_join(ecoregions_geom, phi_lambda_full)
breaks <- classIntervals(c(min(joined_l3_phi_lambda_full$phi_val) - .00001, joined_l3_phi_lambda_full$phi_val), n=5, style = "quantile")
joined_l3_phi_lambda_full <- joined_l3_phi_lambda_full %>% mutate(phi_cat = cut(phi_val, unique(breaks$brks)))

full_phi_onetime_lambda <- ecoregions_geom %>%
  ggplot() +
  geom_sf(size = .1, fill = 'white') + 
  geom_sf(data = joined_l3_phi_lambda_full[joined_l3_phi_lambda_full$timepoint == 99, ], 
          aes(fill=phi_cat), alpha = 0.6, lwd = 0, inherit.aes = FALSE) + 
  facet_grid(type ~ .) +
  theme_minimal() + 
  theme(panel.grid.major = element_line(colour = "lightgrey"))
ggsave(paste0('./sim-study/figures/zip/maps/zip_phi-map_lambda-pi_lambda',
              format(as.POSIXlt(Sys.time(), "America/Denver"), "%d%b%Y_%H%M"),
              ".pdf"),
       plot=full_phi_onetime_lambda,
       device = "pdf")

## pi ----
ecoregions_geom <- ecoregions %>% filter(!NA_L2NAME == "UPPER GILA MOUNTAINS (?)")

pre_phi_pi <- phi_mat_pi %>% as_tibble() %>% 
  rename_with(., ~ reg_cols) %>%
  mutate(timepoint = 1:all_of(t)) %>%
  pivot_longer(cols = c(1:all_of(r)), values_to = "phi_val", names_to = "region") %>%
  left_join(., mod_reg_key) %>% mutate(type = "truth")

median_phi_pi <- apply(post$phi_pi, c(2,3), median)
post_phi_pi <- median_phi_pi %>% as_tibble() %>%
  rename_with(., ~ reg_cols) %>%
  mutate(timepoint = 1:all_of(t)) %>%
  pivot_longer(cols = c(1:all_of(r)), values_to = "phi_val", names_to = "region") %>%
  left_join(., mod_reg_key) %>% mutate(type = "sim")

phi_pi_full <- rbind(pre_phi_pi, post_phi_pi) %>% mutate(type = factor(type, levels = c("truth", "sim")))
joined_l3_phi_pi_full <- right_join(ecoregions_geom, phi_pi_full)
breaks <- classIntervals(c(min(joined_l3_phi_pi_full$phi_val) - .00001, joined_l3_phi_pi_full$phi_val), n=5, style = "quantile")
joined_l3_phi_pi_full <- joined_l3_phi_pi_full %>% mutate(phi_cat = cut(phi_val, unique(breaks$brks)))

full_phi_onetime_pi <- ecoregions_geom %>%
  ggplot() +
  geom_sf(size = .1, fill = 'white') + 
  geom_sf(data = joined_l3_phi_pi_full[joined_l3_phi_pi_full$timepoint == 99, ], 
          aes(fill=phi_cat), alpha = 0.6, lwd = 0, inherit.aes = FALSE) + 
  facet_grid(type ~ .) +
  theme_minimal() + 
  theme(panel.grid.major = element_line(colour = "lightgrey"))
ggsave(paste0('./sim-study/figures/zip/maps/zip_phi-map_lambda-pi_pi',
              format(as.POSIXlt(Sys.time(), "America/Denver"), "%d%b%Y_%H%M"),
              ".pdf"),
       plot=full_phi_onetime_pi,
       device = "pdf")

# save traceplot
MCMCtrace(egpd_fit, params = c("beta_lambda", "beta_pi", "phi_lambda", "phi_pi"),
          ind = TRUE,
          gvals = c(betas_lambda, phi_mat, rho1, rho2),
          open_pdf = FALSE,
          filename = paste0('./sim-study/figures/zip/trace/zip_trace-lambda-pi',
                            format(as.POSIXlt(Sys.time(), "America/Denver"), "%d%b%Y_%H%M"), ".pdf"))

