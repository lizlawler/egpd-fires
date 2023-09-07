library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(tidyverse)
library(stringr)
library(posterior)
library(lubridate)
library(sf)
library(classInt)
library(RColorBrewer)

# following code is for extracting from the actual model fit -------
burn_fits <- paste0("full-model/fire-sims/joint/sigma-ri/csv-fits/",
                    list.files(path = "full-model/fire-sims/joint/sigma-ri/csv-fits",
                               pattern = ".csv", recursive = TRUE))
best_fit <- burn_fits[grepl("theta-time", burn_fits)]
# nfits <- length(burn_fits)/3
# fit_groups <- vector(mode = "list", nfits)
# for(i in 1:nfits) {
#   fit_groups[[i]] <- burn_fits[(3*i-2):(3*i)]
# }

best_group <- list(best_fit)
# burn_names <- lapply(fit_groups, function(x) str_remove(str_remove(basename(x[1]), "_\\d{2}\\w{3}2023_\\d{4}_\\d{1}.csv"), "joint_")) %>% unlist()
# 
extraction <- function(file_group, burn_name) {
  object <- as_cmdstan_fit(file_group)
  # betas <- object$draws(variables = "beta")
  reg <- object$draws(variables = "reg")
  ri_init <- object$draws(variables = "ri_init")
  pi_prob <- object$draws(variables = "pi_prob")
  lambda <- object$draws(variables = "lambda")
  delta <- object$draws(variables = "delta")
  theta <- object$draws(variables = "theta")
  gamma <- object$draws(variables = "gamma")
  phi <- object$draws(variables = "phi")
  temp <- list(reg, ri_init, pi_prob, lambda, delta, theta, gamma, phi)
  names(temp) <- c("reg", "ri_init", "pi_prob", "lambda", "delta", "theta", "gamma", "phi")
  assign(burn_name, temp, parent.frame())
  rm(object)
  gc()
}

extraction(best_group[[1]], "best_fit")

kappa_vals <- best_fit$reg %>% as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>% 
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>% 
  separate_wider_delim(cols = "name", delim = ",", names = c("time", "region")) %>%
  mutate(time = as.numeric(gsub("reg\\[", "", time)),
         region = as.numeric(gsub("\\]", "", region)),
         kappa = exp(value)) %>% select(-value)
# saveRDS(kappa_vals, file = "full-model/figures/paper/mcmc_draws/kappa_vals.RDS")
  
rand_int <- best_fit$ri_init %>% as_draws_df() %>%
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
         value = exp(value)) %>% 
  pivot_wider(names_from = param, values_from = value)
saveRDS(rand_int, file = "full-model/figures/paper/mcmc_draws/rand_int.RDS")

lambda <- best_fit$lambda %>% as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>% 
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>% 
  separate_wider_delim(cols = "name", delim = ",", names = c("time", "region")) %>%
  mutate(time = as.numeric(gsub("lambda\\[", "", time)),
         region = as.numeric(gsub("\\]", "", region)),
         param = "lambda") %>%
  pivot_wider(names_from = param, values_from = value)
saveRDS(lambda, file = "full-model/figures/paper/mcmc_draws/lambda.RDS")

phi <- best_fit$phi %>% as_draws_df() %>%
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
saveRDS(phi, file = "full-model/figures/paper/mcmc_draws/phi.RDS")

pi_prob <- best_fit$pi_prob %>% as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>% 
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  mutate(region = as.character(gsub("pi_prob\\[", "", name)),
         region = as.numeric(gsub("\\]", "", region)),
         name = "pi") %>% 
  pivot_wider(names_from = name, values_from = value)
saveRDS(pi_prob, file = "full-model/figures/paper/mcmc_draws/pi_prob.RDS")

gamma <- best_fit$gamma %>% as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>% 
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  mutate(region = as.character(gsub("gamma\\[", "", name)),
         region = as.numeric(gsub("\\]", "", region)),
         name = "gamma") %>% 
  pivot_wider(names_from = name, values_from = value)
saveRDS(gamma, file = "full-model/figures/paper/mcmc_draws/gamma.RDS")

theta <- best_fit$theta %>% as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>% 
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  mutate(timepoint = as.character(gsub("theta\\[", "", name)),
         timepoint = as.numeric(gsub("\\]", "", timepoint)),
         name = "theta") %>% 
  pivot_wider(names_from = name, values_from = value)
saveRDS(theta, file = "full-model/figures/paper/mcmc_draws/theta.RDS")

## read in region key --------
region_key <- readRDS(file = "./full-model/data/processed/region_key.rds")
full_reg_key <- as_tibble(region_key) %>% 
  mutate(region = c(1:84),
         NA_L2CODE = as.factor(NA_L2CODE),
         NA_L1CODE = as.factor(NA_L1CODE),
         NA_L3CODE = as.factor(NA_L3CODE),
         NA_L1NAME = as.factor(str_to_title(NA_L1NAME)))

date_seq <- seq(as.Date("1995-01-01"), by = "1 month", length.out = 252) %>% as_tibble() %>% rename(date = value)
time_df <- date_seq %>% mutate(time = 1:252)

## create time series plots of phi values for kappa and lambda -------
phi_kappa <- phi %>% select(-lambda) %>% group_by(time, region) %>% summarize(kappa = median(kappa)) %>% ungroup()
p <- phi_kappa %>% filter(time <= 252) %>% left_join(full_reg_key) %>% left_join(time_df) %>%
  ggplot(aes(x = date, y = kappa, color = NA_L2CODE)) + 
  geom_line(linewidth = 0.5) + 
  scale_x_date(name = "Year (1995-2016)", date_breaks = "5 years", date_labels = "%Y") +
  facet_wrap(. ~ NA_L1NAME, nrow = 2) + 
  theme_classic()
ggsave("full-model/figures/paper/phi_kappa_allregs_overtime.pdf", p, width = 15, dpi = 320, bg = "white")

phi_lambda <- phi %>% select(-kappa) %>% group_by(time, region) %>% summarize(lambda = median(lambda)) %>% ungroup()
p <- phi_lambda %>% filter(time <= 252) %>% left_join(full_reg_key) %>% left_join(time_df) %>%
  ggplot(aes(x = date, y = lambda, color = NA_L2CODE)) + 
  geom_line(linewidth = 0.5) + 
  scale_x_date(name = "Year (1995-2016)", date_breaks = "5 years", date_labels = "%Y") +
  facet_wrap(. ~ NA_L1NAME, nrow = 2) + 
  theme_classic()
ggsave("full-model/figures/paper/phi_lambda_allregs_overtime.pdf", p, width = 15, dpi = 320, bg = "white")

phi_etf <- phi %>% left_join(full_reg_key) %>% 
  filter(NA_L1NAME == "Eastern Temperate Forests") %>% 
  group_by(time, region) %>% 
  summarize(lambda = median(lambda), kappa = median(kappa)) %>% 
  ungroup() %>% left_join(full_reg_key)

years <- seq(1, 252, by = 12)[-1]
p <- phi_etf %>% filter(time <= 252) %>% left_join(time_df) %>% 
  pivot_longer(cols = c("lambda", "kappa"), names_to = "param", values_to = "value") %>%
  ggplot(aes(x = date, y = value, color = NA_L2CODE, alpha = 0.5)) + 
  geom_line(linewidth = 0.5) + 
  # geom_vline(xintercept = years, col = "darkgrey", linewidth = 0.5) + 
  scale_x_date(name = "Year (1995-2016)", date_breaks = "5 years", date_labels = "%Y") + 
  facet_wrap(. ~ param, scales = "free_y") +
  theme_classic() + theme(legend.position = "none")
ggsave("full-model/figures/paper/phi_etf_overtime.pdf", p, width = 15, dpi = 320, bg = "white")

p <- phi_etf %>% filter(time <= 252) %>% 
  ggplot(aes(x = time, y = kappa, color = NA_L2CODE, alpha = 0.5)) + 
  geom_line(linewidth = 0.5) + 
  geom_vline(xintercept = years, col = "darkgrey", linewidth = 0.5) +
  theme_classic() + theme(legend.position = "none")
ggsave("full-model/figures/paper/phi_kappa_etf_overtime.pdf", p, width = 15, dpi = 320, bg = "white")

p <- phi_etf %>% filter(time <= 252) %>% 
  ggplot(aes(x = time, y = lambda, color = NA_L2CODE, alpha = 0.5)) + 
  geom_line(linewidth = 0.5) + 
  geom_vline(xintercept = years, col = "darkgrey", linewidth = 0.5) +
  theme_classic() + theme(legend.position = "none")
ggsave("full-model/figures/paper/phi_lambda_etf_overtime.pdf", p, width = 15, dpi = 320, bg = "white")

# create return level plots ------
# egpd functions
pegpd <- function(y, kappa, sigma, xi) (1 - (1 + xi * (y/sigma))^(-1/xi))^kappa
qegpd <- function(p, kappa, sigma, xi) (sigma/xi) * ( (1 - p^(1/kappa) )^-xi - 1)
# expected counts from ZINB params
exp_count <- function(pi, lambda) { # expected counts
  pi_prob <- exp(pi)/(1+exp(pi))
  return((1-pi_prob) * exp(lambda))
}

# exceedance probabilities
rlevel <- function(N, kappa, sigma, xi, eta) {
  p <- ((N-1)/N)^(1/(12 * eta)) * (1-pegpd(1.001, kappa, sigma, xi)) + pegpd(1.001, kappa, sigma, xi)
  return(qegpd(p, kappa, sigma, xi))
}

# high quantiles
high_quant <- function(N, kappa, sigma, xi) {
  p <- ((N-1)/N) * (1-pegpd(1.001, kappa, sigma, xi)) + pegpd(1.001, kappa, sigma, xi)
  return(qegpd(p, kappa, sigma, xi))
}

all_params <- kappa_vals %>% 
  left_join(rand_int) %>% 
  left_join(lambda) %>% 
  left_join(pi_prob) %>%
  mutate(eta = exp_count(pi, lambda))

# calculate expected burn area for the 98th quantile given there are eta fires in that ecoregion
returns <- all_params %>% 
  mutate(yr50 = rlevel(50, kappa, sigma, xi, eta)*1000*0.405) # rescale back to 1000s of acres; convert to hectares
returns_summary <- returns %>% group_by(time, region) %>% 
  summarize(med50 = median(yr50[is.finite(yr50)]), 
            lower = quantile(yr50[is.finite(yr50)], probs = 0.025), 
            upper = quantile(yr50[is.finite(yr50)], probs = 0.975)) %>%
  ungroup()

returns_regional <- returns_summary %>% 
  left_join(full_reg_key) %>% 
  left_join(time_df)

p <- returns_regional %>% ggplot() + 
  geom_ribbon(aes(x=date, ymin=lower, ymax=upper, group = region, fill = NA_L1NAME, alpha = 0.5)) +
  geom_line(aes(x=date, y=med50, group = region, alpha = 0.5), linewidth = 0.5, color = 'darkgrey') + scale_y_log10() +
  scale_x_date(name = "Year (1995-2016)", date_breaks = "5 years", date_labels = "%Y") + 
  ylab("Expected burn area (ha)") +
  facet_wrap(. ~ NA_L1NAME, nrow = 2) +
  theme_classic() + theme(legend.position = "none") + 
  theme(strip.text.x = element_text(size = rel(0.9)))

file_name <- "full-model/figures/paper/50yr_returns.pdf"
ggsave(file_name, p, dpi = 320, bg = "white", width = 17, height = 9)

## calculate 98th quantile expected burn areas ----------
level98_areas <- all_params %>% 
  mutate(yr50 = high_quant(50, kappa, sigma, xi)*1000*0.405) # rescale back to 1000s of acres; convert to hectares
level98_areas_summary <- level98_areas %>% group_by(time, region) %>% 
  summarize(med50 = median(yr50[is.finite(yr50)]), 
            lower = quantile(yr50[is.finite(yr50)], probs = 0.025), 
            upper = quantile(yr50[is.finite(yr50)], probs = 0.975)) %>%
  ungroup()

level98_areas_regional <- level98_areas_summary %>% 
  left_join(full_reg_key) %>% 
  left_join(time_df)

p <- level98_areas_regional %>% ggplot() + 
  geom_ribbon(aes(x=date, ymin=lower, ymax=upper, group = region, fill = NA_L1NAME, alpha = 0.5)) +
  geom_line(aes(x=date, y=med50, group = region), linewidth = 0.5, color = 'darkgrey') + scale_y_log10() +
  scale_x_date(name = "Year (1995-2016)", date_breaks = "5 years", date_labels = "%Y") + 
  ylab("Expected burn area (ha)") +
  facet_wrap(. ~ NA_L1NAME, nrow = 2) +
  theme_classic() + theme(legend.position = "none") + 
  theme(strip.text.x = element_text(size = rel(0.9)))

file_name <- "full-model/figures/paper/98th_quant_burnareas.pdf"
ggsave(file_name, p, dpi = 320, bg = "white", width = 17, height = 9)

## area-weighted average of burn area exceedances -------
ecoregions <- read_rds(file = "ecoregions.RDS")
ecoregions_geom <- ecoregions %>% filter(!NA_L2NAME == "UPPER GILA MOUNTAINS (?)") %>% 
  mutate(NA_L2CODE = as.factor(NA_L2CODE),
         NA_L1CODE = as.factor(NA_L1CODE),
         NA_L3CODE = as.factor(NA_L3CODE),
         NA_L1NAME = as.factor(str_to_title(NA_L1NAME)))

eco_areas <- ecoregions_geom %>% as_tibble() %>% group_by(NA_L3CODE) %>% summarise(area = sum(Shape_Area))

returns_level1 <- returns %>% select(c(draw, time, region, yr50)) %>% 
  left_join(full_reg_key) %>% 
  left_join(eco_areas) %>% 
  filter(yr50 != Inf) %>%
  group_by(NA_L1NAME, time, draw) %>% 
  summarize(wmean = weighted.mean(yr50, area)) %>% 
  ungroup() %>% 
  group_by(time, NA_L1NAME) %>%
  summarize(med50 = median(wmean), 
            lower = quantile(wmean, probs = 0.025), 
            upper = quantile(wmean, probs = 0.975))

p <- returns_level1 %>%
  left_join(time_df) %>%
  ggplot() + 
  geom_ribbon(aes(x=date, ymin=lower, ymax=upper, group = NA_L1NAME, fill = NA_L1NAME, alpha = 0.5)) +
  geom_line(aes(x=date, y=med50, group = NA_L1NAME, alpha = 0.5), linewidth = 0.5, color = 'darkgrey') + scale_y_log10() +
  scale_x_date(name = "Year (1995-2016)", date_breaks = "5 years", date_labels = "%Y") + 
  ylab("Expected burn area (ha)") +
  facet_wrap(. ~ NA_L1NAME, nrow = 2) +
  theme_classic() + theme(legend.position = "none") + 
  theme(strip.text.x = element_text(size = rel(0.9)))

file_name <- "full-model/figures/paper/50yr_returns_level1.pdf"
ggsave(file_name, p, dpi = 320, bg = "white", width = 17, height = 9)

## time-averaged kappa, then high quantile of burn area
kappa_temp_avg <- kappa_vals %>% group_by(region, draw) %>% summarize(kappa = mean(kappa)) %>% ungroup()
all_params_notime <- kappa_temp_avg %>% 
  left_join(rand_int)
level98_timeavg <- all_params_notime %>% 
  mutate(yr50 = high_quant(50, kappa, sigma, xi)*0.405) # convert to hectares (1000s of hectares)
level98_timeavg_summary <- level98_timeavg %>% group_by(region) %>%
  summarize(med50 = median(yr50[is.finite(yr50)]), 
            lower = quantile(yr50[is.finite(yr50)], probs = 0.025), 
            upper = quantile(yr50[is.finite(yr50)], probs = 0.975)) %>%
  ungroup()
level98_timeavg_regions <- level98_timeavg_summary %>% left_join(full_reg_key)

level95_timeavg <- all_params_notime %>% 
  mutate(yr20 = high_quant(20, kappa, sigma, xi)*0.405) # convert to hectares (1000s of hectares)
level95_timeavg_summary <- level95_timeavg %>% group_by(region) %>%
  summarize(med20 = median(yr20[is.finite(yr20)]), 
            lower = quantile(yr20[is.finite(yr20)], probs = 0.025), 
            upper = quantile(yr20[is.finite(yr20)], probs = 0.975)) %>%
  ungroup()
level95_timeavg_regions <- level95_timeavg_summary %>% left_join(full_reg_key)

level75_timeavg <- all_params_notime %>% 
  mutate(yr4 = high_quant(4, kappa, sigma, xi)*0.405)
level75_timeavg_summary <- level75_timeavg %>% group_by(region) %>%
  summarize(med4 = median(yr4[is.finite(yr4)]), 
            lower = quantile(yr4[is.finite(yr4)], probs = 0.025), 
            upper = quantile(yr4[is.finite(yr4)], probs = 0.975)) %>%
  ungroup()
level75_timeavg_regions <- level75_timeavg_summary %>% left_join(full_reg_key)

level50_timeavg <- all_params_notime %>% 
  mutate(yr2 = high_quant(2, kappa, sigma, xi)*0.405) # rescale back to 1000s of acres; convert to hectares
level50_timeavg_summary <- level50_timeavg %>% group_by(region) %>%
  summarize(med2 = median(yr2[is.finite(yr2)]), 
            lower = quantile(yr2[is.finite(yr2)], probs = 0.025), 
            upper = quantile(yr2[is.finite(yr2)], probs = 0.975)) %>%
  ungroup()
level50_timeavg_regions <- level50_timeavg_summary %>% left_join(full_reg_key)

## create maps for each quantile level, using same interval break categories -------
allquantlevels <- c(level50_timeavg_summary$med2, level75_timeavg_summary$med4, level95_timeavg_summary$med20, level98_timeavg_summary$med50)
breaks <- classIntervals(c(min(allquantlevels) - .00001, allquantlevels), style = 'quantile',intervalClosure = 'left')

eco_burns98 <- ecoregions_geom %>% 
  left_join(level98_timeavg_regions) %>% 
  mutate(burns_cat = cut(med50, unique(breaks$brks)))
p <- ecoregions_geom %>%
  ggplot() +
  geom_sf(size = .1, fill = 'white') +
  geom_sf(data = eco_burns98,
          aes(fill=burns_cat), alpha = 0.6, lwd = 0, inherit.aes = FALSE) + 
  scale_fill_brewer(palette = 'YlOrRd') +
  theme_void() 
ggsave("full-model/figures/paper/onekappa_98th_quant_map.pdf", dpi = 320, bg ='white')

eco_burns95 <- ecoregions_geom %>% 
  left_join(level95_timeavg_regions) %>% 
  mutate(burns_cat = cut(med20, unique(breaks$brks)))
p <- ecoregions_geom %>%
  ggplot() +
  geom_sf(size = .1, fill = 'white') +
  geom_sf(data = eco_burns95,
          aes(fill=burns_cat), alpha = 0.6, lwd = 0, inherit.aes = FALSE) + 
  scale_fill_brewer(palette = 'YlOrRd') +
  theme_void() 
ggsave("full-model/figures/paper/onekappa_95th_quant_map.pdf", dpi = 320, bg ='white')

eco_burns75 <- ecoregions_geom %>% 
  left_join(level75_timeavg_regions) %>% 
  mutate(burns_cat = cut(med4, unique(breaks$brks)))
p <- ecoregions_geom %>%
  ggplot() +
  geom_sf(size = .1, fill = 'white') +
  geom_sf(data = eco_burns75,
          aes(fill=burns_cat), alpha = 0.6, lwd = 0, inherit.aes = FALSE) + 
  scale_fill_brewer(palette = 'YlOrRd') +
  theme_void() 
ggsave("full-model/figures/paper/onekappa_75th_quant_map.pdf", dpi = 320, bg ='white')

eco_burns50 <- ecoregions_geom %>% 
  left_join(level50_timeavg_regions) %>% 
  mutate(burns_cat = cut(med2, unique(breaks$brks)))
p <- ecoregions_geom %>%
  ggplot() +
  geom_sf(size = .1, fill = 'white') +
  geom_sf(data = eco_burns50,
          aes(fill=burns_cat), alpha = 0.6, lwd = 0, inherit.aes = FALSE) + 
  scale_fill_brewer(palette = 'YlOrRd') +
  theme_void() 
ggsave("full-model/figures/paper/onekappa_50th_quant_map.pdf", dpi = 320, bg ='white')

## create maps of each parameter (kappa, sigma, xi) using 'quantile' break schema ---------
all_params_notime <- all_params_notime %>% 
  group_by(region) %>% 
  summarize(kappa = mean(kappa),
            sigma = mean(sigma),
            xi = mean(xi)) %>%
  left_join(full_reg_key)
kappa_breaks <- classIntervals(c(min(all_params_notime$kappa) - .00001, all_params_notime$kappa), style = 'quantile', intervalClosure = 'left')
sigma_breaks <- classIntervals(c(min(all_params_notime$sigma) - .00001, all_params_notime$sigma), style = 'quantile', intervalClosure = 'left')
xi_breaks <- classIntervals(c(min(all_params_notime$xi) - .00001, all_params_notime$xi), style = 'quantile',intervalClosure = 'left')
all_params_notime_eco <- ecoregions_geom %>% left_join(all_params_notime) %>%
  mutate(kappa_cat = cut(kappa, unique(kappa_breaks$brks)),
         sigma_cat = cut(sigma, unique(sigma_breaks$brks)),
         xi_cat = cut(xi, unique(xi_breaks$brks)))

p <- ecoregions_geom %>%
  ggplot() +
  geom_sf(size = .1, fill = 'white') +
  geom_sf(data = all_params_notime_eco,
          aes(fill=kappa_cat), alpha = 0.6, lwd = 0, inherit.aes = FALSE) + 
  scale_fill_brewer(palette = 'YlOrRd') +
  theme_void() 
ggsave("full-model/figures/paper/kappa_timeavg_map.pdf", dpi = 320, bg ='white')

p <- ecoregions_geom %>%
  ggplot() +
  geom_sf(size = .1, fill = 'white') +
  geom_sf(data = all_params_notime_eco,
          aes(fill=sigma_cat), alpha = 0.6, lwd = 0, inherit.aes = FALSE) + 
  scale_fill_brewer(palette = 'YlOrRd') +
  theme_void()
ggsave("full-model/figures/paper/sigma_map.pdf", dpi = 320, bg ='white')

p <- ecoregions_geom %>%
  ggplot() +
  geom_sf(size = .1, fill = 'white') +
  geom_sf(data = all_params_notime_eco,
          aes(fill=xi_cat), alpha = 0.6, lwd = 0, inherit.aes = FALSE) + 
  scale_fill_brewer(palette = 'YlOrRd') +
  theme_void()
ggsave("full-model/figures/paper/xi_cat.pdf", dpi = 320, bg ='white')


## create boxplot of predicted counts for entire US annually for all timepoints (holdout and training) and overlay truth ------
betas <- readRDS("~/research/egpd-fires/full-model/fire-sims/model_comparison/extracted_values/joint_sigma-ri_theta-time_gamma-ri_betas.RDS")
stan_data <- readRDS("full-model/data/stan_data_joint.RDS")
X_full_count <- stan_data$X_full_count
area_offset <- stan_data$area_offset %>% as_tibble() %>% mutate(region = as.numeric(1:84)) %>% rename(area = value)

beta_count <- betas$beta_count %>% as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>% 
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>% 
  separate_wider_delim(cols = "name", delim = ",", names = c("covar", "region")) %>%
  mutate(covar = as.numeric(gsub("beta_count\\[", "", covar)),
         region = as.numeric(gsub("\\]", "", region)))

beta_by_reg <- split(beta_count, beta_count$region) %>% 
  lapply(., function(x) {
  x %>% 
    select(-region) %>% 
    pivot_wider(names_from = draw, values_from = value) %>%
    select(-covar) %>%
    as.matrix()})

X_by_reg <- asplit(X_full_count, 1)
reg_comp <- mapply(function(x,y) x %*% y, x = X_by_reg, y = beta_by_reg, SIMPLIFY = FALSE)
reg_comp <- lapply(reg_comp, function(x) {
  x %>% 
  as_tibble() %>% 
  rowid_to_column() %>% 
  rename(time = rowid) %>% 
  pivot_longer(!time, names_to = "draw", values_to = "value") %>%
  mutate(draw = as.numeric(draw))})
names(reg_comp) <- paste0(1:84, "")
reg_comp_allregions <- bind_rows(reg_comp, .id = "region") %>% 
  mutate(region = as.numeric(region),
         time = as.numeric(time),
         draw = as.integer(draw))

count_preds <- phi %>% select(-kappa) %>% 
  rename(phi = lambda) %>% 
  left_join(reg_comp_allregions) %>% rename(reg = value) %>%
  left_join(theta %>% rename(time = timepoint)) %>% 
  left_join(pi_prob) %>% left_join(area_offset) %>% 
  mutate(lambda = reg + phi + area + theta) %>% 
  mutate(preds = exp_count(pi, lambda))

preds_only <- count_preds %>% select(c(draw, time, region, preds))
date_seq <- seq(as.Date("1990-01-01"), by = "1 month", length.out = 372) %>% as_tibble() %>% rename(date = value)
time_df <- date_seq %>% mutate(time = 1:372)
preds_only <- preds_only %>% left_join(time_df) %>% mutate(year = year(date))

all_years <- 1990:2020
first_five <- 1990:1994
last_five <- 2020:2016
test_years <- sort(c(first_five, last_five))
train_years <- setdiff(all_years, test_years)

preds_annual <- preds_only %>% group_by(year, draw) %>% summarize(total_fires = sum(preds)) %>% ungroup()
preds_annual <- preds_annual %>% mutate(train = 
                                          case_when(
                                            year %in% train_years ~ TRUE,
                                            year %in% test_years ~ FALSE))

# read in train counts
y_train_count <- stan_data$y_train_count
date_seq <- seq(as.Date("1995-01-01"), by = "1 month", length.out = 252) %>% as_tibble() %>% rename(date = value)
time_df <- date_seq %>% mutate(time = 1:252)
y_train_count <- y_train_count %>% 
  as_tibble() %>% 
  rowid_to_column(var = "time") %>% 
  pivot_longer(!time, names_to = "region", values_to = "value") %>%
  mutate(region = as.numeric(gsub("V", "", region))) %>%
  left_join(time_df) %>% 
  mutate(year = year(date)) %>%
  group_by(year) %>%
  summarize(true_count = sum(value)) %>%
  ungroup()

# read in holdout counts
y_hold_count <- stan_data$y_hold_count
date_seq <- bind_rows(seq(as.Date("1990-01-01"), by = "1 month", length.out = 60) %>% as_tibble() %>% rename(date = value), 
              seq(as.Date("2016-01-01"), by = "1 month", length.out = 60) %>% as_tibble() %>% rename(date = value))
time_df <- date_seq %>% mutate(time = 1:120)
y_hold_count <- y_hold_count %>% 
  as_tibble() %>% 
  rowid_to_column(var = "time") %>% 
  pivot_longer(!time, names_to = "region", values_to = "value") %>%
  mutate(region = as.numeric(gsub("V", "", region))) %>%
  left_join(time_df) %>% 
  mutate(year = year(date)) %>%
  group_by(year) %>%
  summarize(true_count = sum(value)) %>%
  ungroup()

true_counts <- bind_rows(y_train_count, y_hold_count)

p <- preds_annual %>% 
  ggplot(aes(x = year, y = total_fires, group = year, color = train,)) + 
  geom_boxplot(outlier.size = 0.2) + scale_color_grey(start = 0.4, end = 0.6) +
  geom_point(inherit.aes = FALSE, data = true_counts, aes(x = year, y = true_count), col = "red", size = 0.35) +
  geom_line(inherit.aes = FALSE, data = true_counts, aes(x = year, y = true_count), col = "red", linewidth = 0.35) +
  theme_classic() + theme(legend.position = "none")
ggsave("full-model/figures/paper/counts_preds-vs-truth.pdf")

p <- preds_annual %>% 
  ggplot(aes(x = year, y = total_fires, group = year, color = train,)) + 
  geom_boxplot(outlier.shape = NA) + ylim(NA, 1500) + 
  scale_color_grey(start = 0.4, end = 0.6) +
  geom_point(inherit.aes = FALSE, data = true_counts, aes(x = year, y = true_count), col = "red", size = 0.35) +
  geom_line(inherit.aes = FALSE, data = true_counts, aes(x = year, y = true_count), col = "red", linewidth = 0.35) +
  theme_classic() + theme(legend.position = "none")
ggsave("full-model/figures/paper/counts_preds-vs-truth_no-outliers.pdf")

## time series of theta and map of gamma ---------
date_seq <- seq(as.Date("1990-01-01"), by = "1 month", length.out = 372) %>% as_tibble() %>% rename(date = value)
time_df <- date_seq %>% mutate(time = 1:372)
theta_gamma <- theta %>% 
  full_join(gamma) %>% 
  mutate(burn_effect = gamma * theta) %>% 
  rename(time = timepoint) %>%
  left_join(time_df)

theta_meds <- theta %>% group_by(timepoint) %>% summarize(theta = median(theta)) %>% ungroup() %>% rename(time = timepoint)

burn_effect_meds <- theta_gamma %>% 
  group_by(time, region) %>% 
  summarize(burn_effect = median(burn_effect)) %>%
  ungroup() %>% left_join(time_df) %>% left_join(full_reg_key)

p <- burn_effect_meds %>% filter(time <= 252) %>% ggplot(aes(x = time, y = burn_effect, group = region, color = NA_L2CODE)) + 
  geom_line(linewidth = 0.35, alpha = 0.5) +
  geom_line(inherit.aes = FALSE, data = (theta_meds %>% filter(time <= 252)), aes(x=time, y = theta), col = "darkgrey", linewidth = 0.75, alpha = 1) + 
  theme_classic() + theme(legend.position = "none")
ggsave("full-model/figures/paper/time_series_shared_effect.pdf")

years <- seq(1, 252, by = 12)[-1]
p <- theta_meds %>% filter(time <= 252) %>% ggplot(aes(x = time, y = theta)) + geom_line(col = "blue") + 
  geom_vline(xintercept = years, col = "darkgrey", linewidth = 0.5) +
  theme_classic() + theme(legend.position = "none")
ggsave("full-model/figures/paper/time_series_shared_effect_noregions.pdf")


gamma_map <- gamma %>% group_by(region) %>% summarize(gamma = median(gamma))
breaks <- classIntervals(c(min(gamma_map$gamma) - .00001, gamma_map$gamma), style = 'fixed', 
                         fixedBreaks = c(-1, -0.75, -.50, -0.25, 0, 0.25, 0.5, 0.75, 1, 1.25, 1.5), intervalClosure = 'left')
eco_gamma <- ecoregions_geom %>% left_join(gamma_map %>% left_join(full_reg_key)) %>% mutate(gamma_cat = cut(gamma, unique(breaks$brks)))
p <- ecoregions_geom %>%
  ggplot() +
  geom_sf(size = .1, fill = 'white') +
  geom_sf(data = eco_gamma,
          aes(fill=gamma_cat), alpha = 0.6, lwd = 0, inherit.aes = FALSE) +
  theme_void() + scale_fill_brewer(palette = "Spectral")
ggsave("full-model/figures/paper/gamma_map.pdf", dpi = 320, bg ='white')

p <- 

# 
# er_map_l1 <- ecoregions_geom %>%
#   ggplot() +
#   geom_sf(size = .2, aes(fill = NA_L1NAME)) + 
#   theme_void() +
#   coord_sf(ndiscr = FALSE)
# 
# er_map_l1 <- ecoregions_geom %>%
#   ggplot() +
#   geom_sf(size = .2, fill = "white") +
#   geom_sf(data = ecoregions_geom,
#           aes(fill = NA_L2CODE), alpha = 0.6, lwd = 0, inherit.aes = FALSE, show.legend = FALSE) +
#   geom_sf(data = ecoregions_geom %>% group_by(NA_L1CODE) %>% summarise(),
#           fill = "transparent", lwd = 1, color = "gray20", inherit.aes = FALSE, show.legend = FALSE) +
#   theme_void() +
#   coord_sf(ndiscr = FALSE)
# 
# ggsave("test_map.png", er_map_l1, dpi = 320, bg = "white")
# 
# ri_map <- `sigma_ri_theta-ri_gamma-ri`$ri_matrix %>% as_draws_df() %>%
#   select(-c(".iteration", ".chain")) %>% 
#   pivot_longer(cols = !".draw") %>%
#   rename(draw = ".draw") %>% 
#   separate_wider_delim(cols = "name", delim = ",", names = c("param", "time", "region")) %>% select(-time) %>% distinct() %>%
#   mutate(param = as.character(gsub("ri_matrix\\[", "", param)),
#          region = as.numeric(gsub("\\]", "", region)),
#          param = case_when(param == "1" ~ "sigma",
#                            param == "2" ~ "xi",
#                            TRUE ~ param))
# 
# ri_map <- ri_map %>% group_by(param, region) %>% summarize(med_val = median(value))
# ri_map <- ri_map %>% ungroup()
# xi_map <- ri_map %>% filter(param == "xi") %>% select(-param) %>% mutate(med_val = exp(med_val))
# sigma_map <- ri_map %>% filter(param == "sigma") %>% select(-param) %>% mutate(med_val = exp(med_val))
# 
# gamma_map <- gamma %>% group_by(region) %>% summarize(med_gamma = median(gamma)) %>% ungroup() %>% left_join(full_reg_key)
# eco_gamma <- ecoregions_geom %>% left_join(gamma_map)
# breaks <- classIntervals(c(min(eco_gamma$med_gamma) - .00001, eco_gamma$med_gamma), style = 'fixed', 
#                          fixedBreaks = c(-0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4), intervalClosure = 'left')
# eco_gamma <- eco_gamma %>% mutate(gamma_cat = cut(med_gamma, unique(breaks$brks)))
# p <- ecoregions_geom %>%
#   ggplot() +
#   geom_sf(size = .1, fill = 'white') +
#   geom_sf(data = eco_gamma,
#           aes(fill=gamma_cat), alpha = 0.6, lwd = 0, inherit.aes = FALSE) +
#   theme_void() + scale_fill_brewer(palette = "Spectral")
# ggsave("full-model/figures/paper/gamma_map.png", dpi = 320, bg ='white')
# 
# theta_map <- theta %>% group_by(region) %>% summarize(med_theta = median(theta)) %>% ungroup() %>% left_join(full_reg_key)
# eco_theta <- ecoregions_geom %>% left_join(theta_map)
# breaks <- classIntervals(c(min(eco_theta$med_theta) - .00001, eco_theta$med_theta), style = 'fixed', 
#                          fixedBreaks = c(-3, -2, -1, 0, 1, 2), intervalClosure = 'left')
# eco_theta <- eco_theta %>% mutate(theta_cat = cut(med_theta, unique(breaks$brks)))
# p <- ecoregions_geom %>%
#   ggplot() +
#   geom_sf(size = .1, fill = 'white') +
#   geom_sf(data = eco_theta,
#           aes(fill=theta_cat), alpha = 0.6, lwd = 0, inherit.aes = FALSE) +
#   theme_void() + scale_fill_brewer(palette = "Spectral")
# ggsave("full-model/figures/paper/theta_map.png", dpi = 320, bg ='white')
# 
# 
# 
# xi_vals_regions <- xi_map %>% left_join(full_reg_key)
# eco_xi <- ecoregions_geom %>% 
#   mutate(NA_L2CODE = as.factor(NA_L2CODE), 
#          NA_L1CODE = as.factor(NA_L1CODE), 
#          NA_L3CODE = as.factor(NA_L3CODE)) %>% 
#   left_join(xi_vals_regions)
# 
# breaks <- classIntervals(c(min(eco_xi$med_val) - .00001, eco_xi$med_val), style = 'fixed', 
#                          fixedBreaks = c(0, 0.2, 0.4, 0.6, 0.8, 2.0), intervalClosure = 'left')
# 
# eco_xi_cat <- eco_xi %>% 
#   mutate(xi_cat = cut(med_val, unique(breaks$brks)))
# 
# p <- ecoregions_geom %>%
#   ggplot() +
#   geom_sf(size = .1, fill = 'white') +
#   geom_sf(data = eco_xi_cat,
#           aes(fill=xi_cat), alpha = 0.6, lwd = 0, inherit.aes = FALSE) +
#   theme_void() + scale_fill_brewer(palette = 'YlOrRd')
# ggsave("full-model/figures/paper/xi_map.png", dpi = 320, bg ='white')
# 
# 
# sigma_vals_regions <- sigma_map %>% left_join(full_reg_key)
# eco_sigma <- ecoregions_geom %>% 
#   mutate(NA_L2CODE = as.factor(NA_L2CODE), 
#          NA_L1CODE = as.factor(NA_L1CODE), 
#          NA_L3CODE = as.factor(NA_L3CODE)) %>% 
#   left_join(sigma_vals_regions)
# 
# breaks <- classIntervals(c(min(eco_sigma$med_val) - .00001, eco_sigma$med_val), style = 'fisher', n=5, intervalClosure = 'left')
# eco_sigma_cat <- eco_sigma %>% 
#   mutate(sigma_cat = cut(med_val, unique(breaks$brks)))
# 
# p <- ecoregions_geom %>%
#   ggplot() +
#   geom_sf(size = .1, fill = 'white') +
#   geom_sf(data = eco_sigma_cat,
#           aes(fill=sigma_cat), alpha = 0.6, lwd = 0, inherit.aes = FALSE) +
#   theme_void() + scale_fill_brewer(palette = 'YlOrRd')
# ggsave("full-model/figures/paper/sigma_map.png", dpi = 320, bg ='white')
# 
# 
