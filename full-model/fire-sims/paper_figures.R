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
best_fit <- burn_fits[grepl("theta-ri_gamma-ri", burn_fits)]
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

extraction(best_group[[1]], "sigma_ri_theta-ri_gamma-ri")

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
high_quant <- rlevel <- function(N, kappa, sigma, xi) {
  p <- ((N-1)/N) * (1-pegpd(1.001, kappa, sigma, xi)) + pegpd(1.001, kappa, sigma, xi)
  return(qegpd(p, kappa, sigma, xi))
}

kappa_vals <- `sigma_ri_theta-ri_gamma-ri`$reg %>% as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>% 
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>% 
  separate_wider_delim(cols = "name", delim = ",", names = c("time", "region")) %>%
  mutate(time = as.numeric(gsub("reg\\[", "", time)),
         region = as.numeric(gsub("\\]", "", region)),
         kappa = exp(value)) %>% select(-value)
# saveRDS(kappa_vals, file = "full-model/figures/paper/mcmc_draws/kappa_vals.RDS")
  
rand_int <- `sigma_ri_theta-ri_gamma-ri`$ri_init %>% as_draws_df() %>%
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

lambda <- `sigma_ri_theta-ri_gamma-ri`$lambda %>% as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>% 
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>% 
  separate_wider_delim(cols = "name", delim = ",", names = c("time", "region")) %>%
  mutate(time = as.numeric(gsub("lambda\\[", "", time)),
         region = as.numeric(gsub("\\]", "", region)),
         param = "lambda") %>%
  pivot_wider(names_from = param, values_from = value)

phi <- `sigma_ri_theta-ri_gamma-ri`$phi %>% as_draws_df() %>%
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

pi_prob <- `sigma_ri_theta-ri_gamma-ri`$pi_prob %>% as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>% 
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  mutate(region = as.character(gsub("pi_prob\\[", "", name)),
         region = as.numeric(gsub("\\]", "", region)),
         name = "pi") %>% 
  pivot_wider(names_from = name, values_from = value)

gamma <- `sigma_ri_theta-ri_gamma-ri`$gamma %>% as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>% 
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  mutate(region = as.character(gsub("gamma\\[", "", name)),
         region = as.numeric(gsub("\\]", "", region)),
         name = "gamma") %>% 
  pivot_wider(names_from = name, values_from = value)

theta <- `sigma_ri_theta-ri_gamma-ri`$theta %>% as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>% 
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  mutate(region = as.character(gsub("theta\\[", "", name)),
         region = as.numeric(gsub("\\]", "", region)),
         name = "theta") %>% 
  pivot_wider(names_from = name, values_from = value)

phi_onereg <- phi %>% filter(region == 10) %>% group_by(time) %>% summarize(kappa = median(kappa),
                                                                            medlambda = median(lambda)) %>% ungroup()
p <- phi_onereg %>% rename(kappa = medkappa, lambda = medlambda) %>% pivot_longer(cols = !time, names_to = "param", values_to = "value") %>%
  ggplot(aes(x = time, y = value, color = param)) + geom_line() + theme_classic()

ggsave("full-model/figures/paper/phi_reg10_overtime.png", p, dpi = 320, bg = "white")


saveRDS(kappa_vals, "full-model/figures/paper/mcmc_draws/gamma.RDS")
saveRDS(gamma, "full-model/figures/paper/mcmc_draws/gamma.RDS")
saveRDS(theta, "full-model/figures/paper/mcmc_draws/theta.RDS")
saveRDS(lambda, "full-model/figures/paper/mcmc_draws/lambda.RDS")
saveRDS(rand_int, "full-model/figures/paper/mcmc_draws/rand_int.RDS")
saveRDS(phi, "full-model/figures/paper/mcmc_draws/phi.RDS")
saveRDS(pi_prob, "full-model/figures/paper/mcmc_draws/pi_prob.RDS")


all_params <- kappa_vals %>% 
  left_join(rand_int) %>% 
  left_join(lambda) %>% 
  left_join(pi_prob) %>%
  mutate(eta = exp_count(pi, lambda))

returns <- all_params %>% mutate(yr50 = rlevel(50, kappa, sigma, xi, eta)*1000*0.405) # rescale back to 1000s of acres; convert to hectares
returns_summary <- returns %>% group_by(time, region) %>% 
  summarize(med50 = median(yr50[is.finite(yr50)]), 
            lower = quantile(yr50[is.finite(yr50)], probs = 0.025), 
            upper = quantile(yr50[is.finite(yr50)], probs = 0.975)) %>%
  ungroup()

region_key <- readRDS(file = "./full-model/data/processed/region_key.rds")
full_reg_key <- as_tibble(region_key) %>% 
  mutate(region = c(1:84),
         NA_L2CODE = as.factor(NA_L2CODE),
         NA_L1CODE = as.factor(NA_L1CODE),
         NA_L3CODE = as.factor(NA_L3CODE),
         NA_L1NAME = as.factor(str_to_title(NA_L1NAME)))
# reg_cols <- full_reg_key$region
returns_regional <- returns_summary %>% left_join(full_reg_key)

date_seq <- seq(as.Date("1995-01-01"), by = "1 month", length.out = 252) %>% as_tibble() %>% rename(date = value)
time_df <- date_seq %>% mutate(time = 1:252)
returns_regional <- returns_regional %>% left_join(time_df)

p <- returns_regional %>% ggplot() + 
  geom_ribbon(aes(x=date, ymin=lower, ymax=upper, group = region, fill = NA_L1NAME, alpha = 0.5)) +
  geom_line(aes(x=date, y=med50, group = region, alpha = 0.5), linewidth = 0.3, color = 'darkgrey') + scale_y_log10() +
  scale_x_date(name = "Year (1995-2016)", date_breaks = "5 years", date_labels = "%Y") + 
  ylab("Expected burn area (ha)") +
  facet_wrap(. ~ NA_L1NAME, nrow = 2) +
  theme_classic() + theme(legend.position = "none") + 
  theme(strip.text.x = element_text(size = rel(0.9)))

file_name <- "full-model/figures/paper/50yr_returns.png"
ggsave(file_name, p, dpi = 320, bg = "white", width = 17, height = 9)


ecoregions <- read_rds(file = "ecoregions.RDS")
ecoregions_geom <- ecoregions %>% filter(!NA_L2NAME == "UPPER GILA MOUNTAINS (?)") %>% 
  mutate(NA_L2CODE = as.factor(NA_L2CODE),
         NA_L1CODE = as.factor(NA_L1CODE),
         NA_L3CODE = as.factor(NA_L3CODE),
         NA_L1NAME = as.factor(str_to_title(NA_L1NAME)))

er_map_l1 <- ecoregions_geom %>%
  ggplot() +
  geom_sf(size = .2, aes(fill = NA_L1NAME)) + 
  theme_void() +
  coord_sf(ndiscr = FALSE)

er_map_l1 <- ecoregions_geom %>%
  ggplot() +
  geom_sf(size = .2, fill = "white") +
  geom_sf(data = ecoregions_geom,
          aes(fill = NA_L2CODE), alpha = 0.6, lwd = 0, inherit.aes = FALSE, show.legend = FALSE) +
  geom_sf(data = ecoregions_geom %>% group_by(NA_L1CODE) %>% summarise(),
          fill = "transparent", lwd = 1, color = "gray20", inherit.aes = FALSE, show.legend = FALSE) +
  theme_void() +
  coord_sf(ndiscr = FALSE)

ggsave("test_map.png", er_map_l1, dpi = 320, bg = "white")

ri_map <- `sigma_ri_theta-ri_gamma-ri`$ri_matrix %>% as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>% 
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>% 
  separate_wider_delim(cols = "name", delim = ",", names = c("param", "time", "region")) %>% select(-time) %>% distinct() %>%
  mutate(param = as.character(gsub("ri_matrix\\[", "", param)),
         region = as.numeric(gsub("\\]", "", region)),
         param = case_when(param == "1" ~ "sigma",
                           param == "2" ~ "xi",
                           TRUE ~ param))

ri_map <- ri_map %>% group_by(param, region) %>% summarize(med_val = median(value))
ri_map <- ri_map %>% ungroup()
xi_map <- ri_map %>% filter(param == "xi") %>% select(-param) %>% mutate(med_val = exp(med_val))
sigma_map <- ri_map %>% filter(param == "sigma") %>% select(-param) %>% mutate(med_val = exp(med_val))

gamma_map <- gamma %>% group_by(region) %>% summarize(med_gamma = median(gamma)) %>% ungroup() %>% left_join(full_reg_key)
eco_gamma <- ecoregions_geom %>% left_join(gamma_map)
breaks <- classIntervals(c(min(eco_gamma$med_gamma) - .00001, eco_gamma$med_gamma), style = 'fixed', 
                         fixedBreaks = c(-0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4), intervalClosure = 'left')
eco_gamma <- eco_gamma %>% mutate(gamma_cat = cut(med_gamma, unique(breaks$brks)))
p <- ecoregions_geom %>%
  ggplot() +
  geom_sf(size = .1, fill = 'white') +
  geom_sf(data = eco_gamma,
          aes(fill=gamma_cat), alpha = 0.6, lwd = 0, inherit.aes = FALSE) +
  theme_void() + scale_fill_brewer(palette = "Spectral")
ggsave("full-model/figures/paper/gamma_map.png", dpi = 320, bg ='white')

theta_map <- theta %>% group_by(region) %>% summarize(med_theta = median(theta)) %>% ungroup() %>% left_join(full_reg_key)
eco_theta <- ecoregions_geom %>% left_join(theta_map)
breaks <- classIntervals(c(min(eco_theta$med_theta) - .00001, eco_theta$med_theta), style = 'fixed', 
                         fixedBreaks = c(-3, -2, -1, 0, 1, 2), intervalClosure = 'left')
eco_theta <- eco_theta %>% mutate(theta_cat = cut(med_theta, unique(breaks$brks)))
p <- ecoregions_geom %>%
  ggplot() +
  geom_sf(size = .1, fill = 'white') +
  geom_sf(data = eco_theta,
          aes(fill=theta_cat), alpha = 0.6, lwd = 0, inherit.aes = FALSE) +
  theme_void() + scale_fill_brewer(palette = "Spectral")
ggsave("full-model/figures/paper/theta_map.png", dpi = 320, bg ='white')



xi_vals_regions <- xi_map %>% left_join(full_reg_key)
eco_xi <- ecoregions_geom %>% 
  mutate(NA_L2CODE = as.factor(NA_L2CODE), 
         NA_L1CODE = as.factor(NA_L1CODE), 
         NA_L3CODE = as.factor(NA_L3CODE)) %>% 
  left_join(xi_vals_regions)

breaks <- classIntervals(c(min(eco_xi$med_val) - .00001, eco_xi$med_val), style = 'fixed', 
                         fixedBreaks = c(0, 0.2, 0.4, 0.6, 0.8, 2.0), intervalClosure = 'left')

eco_xi_cat <- eco_xi %>% 
  mutate(xi_cat = cut(med_val, unique(breaks$brks)))

p <- ecoregions_geom %>%
  ggplot() +
  geom_sf(size = .1, fill = 'white') +
  geom_sf(data = eco_xi_cat,
          aes(fill=xi_cat), alpha = 0.6, lwd = 0, inherit.aes = FALSE) +
  theme_void() + scale_fill_brewer(palette = 'YlOrRd')
ggsave("full-model/figures/paper/xi_map.png", dpi = 320, bg ='white')


sigma_vals_regions <- sigma_map %>% left_join(full_reg_key)
eco_sigma <- ecoregions_geom %>% 
  mutate(NA_L2CODE = as.factor(NA_L2CODE), 
         NA_L1CODE = as.factor(NA_L1CODE), 
         NA_L3CODE = as.factor(NA_L3CODE)) %>% 
  left_join(sigma_vals_regions)

breaks <- classIntervals(c(min(eco_sigma$med_val) - .00001, eco_sigma$med_val), style = 'fisher', n=5, intervalClosure = 'left')
eco_sigma_cat <- eco_sigma %>% 
  mutate(sigma_cat = cut(med_val, unique(breaks$brks)))

p <- ecoregions_geom %>%
  ggplot() +
  geom_sf(size = .1, fill = 'white') +
  geom_sf(data = eco_sigma_cat,
          aes(fill=sigma_cat), alpha = 0.6, lwd = 0, inherit.aes = FALSE) +
  theme_void() + scale_fill_brewer(palette = 'YlOrRd')
ggsave("full-model/figures/paper/sigma_map.png", dpi = 320, bg ='white')


