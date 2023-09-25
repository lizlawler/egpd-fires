library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(tidyverse)
library(stringr)
library(posterior)
library(lubridate)
library(sf)
library(classInt)
library(RColorBrewer)

## read in region key 
region_key <- readRDS(file = "./full-model/data/processed/region_key.rds")
full_reg_key <- as_tibble(region_key) %>% 
  mutate(region = c(1:84),
         NA_L2CODE = as.factor(NA_L2CODE),
         NA_L1CODE = as.factor(NA_L1CODE),
         NA_L3CODE = as.factor(NA_L3CODE),
         NA_L1NAME = as.factor(str_to_title(NA_L1NAME)))

## read in ecoregion shapes
ecoregions <- read_rds(file = "ecoregions.RDS")
ecoregions_geom <- ecoregions %>% filter(!NA_L2NAME == "UPPER GILA MOUNTAINS (?)") %>% 
  mutate(NA_L2CODE = as.factor(NA_L2CODE),
         NA_L1CODE = as.factor(NA_L1CODE),
         NA_L3CODE = as.factor(NA_L3CODE),
         NA_L1NAME = as.factor(str_to_title(NA_L1NAME)))

# egpd functions
pegpd <- function(y, kappa, sigma, xi) (1 - (1 + xi * (y/sigma))^(-1/xi))^kappa
qegpd <- function(p, kappa, sigma, xi) (sigma/xi) * ( (1 - p^(1/kappa) )^-xi - 1)
regpd <- function(n, kappa, sigma, xi) { # truncated version
  u <- runif(n)
  u_adj <- u * (1 - pegpd(1.001, kappa, sigma, xi)) + pegpd(1.001, kappa, sigma, xi)
  return(qegpd(u_adj, kappa, sigma, xi))
}
exp_egpd <- function(n, kappa, sigma, xi) {
  mcmc_rng <- regpd(n, kappa, sigma, xi) 
  return(mean(mcmc_rng[is.finite(mcmc_rng)]))
}

# expected counts from ZINB params
exp_count <- function(pi, lambda) { # expected counts
  pi_prob <- exp(pi)/(1+exp(pi))
  return((1-pi_prob) * exp(lambda))
}

# probability of 1 or more fire occurring
fire_occur <- function(pi, delta, lambda) {
  pi_prob <- exp(pi)/(1+exp(pi))
  neg_bin <- (delta / (exp(lambda) + delta)) ^ delta
  return(1 - (1-pi_prob)*neg_bin)
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

# cleaning up delta tibble
delta <- delta %>% 
  mutate(name = gsub("delta\\[", "", name),
         name = as.integer(gsub("\\]", "", name))) %>%
  rename(region = name, delta = value)
  
# read in extracted mcmc draws ------
files <- paste0("full-model/figures/paper/mcmc_draws/theta-time_gamma-cst_2000iter/",
                list.files(path = "full-model/figures/paper/mcmc_draws/theta-time_gamma-cst_2000iter/",
                           pattern = ".RDS"))
object_names <- str_remove(basename(files), ".RDS")
for(i in seq_along(object_names)) {
  assign(object_names[i], readRDS(files[i]))
}

## create time series plots of phi values for kappa and lambda -------
date_seq <- seq(as.Date("1990-01-01"), by = "1 month", length.out = 372) %>% as_tibble() %>% rename(date = value)
time_df <- date_seq %>% mutate(time = 1:372)

phi_kappa <- phi %>% select(-lambda) %>% group_by(time, region) %>% summarize(kappa = median(kappa)) %>% ungroup()
p <- phi_kappa %>% left_join(full_reg_key) %>% left_join(time_df) %>%
  ggplot(aes(x = date, y = kappa, color = NA_L2CODE)) + 
  geom_line(linewidth = 0.5) + 
  scale_x_date(name = "Year (1990-2020)", date_breaks = "5 years", date_labels = "%Y") +
  facet_wrap(. ~ NA_L1NAME, nrow = 2) + 
  theme_classic()
ggsave("full-model/figures/paper/phi_kappa_allregs_overtime.pdf", p, width = 15, dpi = 320, bg = "white")

phi_lambda <- phi %>% select(-kappa) %>% group_by(time, region) %>% summarize(lambda = median(lambda)) %>% ungroup()
p <- phi_lambda %>% left_join(full_reg_key) %>% left_join(time_df) %>%
  ggplot(aes(x = date, y = lambda, color = NA_L2CODE)) + 
  geom_line(linewidth = 0.5) + 
  scale_x_date(name = "Year (1990-2020)", date_breaks = "5 years", date_labels = "%Y") +
  facet_wrap(. ~ NA_L1NAME, nrow = 2) + 
  theme_classic()
ggsave("full-model/figures/paper/phi_lambda_allregs_overtime.pdf", p, width = 15, dpi = 320, bg = "white")

phi_med_cal <- phi %>% left_join(full_reg_key) %>% 
  filter(NA_L1NAME == "Mediterranean California") %>% 
  group_by(time, region) %>% 
  summarize(lambda = median(lambda), kappa = median(kappa)) %>% 
  ungroup() %>% left_join(full_reg_key) %>% filter(region == 10)

years <- seq(1, 372, by = 12)[-1]
p <- phi_med_cal %>% left_join(time_df) %>% 
  pivot_longer(cols = c("lambda", "kappa"), names_to = "param", values_to = "value") %>%
  ggplot(aes(x = date, y = value)) + 
  geom_line(linewidth = 0.5, col = "blue") +
  scale_x_date(name = "Year (1990-2020)", date_breaks = "5 years", date_labels = "%Y") + 
  facet_wrap(. ~ param, scales = "free_y") +
  theme_classic() + theme(legend.position = "none")
ggsave("full-model/figures/paper/phi_med_cal_overtime.pdf", p, width = 15, dpi = 320, bg = "white")

p <- phi_med_cal %>%
  ggplot(aes(x = time, y = kappa)) +
  geom_line(linewidth = 0.5, col = "blue") +
  geom_vline(xintercept = years, col = "darkgrey", linewidth = 0.5) +
  theme_classic() + theme(legend.position = "none")
ggsave("full-model/figures/paper/phi_kappa_med_cal_overtime.pdf", p, width = 15, dpi = 320, bg = "white")

p <- phi_med_cal %>%
  ggplot(aes(x = time, y = lambda)) +
  geom_line(linewidth = 0.5, col = "blue") +
  geom_vline(xintercept = years, col = "darkgrey", linewidth = 0.5) +
  theme_classic() + theme(legend.position = "none")
ggsave("full-model/figures/paper/phi_lambda_med_cal_overtime.pdf", p, width = 15, dpi = 320, bg = "white")

# create return level plots ------
all_params <- kappa %>% 
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
  scale_x_date(name = "Year (1990-2020)", date_breaks = "5 years", date_labels = "%Y") + 
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
  scale_x_date(name = "Year (1990-2020)", date_breaks = "5 years", date_labels = "%Y") + 
  ylab("Expected burn area (ha) given fire occurrence") +
  facet_wrap(. ~ NA_L1NAME, nrow = 2) +
  theme_classic() + theme(legend.position = "none") + 
  theme(strip.text.x = element_text(size = rel(0.9)))

file_name <- "full-model/figures/paper/98th_quant_burnareas.pdf"
ggsave(file_name, p, dpi = 320, bg = "white", width = 17, height = 9)

## area-weighted average of burn area exceedances -------
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
  scale_x_date(name = "Year (1990-2020)", date_breaks = "5 years", date_labels = "%Y") + 
  ylab("Expected burn area (ha)") +
  facet_wrap(. ~ NA_L1NAME, nrow = 2) +
  theme_classic() + theme(legend.position = "none") + 
  theme(strip.text.x = element_text(size = rel(0.9)))

file_name <- "full-model/figures/paper/50yr_returns_level1.pdf"
ggsave(file_name, p, dpi = 320, bg = "white", width = 17, height = 9)

## time-averaged kappa, then high quantile of burn area --------
# kappa_temp_avg <- kappa_vals %>% group_by(region, draw) %>% summarize(kappa = mean(kappa)) %>% ungroup()
# all_params_notime <- kappa_temp_avg %>% 
#   left_join(rand_int)
# level98_timeavg <- all_params_notime %>% 
#   mutate(yr50 = high_quant(50, kappa, sigma, xi)*0.405) # convert to hectares (1000s of hectares)
# level98_timeavg_summary <- level98_timeavg %>% group_by(region) %>%
#   summarize(med50 = median(yr50[is.finite(yr50)]), 
#             lower = quantile(yr50[is.finite(yr50)], probs = 0.025), 
#             upper = quantile(yr50[is.finite(yr50)], probs = 0.975)) %>%
#   ungroup()
# level98_timeavg_regions <- level98_timeavg_summary %>% left_join(full_reg_key)

# level95_timeavg <- all_params_notime %>% 
#   mutate(yr20 = high_quant(20, kappa, sigma, xi)*0.405) # convert to hectares (1000s of hectares)
# level95_timeavg_summary <- level95_timeavg %>% group_by(region) %>%
#   summarize(med20 = median(yr20[is.finite(yr20)]), 
#             lower = quantile(yr20[is.finite(yr20)], probs = 0.025), 
#             upper = quantile(yr20[is.finite(yr20)], probs = 0.975)) %>%
#   ungroup()
# level95_timeavg_regions <- level95_timeavg_summary %>% left_join(full_reg_key)
# 
# level75_timeavg <- all_params_notime %>% 
#   mutate(yr4 = high_quant(4, kappa, sigma, xi)*0.405)
# level75_timeavg_summary <- level75_timeavg %>% group_by(region) %>%
#   summarize(med4 = median(yr4[is.finite(yr4)]), 
#             lower = quantile(yr4[is.finite(yr4)], probs = 0.025), 
#             upper = quantile(yr4[is.finite(yr4)], probs = 0.975)) %>%
#   ungroup()
# level75_timeavg_regions <- level75_timeavg_summary %>% left_join(full_reg_key)
# 
# level50_timeavg <- all_params_notime %>% 
#   mutate(yr2 = high_quant(2, kappa, sigma, xi)*0.405) # rescale back to 1000s of acres; convert to hectares
# level50_timeavg_summary <- level50_timeavg %>% group_by(region) %>%
#   summarize(med2 = median(yr2[is.finite(yr2)]), 
#             lower = quantile(yr2[is.finite(yr2)], probs = 0.025), 
#             upper = quantile(yr2[is.finite(yr2)], probs = 0.975)) %>%
#   ungroup()
# level50_timeavg_regions <- level50_timeavg_summary %>% left_join(full_reg_key)

# ## create maps for each quantile level, using same interval break categories -------
# allquantlevels <- c(level50_timeavg_summary$med2, level75_timeavg_summary$med4, level95_timeavg_summary$med20, level98_timeavg_summary$med50)
# breaks <- classIntervals(c(min(allquantlevels) - .00001, allquantlevels), style = 'quantile',intervalClosure = 'left')
# 
# eco_burns98 <- ecoregions_geom %>% 
#   left_join(level98_timeavg_regions) %>% 
#   mutate(burns_cat = cut(med50, unique(breaks$brks)))
# p <- ecoregions_geom %>%
#   ggplot() +
#   geom_sf(size = .1, fill = 'white') +
#   geom_sf(data = eco_burns98,
#           aes(fill=burns_cat), alpha = 0.6, lwd = 0, inherit.aes = FALSE) + 
#   scale_fill_brewer(palette = 'YlOrRd') +
#   theme_void() 
# ggsave("full-model/figures/paper/onekappa_98th_quant_map.pdf", dpi = 320, bg ='white')

# eco_burns95 <- ecoregions_geom %>% 
#   left_join(level95_timeavg_regions) %>% 
#   mutate(burns_cat = cut(med20, unique(breaks$brks)))
# p <- ecoregions_geom %>%
#   ggplot() +
#   geom_sf(size = .1, fill = 'white') +
#   geom_sf(data = eco_burns95,
#           aes(fill=burns_cat), alpha = 0.6, lwd = 0, inherit.aes = FALSE) + 
#   scale_fill_brewer(palette = 'YlOrRd') +
#   theme_void() 
# ggsave("full-model/figures/paper/onekappa_95th_quant_map.pdf", dpi = 320, bg ='white')
# 
# eco_burns75 <- ecoregions_geom %>% 
#   left_join(level75_timeavg_regions) %>% 
#   mutate(burns_cat = cut(med4, unique(breaks$brks)))
# p <- ecoregions_geom %>%
#   ggplot() +
#   geom_sf(size = .1, fill = 'white') +
#   geom_sf(data = eco_burns75,
#           aes(fill=burns_cat), alpha = 0.6, lwd = 0, inherit.aes = FALSE) + 
#   scale_fill_brewer(palette = 'YlOrRd') +
#   theme_void() 
# ggsave("full-model/figures/paper/onekappa_75th_quant_map.pdf", dpi = 320, bg ='white')
# 
# eco_burns50 <- ecoregions_geom %>% 
#   left_join(level50_timeavg_regions) %>% 
#   mutate(burns_cat = cut(med2, unique(breaks$brks)))
# p <- ecoregions_geom %>%
#   ggplot() +
#   geom_sf(size = .1, fill = 'white') +
#   geom_sf(data = eco_burns50,
#           aes(fill=burns_cat), alpha = 0.6, lwd = 0, inherit.aes = FALSE) + 
#   scale_fill_brewer(palette = 'YlOrRd') +
#   theme_void() 
# ggsave("full-model/figures/paper/onekappa_50th_quant_map.pdf", dpi = 320, bg ='white')

# ## create maps of each parameter (kappa, sigma, xi) using 'quantile' break schema ---------
# all_params_notime <- all_params_notime %>% 
#   group_by(region) %>% 
#   summarize(kappa = mean(kappa),
#             sigma = mean(sigma),
#             xi = mean(xi)) %>%
#   left_join(full_reg_key)
# kappa_breaks <- classIntervals(c(min(all_params_notime$kappa) - .00001, all_params_notime$kappa), style = 'quantile', intervalClosure = 'left')
# sigma_breaks <- classIntervals(c(min(all_params_notime$sigma) - .00001, all_params_notime$sigma), style = 'quantile', intervalClosure = 'left')
# xi_breaks <- classIntervals(c(min(all_params_notime$xi) - .00001, all_params_notime$xi), style = 'quantile',intervalClosure = 'left')
# all_params_notime_eco <- ecoregions_geom %>% left_join(all_params_notime) %>%
#   mutate(kappa_cat = cut(kappa, unique(kappa_breaks$brks)),
#          sigma_cat = cut(sigma, unique(sigma_breaks$brks)),
#          xi_cat = cut(xi, unique(xi_breaks$brks)))
# 
# p <- ecoregions_geom %>%
#   ggplot() +
#   geom_sf(size = .1, fill = 'white') +
#   geom_sf(data = all_params_notime_eco,
#           aes(fill=kappa_cat), alpha = 0.6, lwd = 0, inherit.aes = FALSE) + 
#   scale_fill_brewer(palette = 'YlOrRd') +
#   theme_void() 
# ggsave("full-model/figures/paper/kappa_timeavg_map.pdf", dpi = 320, bg ='white')
# 
# p <- ecoregions_geom %>%
#   ggplot() +
#   geom_sf(size = .1, fill = 'white') +
#   geom_sf(data = all_params_notime_eco,
#           aes(fill=sigma_cat), alpha = 0.6, lwd = 0, inherit.aes = FALSE) + 
#   scale_fill_brewer(palette = 'YlOrRd') +
#   theme_void()
# ggsave("full-model/figures/paper/sigma_map.pdf", dpi = 320, bg ='white')
# 
# p <- ecoregions_geom %>%
#   ggplot() +
#   geom_sf(size = .1, fill = 'white') +
#   geom_sf(data = all_params_notime_eco,
#           aes(fill=xi_cat), alpha = 0.6, lwd = 0, inherit.aes = FALSE) + 
#   scale_fill_brewer(palette = 'YlOrRd') +
#   theme_void()
# ggsave("full-model/figures/paper/xi_cat.pdf", dpi = 320, bg ='white')


## create boxplot of predicted counts for entire US annually for all timepoints (holdout and training) and overlay truth ------
count_preds <- lambda %>% 
  left_join(pi_prob) %>% 
  mutate(preds = exp_count(pi, lambda)) %>% 
  select(c(draw, time, region, preds)) %>%
  left_join(full_reg_key) %>% 
  left_join(time_df) %>% 
  mutate(year = year(date))

all_years <- 1990:2020
first_five <- 1990:1994
last_five <- 2020:2016
test_years <- sort(c(first_five, last_five))
train_years <- setdiff(all_years, test_years)

preds_annual <- count_preds %>% 
  group_by(NA_L1NAME, year, draw) %>% 
  summarize(total_fires = sum(preds)) %>% 
  ungroup() %>% 
  mutate(train = case_when(year %in% train_years ~ TRUE,
                           year %in% test_years ~ FALSE))

# read in train counts
stan_data <- readRDS("full-model/data/stan_data_joint.RDS")
y_train_count <- stan_data$y_train_count %>% 
  as_tibble() %>% 
  rowid_to_column(var = "time") %>% 
  pivot_longer(!time, names_to = "region", values_to = "value") %>%
  mutate(region = as.numeric(gsub("V", "", region)),
         time = time + 60) %>%
  left_join(full_reg_key) %>%
  left_join(time_df) %>% 
  mutate(year = year(date)) %>%
  group_by(NA_L1NAME, year) %>%
  summarize(true_count = sum(value)) %>%
  ungroup()

# read in holdout counts
y_hold_count <- stan_data$y_hold_count %>% 
  as_tibble() %>% 
  rowid_to_column(var = "time") %>% 
  pivot_longer(!time, names_to = "region", values_to = "value") %>%
  mutate(region = as.numeric(gsub("V", "", region)),
         time = case_when(time > 60 ~ time + 252,
                          TRUE ~ time)) %>%
  left_join(full_reg_key) %>%
  left_join(time_df) %>% 
  mutate(year = year(date)) %>%
  group_by(NA_L1NAME, year) %>%
  summarize(true_count = sum(value)) %>%
  ungroup()

true_counts <- bind_rows(y_train_count, y_hold_count) %>% arrange(year, NA_L1NAME, locale = ".en")

preds_winsor_limits <- preds_annual %>% group_by(NA_L1NAME, year) %>% 
  summarize(lower = quantile(total_fires, 0.05),
            upper = quantile(total_fires, 0.95)) %>%
  ungroup()
preds_winsor <- preds_annual %>% 
  left_join(preds_winsor_limits) %>%
  mutate(winsor_total = case_when(total_fires <= lower ~ lower,
                                  total_fires >= upper ~ upper,
                                  .default = total_fires))

p <- preds_annual %>% 
  ggplot(aes(x = year, y = total_fires, group = year, color = train)) + 
  geom_boxplot(outlier.size = 0.2) + scale_color_grey(start = 0.4, end = 0.6) +
  geom_point(inherit.aes = FALSE, data = true_counts, aes(x = year, y = true_count), col = "red", size = 0.35) +
  geom_line(inherit.aes = FALSE, data = true_counts, aes(x = year, y = true_count), col = "red", linewidth = 0.35) +
  facet_wrap(. ~ NA_L1NAME, scales = "free_y", nrow = 2) + 
  theme_classic() + theme(legend.position = "none")
ggsave("full-model/figures/paper/counts_preds-vs-truth_bylevel1.pdf", width = 15)

p <- preds_winsor %>% 
  ggplot(aes(x = year, y = winsor_total, group = year, color = train)) + 
  geom_boxplot(outlier.size = 0.2) + scale_color_grey(start = 0.4, end = 0.6) +
  geom_point(inherit.aes = FALSE, data = true_counts, aes(x = year, y = true_count), col = "red", size = 0.35) +
  geom_line(inherit.aes = FALSE, data = true_counts, aes(x = year, y = true_count), col = "red", linewidth = 0.35) +
  facet_wrap(. ~ NA_L1NAME, scales = "free_y", nrow = 2) + 
  theme_classic() + theme(legend.position = "none")
ggsave("full-model/figures/paper/counts_preds-vs-truth_bylevel1_winsor.pdf", width = 15)

# over entire US (instead of by level 1)
true_count_entireUS <- true_counts %>% group_by(year) %>% summarize(true_count = sum(true_count)) %>% ungroup()
preds_annual_US <- preds_annual %>% group_by(year, draw, train) %>% summarize(total_fires = sum(total_fires)) %>% ungroup()
p <- preds_annual_US %>%
  ggplot(aes(x = year, y = total_fires, group = year, color = train)) + 
  geom_boxplot(outlier.size = 0.2) + scale_color_grey(start = 0.4, end = 0.6) +
  geom_point(inherit.aes = FALSE, data = true_count_entireUS, aes(x = year, y = true_count), col = "red", size = 0.35) +
  geom_line(inherit.aes = FALSE, data = true_count_entireUS, aes(x = year, y = true_count), col = "red", linewidth = 0.35) +
  theme_classic() + theme(legend.position = "none")
ggsave("full-model/figures/paper/counts_preds-vs-truth_entireUS.pdf", width = 15)

preds_winsor_limits_US <- preds_annual_US %>% group_by(year) %>%
  summarize(lower = quantile(total_fires, 0.05),
            upper = quantile(total_fires, 0.95)) %>%
  ungroup()

preds_winsor_US <- preds_annual_US %>% 
  left_join(preds_winsor_limits_US) %>%
  mutate(winsor_total = case_when(total_fires <= lower ~ lower,
                                  total_fires >= upper ~ upper,
                                  .default = total_fires))
p <- preds_winsor_US %>%
  ggplot(aes(x = year, y = winsor_total, group = year, color = train)) + 
  geom_boxplot(outlier.size = 0.2) + scale_color_grey(start = 0.4, end = 0.6) +
  geom_point(inherit.aes = FALSE, data = true_count_entireUS, aes(x = year, y = true_count), col = "red", size = 0.35) +
  geom_line(inherit.aes = FALSE, data = true_count_entireUS, aes(x = year, y = true_count), col = "red", linewidth = 0.35) +
  theme_classic() + theme(legend.position = "none")
ggsave("full-model/figures/paper/counts_preds-vs-truth_entireUS_winsor.pdf", width = 15)

## create boxplot of predicted burn area for entire US annually for all timepoints (holdout and training) and overlay truth ------
burn_params <- kappa %>% 
  left_join(rand_int)

burn_preds <- readRDS("full-model/figures/paper/burn_preds_df.RDS")
burn_preds_negbin <- burn_preds %>% left_join(delta) %>% left_join(lambda) %>% left_join(pi_prob)
burn_preds_negbin <- burn_preds_negbin %>% select(-NA_L3NAME)

test_preds <- burn_preds_negbin[50:65,]
test_preds_negbin <- test_preds %>% mutate(true_preds = preds * fire_occur(pi, delta, lambda))

burn_preds_negbin <- burn_preds_negbin %>% mutate(condl_preds = preds * fire_occur(pi, delta, lambda))
burn_preds <- burn_preds_negbin %>% select(-c("preds", "kappa", "sigma", "xi", "delta", "lambda", "pi")) %>%
  rename(preds = condl_preds)
# burn_preds <- burn_params %>%
#   mutate(preds = purrr::pmap_dbl(list(500, kappa, sigma, xi), med_egpd)) %>%
#   left_join(full_reg_key) %>%
#   left_join(time_df) %>%
#   mutate(year = year(date))
# saveRDS(burn_preds, file = "full-model/figures/paper/burn_preds_df.RDS")

burn_preds_annual <- burn_preds %>% 
  group_by(NA_L1NAME, draw, year) %>% 
  summarize(total_area = sum(preds[is.finite(preds)])) %>% 
  ungroup() %>% 
  mutate(train = case_when(year %in% train_years ~ TRUE,
                           year %in% test_years ~ FALSE))

true_burns <- readRDS("full-model/data/burn_area_level1.RDS") %>% 
  mutate(NA_L1CODE = as.factor(NA_L1CODE), total_burns = total_burns) %>%
  rename(year = fire_yr, true_area = total_burns)

true_burns_full <- true_counts %>% 
  left_join(full_reg_key %>% select(c("NA_L1CODE", "NA_L1NAME")) %>% distinct()) %>% 
  left_join(true_burns) %>%
  select(-c("true_count", "NA_L1CODE")) %>%
  mutate(true_area = true_area/1000)

burn_preds_winsor_limits <- burn_preds_annual %>% group_by(NA_L1NAME, year) %>% 
  summarize(lower = quantile(total_area, 0.05, na.rm = TRUE),
            upper = quantile(total_area, 0.95, na.rm = TRUE)) %>%
  ungroup()
burn_preds_winsor <- burn_preds_annual %>% 
  left_join(burn_preds_winsor_limits) %>%
  mutate(winsor_total = case_when(total_area <= lower ~ lower,
                                  total_area >= upper ~ upper,
                                  .default = total_area))

p <- burn_preds_annual %>% 
  ggplot(aes(x = year, y = total_area, group = year, color = train)) + 
  geom_boxplot(outlier.size = 0.2) + scale_color_grey(start = 0.4, end = 0.6) +
  geom_point(inherit.aes = FALSE, data = true_burns_full, aes(x = year, y = true_area), col = "red", size = 0.35) +
  geom_line(inherit.aes = FALSE, data = true_burns_full, aes(x = year, y = true_area), col = "red", linewidth = 0.35) +
  facet_wrap(. ~ NA_L1NAME, scales = "free_y", nrow = 2) + scale_y_log10() +
  theme_classic() + theme(legend.position = "none")
ggsave("full-model/figures/paper/burns_preds-vs-truth_bylevel1_v2.pdf", width = 15)

p <- burn_preds_winsor %>% 
  ggplot(aes(x = year, y = winsor_total, group = year, color = train)) + 
  geom_boxplot(outlier.size = 0.2) + scale_color_grey(start = 0.4, end = 0.6) +
  geom_point(inherit.aes = FALSE, data = true_burns_full, aes(x = year, y = true_area), col = "red", size = 0.35) +
  geom_line(inherit.aes = FALSE, data = true_burns_full, aes(x = year, y = true_area), col = "red", linewidth = 0.35) +
  facet_wrap(. ~ NA_L1NAME, scales = "free_y", nrow = 2) + scale_y_log10() +
  theme_classic() + theme(legend.position = "none")
ggsave("full-model/figures/paper/burns_preds-vs-truth_bylevel1_winsor_v2.pdf", width = 15)

# over entire US (instead of by level 1)
true_count_entireUS <- true_counts %>% group_by(year) %>% summarize(true_count = sum(true_count)) %>% ungroup()
preds_annual_US <- preds_annual %>% group_by(year, draw, train) %>% summarize(total_fires = sum(total_fires)) %>% ungroup()
p <- preds_annual_US %>%
  ggplot(aes(x = year, y = total_fires, group = year, color = train)) + 
  geom_boxplot(outlier.size = 0.2) + scale_color_grey(start = 0.4, end = 0.6) +
  geom_point(inherit.aes = FALSE, data = true_count_entireUS, aes(x = year, y = true_count), col = "red", size = 0.35) +
  geom_line(inherit.aes = FALSE, data = true_count_entireUS, aes(x = year, y = true_count), col = "red", linewidth = 0.35) +
  theme_classic() + theme(legend.position = "none")
ggsave("full-model/figures/paper/counts_preds-vs-truth_entireUS.pdf", width = 15)

preds_winsor_limits_US <- preds_annual_US %>% group_by(year) %>%
  summarize(lower = quantile(total_fires, 0.05),
            upper = quantile(total_fires, 0.95)) %>%
  ungroup()

preds_winsor_US <- preds_annual_US %>% 
  left_join(preds_winsor_limits_US) %>%
  mutate(winsor_total = case_when(total_fires <= lower ~ lower,
                                  total_fires >= upper ~ upper,
                                  .default = total_fires))
p <- preds_winsor_US %>%
  ggplot(aes(x = year, y = winsor_total, group = year, color = train)) + 
  geom_boxplot(outlier.size = 0.2) + scale_color_grey(start = 0.4, end = 0.6) +
  geom_point(inherit.aes = FALSE, data = true_count_entireUS, aes(x = year, y = true_count), col = "red", size = 0.35) +
  geom_line(inherit.aes = FALSE, data = true_count_entireUS, aes(x = year, y = true_count), col = "red", linewidth = 0.35) +
  theme_classic() + theme(legend.position = "none")
ggsave("full-model/figures/paper/counts_preds-vs-truth_entireUS_winsor.pdf", width = 15)


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
