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

all_years <- 1990:2020
first_five <- 1990:1994
last_five <- 2020:2016
test_years <- sort(c(first_five, last_five))
train_years <- setdiff(all_years, test_years)

date_seq <- seq(as.Date("1990-01-01"), by = "1 month", length.out = 372) %>% as_tibble() %>% rename(date = value)
time_df <- date_seq %>% mutate(time = 1:372)

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
# 
# # cleaning up delta tibble
# delta <- delta %>% 
#   mutate(name = gsub("delta\\[", "", name),
#          name = as.integer(gsub("\\]", "", name))) %>%
#   rename(region = name, delta = value)
# 
# saveRDS(delta, "full-model/figures/paper/mcmc_draws/theta-time_gamma-ri_1000iter/delta.RDS")
  
# read in extracted mcmc draws ------
files <- paste0("full-model/figures/paper/mcmc_draws/theta-time_gamma-ri_1000iter/",
                list.files(path = "full-model/figures/paper/mcmc_draws/theta-time_gamma-ri_1000iter/",
                           pattern = ".RDS"))
files <- files[!grepl("old", files)]
files <- files[!grepl("beta", files)]
files <- files[!grepl("phi", files)]
files <- files[!grepl("kappa", files)]
object_names <- str_remove(basename(files), ".RDS")
for(i in seq_along(object_names)) {
  assign(object_names[i], readRDS(files[i]))
}

## create time series plots of phi values for kappa and lambda ------
# phi_kappa <- phi %>% select(-lambda) %>% group_by(time, region) %>% summarize(kappa = median(kappa)) %>% ungroup()
# p <- phi_kappa %>% left_join(full_reg_key) %>% left_join(time_df) %>%
#   ggplot(aes(x = date, y = kappa, color = NA_L2CODE)) + 
#   geom_line(linewidth = 0.5) + 
#   scale_x_date(name = "Year (1990-2020)", date_breaks = "5 years", date_labels = "%Y") +
#   facet_wrap(. ~ NA_L1NAME, nrow = 2) + 
#   theme_classic()
# ggsave("full-model/figures/paper/phi_kappa_allregs_overtime.pdf", p, width = 15, dpi = 320, bg = "white")
# 
# ## create time series for theta ------
# years <- seq(1, 372, by = 12)[-1]
# theta_time <- theta %>% group_by(timepoint) %>% summarize(theta = median(theta)) %>% ungroup() %>% rename(time = timepoint)
# p <- theta_time %>% left_join(time_df) %>%
#   ggplot(aes(x = time, y = theta)) + 
#   geom_line(linewidth = 0.5, col = "blue") + 
#   geom_vline(xintercept = years, col = "darkgrey", linewidth = 0.5) +
#   # scale_x_date(name = "Year (1990-2020)", date_breaks = "5 years", date_labels = "%Y") +
#   theme_classic()
# ggsave("full-model/figures/paper/theta_overtime.pdf", p, width = 15, dpi = 320, bg = "white")

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
saveRDS(returns_regional, "full-model/figures/paper/tibbles/returns_regional_gamma-cst.RDS")

p <- returns_regional %>% ggplot() + 
  geom_ribbon(aes(x=date, ymin=lower, ymax=upper, group = region, fill = NA_L1NAME, alpha = 0.5)) +
  geom_line(aes(x=date, y=med50, group = region, alpha = 0.5), linewidth = 0.5, color = 'darkgrey') + scale_y_log10() +
  scale_x_date(name = "Year (1990-2020)", date_breaks = "5 years", date_labels = "%Y") + 
  # scale_fill_brewer(palette = 'Set3') +
  ylab("Expected burn area (ha)") +
  facet_wrap(. ~ NA_L1NAME, nrow = 2) +
  theme_classic() + theme(legend.position = "none") + 
  theme(strip.text.x = element_text(size = rel(1.1)))

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
saveRDS(level98_areas_regional, "full-model/figures/paper/tibbles/level98_areas_regional_gamma-cst.RDS")

p <- level98_areas_regional %>% ggplot() + 
  geom_ribbon(aes(x=date, ymin=lower, ymax=upper, group = region, fill = NA_L1NAME, alpha = 0.5)) +
  geom_line(aes(x=date, y=med50, group = region), linewidth = 0.5, color = 'darkgrey') + scale_y_log10() +
  scale_x_date(name = "Year (1990-2020)", date_breaks = "5 years", date_labels = "%Y") + 
  ylab("Expected burn area (ha) given fire occurrence") +
  facet_wrap(. ~ NA_L1NAME, nrow = 2) +
  theme_classic() + theme(legend.position = "none") + 
  theme(strip.text.x = element_text(size = rel(1.1)))

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
  mutate(year = year(date)) %>% 
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

preds_winsor_limits <- count_preds %>% group_by(NA_L1NAME, year) %>% 
  summarize(lower = quantile(total_fires, 0.025),
            upper = quantile(total_fires, 0.975)) %>%
  ungroup()
count_preds_winsor <- count_preds %>% 
  left_join(preds_winsor_limits) %>%
  mutate(winsor_total = case_when(total_fires <= lower ~ lower,
                                  total_fires >= upper ~ upper,
                                  .default = total_fires))

p <- count_preds %>% 
  ggplot(aes(x = year, y = total_fires, group = year, color = train)) + 
  geom_boxplot(outlier.size = 0.2) + 
  scale_color_grey(start = 0.4, end = 0.6) +
  geom_point(inherit.aes = FALSE, data = true_counts, aes(x = year, y = true_count), col = "red", size = 0.35) +
  geom_line(inherit.aes = FALSE, data = true_counts, aes(x = year, y = true_count), col = "red", linewidth = 0.35) +
  facet_wrap(. ~ NA_L1NAME, scales = "free_y", nrow = 2) + 
  xlab("Year (1990-2020)") + 
  ylab("Expected number of fires") +
  theme_classic() + 
  theme(legend.position = "none")
ggsave("full-model/figures/paper/counts_preds-vs-truth_gamma-ri.pdf", width = 15)

p <- count_preds_winsor %>% 
  ggplot(aes(x = year, y = winsor_total, group = year, color = train)) + 
  geom_boxplot(outlier.size = 0.2) + 
  scale_color_grey(start = 0.4, end = 0.6) +
  geom_point(inherit.aes = FALSE, data = true_counts, aes(x = year, y = true_count), col = "red", size = 0.35) +
  geom_line(inherit.aes = FALSE, data = true_counts, aes(x = year, y = true_count), col = "red", linewidth = 0.35) +
  facet_wrap(. ~ NA_L1NAME, scales = "free_y", nrow = 2) + 
  xlab("Year (1990-2020)") + 
  ylab("Expected number of fires") +
  theme_classic() + 
  theme(legend.position = "none")
ggsave("full-model/figures/paper/counts_preds-vs-truth_winsor_gamma-ri.pdf", width = 15)


## create boxplot of predicted burn area for entire US annually for all timepoints (holdout and training) and overlay truth ------
burn_preds_gamma_ri <- readRDS("full-model/figures/paper/mcmc_draws/theta-time_gamma-ri_1000iter/burn_pred.RDS")
burn_preds_gamma_cst <- readRDS("full-model/figures/paper/mcmc_draws/theta-time_gamma-cst_1000iter/burn_pred.RDS")

# first create for 'gamma-cst' model
burn_preds_gamma_cst <- burn_preds_gamma_cst %>% 
  left_join(full_reg_key) %>% 
  left_join(time_df) %>% 
  mutate(year = year(date)) %>% 
  rename(preds = value) %>%
  group_by(NA_L1NAME, year, draw) %>% 
  summarize(total_area = sum(preds[is.finite(preds)])*1000*0.405) %>% 
  ungroup() %>% 
  mutate(train = case_when(year %in% train_years ~ TRUE,
                           year %in% test_years ~ FALSE))

# create for 'gamma-ri' model
burn_preds_gamma_ri <- burn_preds_gamma_ri %>% 
  left_join(full_reg_key) %>% 
  left_join(time_df) %>% 
  mutate(year = year(date)) %>% rename(preds = value) %>%
  group_by(NA_L1NAME, year, draw) %>% 
  summarize(total_area = sum(preds[is.finite(preds)])*1000*0.405) %>% 
  ungroup() %>% 
  mutate(train = case_when(year %in% train_years ~ TRUE,
                           year %in% test_years ~ FALSE))

burn_preds_both <- burn_preds_gamma_cst %>% 
  mutate(model = 'cst') %>% 
  rbind(burn_preds_gamma_ri %>% mutate(model = 'ri')) %>%
  mutate(model = as.factor(model))
# G1 burn preds ------
# g1_burns_full <- g1_burn_preds_1 %>% 
#   rbind(g1_burn_preds_2 %>% mutate(draw = 1000 + draw)) %>%
#   rbind(g1_burn_preds_3 %>% mutate(draw = 2000 + draw))
# 
# g1_burns_full <- g1_burns_full %>% left_join(full_reg_key) %>% left_join(time_df) %>% mutate(year = year(date))
# g1_burns_full_annual <- g1_burns_full %>%
#   group_by(NA_L1NAME, year, draw) %>% 
#   summarize(total_area = sum(preds[is.finite(preds)])*1000*0.405) %>% 
#   ungroup() %>% 
#   mutate(train = case_when(year %in% train_years ~ TRUE,
#                            year %in% test_years ~ FALSE))

# read in true burn areas -------
true_burns <- readRDS("full-model/data/burn_df_agg.RDS") %>% 
  mutate(NA_L2NAME = case_when(NA_L2NAME == "UPPER GILA MOUNTAINS (?)" ~ "UPPER GILA MOUNTAINS",
                               TRUE ~ NA_L2NAME),
         NA_L1NAME = as.factor(str_to_title(NA_L1NAME)))
true_burns_level1 <- true_burns %>% group_by(NA_L1NAME, fire_yr) %>%
  summarize(true_area = sum(total_burns)*0.405) %>% ungroup()

true_burns_level1_full <- true_burns_level1 %>% complete(NA_L1NAME, fire_yr) %>% rename(year = fire_yr)

# burn_preds_both_winsor_limits <- burn_preds_both %>% group_by(NA_L1NAME, model, year) %>% 
#   summarize(lower = quantile(total_area, 0.05, na.rm = TRUE),
#             upper = quantile(total_area, 0.95, na.rm = TRUE)) %>%
#   ungroup()
# burn_preds_both_winsor <- burn_preds_both %>%
#   left_join(burn_preds_both_winsor_limits) %>%
#   mutate(winsor_total = case_when(total_area <= lower ~ lower,
#                                   total_area >= upper ~ upper,
#                                   .default = total_area))

# g1_burn_preds_winsor_limits <- g1_burns_full_annual %>% group_by(NA_L1NAME, year) %>%
#   summarize(lower = quantile(total_area, 0.05, na.rm = TRUE),
#             upper = quantile(total_area, 0.95, na.rm = TRUE)) %>%
#   ungroup()
# g1_burn_preds_winsor <- g1_burns_full_annual %>%
#   left_join(g1_burn_preds_winsor_limits) %>%
#   mutate(winsor_total = case_when(total_area <= lower ~ lower,
#                                   total_area >= upper ~ upper,
#                                   .default = total_area))

# p <- g1_burns_full_annual %>% 
#   ggplot(aes(x = year, y = total_area, group = year, color = train)) + 
#   geom_boxplot(outlier.size = 0.2) + scale_color_grey(start = 0.4, end = 0.6) +
#   geom_point(inherit.aes = FALSE, data = true_burns_level1_full, aes(x = year, y = true_area), col = "red", size = 0.35) +
#   geom_line(inherit.aes = FALSE, data = true_burns_level1_full, aes(x = year, y = true_area), col = "red", linewidth = 0.35) +
#   facet_wrap(. ~ NA_L1NAME, scales = "free_y", nrow = 2) +
#   scale_y_log10() +
#   xlab("Year (1990-2020)") + 
#   ylab("Expected burn area (ha)") +
#   theme_classic() + theme(legend.position = "none")
# 
# p <- g1_burn_preds_winsor %>%
#   ggplot(aes(x = year, y = winsor_total, group = year, color = train)) +
#   geom_boxplot(outlier.size = 0.2) + scale_color_grey(start = 0.4, end = 0.6) +
#   geom_point(inherit.aes = FALSE, data = true_burns_level1_full, aes(x = year, y = true_area), col = "red", size = 0.35) +
#   geom_line(inherit.aes = FALSE, data = true_burns_level1_full, aes(x = year, y = true_area), col = "red", linewidth = 0.35) +
#   facet_wrap(. ~ NA_L1NAME, scales = "free_y", nrow = 2) +
#   scale_y_log10() +
#   xlab("Year (1990-2020)") +
#   ylab("Expected burn area (ha)") +
#   theme_classic() + theme(legend.position = "none")
# ggsave("full-model/figures/paper/burns_preds-vs-truth_winsor_g1.pdf", width = 15)

burn_preds_gamma_ri <- readRDS("full-model/figures/paper/mcmc_draws/theta-time_gamma-ri_1000iter/burn_pred_old.RDS")
burn_preds_gamma_cst <- readRDS("full-model/figures/paper/mcmc_draws/theta-time_gamma-cst_1000iter/burn_pred_old.RDS")

# first create for 'gamma-cst' model
burn_preds_gamma_cst <- burn_preds_gamma_cst %>% 
  left_join(full_reg_key) %>% 
  left_join(time_df) %>% 
  mutate(year = year(date)) %>% 
  rename(preds = value) %>%
  group_by(NA_L1NAME, year, draw) %>% 
  summarize(total_area = sum(preds[is.finite(preds)])*1000*0.405) %>% 
  ungroup() %>% 
  mutate(train = case_when(year %in% train_years ~ TRUE,
                           year %in% test_years ~ FALSE))

burn_preds_gamma_cst_winsor_limits <- burn_preds_gamma_cst %>% group_by(NA_L1NAME, year) %>%
  summarize(lower = quantile(total_area, 0.05, na.rm = TRUE),
            upper = quantile(total_area, 0.95, na.rm = TRUE)) %>%
  ungroup()
burn_preds_gamma_cst_winsor <- burn_preds_gamma_cst %>%
  left_join(burn_preds_gamma_cst_winsor_limits) %>%
  mutate(winsor_total = case_when(total_area <= lower ~ lower,
                                  total_area >= upper ~ upper,
                                  .default = total_area))

# create for 'gamma-ri' model
burn_preds_gamma_ri <- burn_preds_gamma_ri %>% 
  left_join(full_reg_key) %>% 
  left_join(time_df) %>% 
  mutate(year = year(date)) %>% rename(preds = value) %>%
  group_by(NA_L1NAME, year, draw) %>% 
  summarize(total_area = sum(preds[is.finite(preds)])*1000*0.405) %>% 
  ungroup() %>% 
  mutate(train = case_when(year %in% train_years ~ TRUE,
                           year %in% test_years ~ FALSE))

burn_preds_gamma_ri_winsor_limits <- burn_preds_gamma_ri %>% group_by(NA_L1NAME, year) %>%
  summarize(lower = quantile(total_area, 0.05, na.rm = TRUE),
            upper = quantile(total_area, 0.95, na.rm = TRUE)) %>%
  ungroup()
burn_preds_gamma_ri_winsor <- burn_preds_gamma_ri %>%
  left_join(burn_preds_gamma_ri_winsor_limits) %>%
  mutate(winsor_total = case_when(total_area <= lower ~ lower,
                                  total_area >= upper ~ upper,
                                  .default = total_area))

p <- burn_preds_gamma_cst %>% 
  ggplot(aes(x = year, y = total_area, group = year, color = train)) + 
  geom_boxplot(outlier.size = 0.2) + 
  scale_color_grey(start = 0.4, end = 0.6) +
  geom_point(inherit.aes = FALSE, data = true_burns_level1_full, aes(x = year, y = true_area), col = "red", size = 0.35) +
  geom_line(inherit.aes = FALSE, data = true_burns_level1_full, aes(x = year, y = true_area), col = "red", linewidth = 0.35) +
  facet_wrap(. ~ NA_L1NAME, scales = "free_y", nrow = 2) +
  scale_y_log10() +
  xlab("Year (1990-2020)") + 
  ylab("Expected burn area (ha)") +
  theme_classic() + theme(legend.position = "none")
ggsave("full-model/figures/paper/burns_preds-vs-truth_gamma-cst_old.pdf", width = 15)

p <- burn_preds_gamma_cst_winsor %>% 
  ggplot(aes(x = year, y = winsor_total, group = year, color = train)) + 
  geom_boxplot(outlier.size = 0.2) + 
  scale_color_grey(start = 0.4, end = 0.6) +
  geom_point(inherit.aes = FALSE, data = true_burns_level1_full, aes(x = year, y = true_area), col = "red", size = 0.35) +
  geom_line(inherit.aes = FALSE, data = true_burns_level1_full, aes(x = year, y = true_area), col = "red", linewidth = 0.35) +
  facet_wrap(. ~ NA_L1NAME, scales = "free_y", nrow = 2) +
  scale_y_log10() +
  xlab("Year (1990-2020)") + 
  ylab("Expected burn area (ha)") +
  theme_classic() + theme(legend.position = "none")
ggsave("full-model/figures/paper/burns_preds-vs-truth_winsor_gamma-cst_old.pdf", width = 15)

p <- burn_preds_gamma_ri %>% 
  ggplot(aes(x = year, y = total_area, group = year, color = train)) + 
  geom_boxplot(outlier.size = 0.2) + 
  scale_color_grey(start = 0.4, end = 0.6) +
  geom_point(inherit.aes = FALSE, data = true_burns_level1_full, aes(x = year, y = true_area), col = "red", size = 0.35) +
  geom_line(inherit.aes = FALSE, data = true_burns_level1_full, aes(x = year, y = true_area), col = "red", linewidth = 0.35) +
  facet_wrap(. ~ NA_L1NAME, scales = "free_y", nrow = 2) +
  scale_y_log10() +
  xlab("Year (1990-2020)") + 
  ylab("Expected burn area (ha)") +
  theme_classic() + theme(legend.position = "none")
ggsave("full-model/figures/paper/burns_preds-vs-truth_gamma-ri_old.pdf", width = 15)

p <- burn_preds_gamma_ri_winsor %>% 
  ggplot(aes(x = year, y = winsor_total, group = year, color = train)) + 
  geom_boxplot(outlier.size = 0.2) + 
  scale_color_grey(start = 0.4, end = 0.6) +
  geom_point(inherit.aes = FALSE, data = true_burns_level1_full, aes(x = year, y = true_area), col = "red", size = 0.35) +
  geom_line(inherit.aes = FALSE, data = true_burns_level1_full, aes(x = year, y = true_area), col = "red", linewidth = 0.35) +
  facet_wrap(. ~ NA_L1NAME, scales = "free_y", nrow = 2) +
  scale_y_log10() +
  xlab("Year (1990-2020)") + 
  ylab("Expected burn area (ha)") +
  theme_classic() + theme(legend.position = "none")
ggsave("full-model/figures/paper/burns_preds-vs-truth_winsor_gamma-ri_old.pdf", width = 15)


# p <- burn_preds_both %>% 
#   ggplot(aes(x = year, y = total_area, group = interaction(year, model), color = model)) + 
#   geom_boxplot(outlier.size = 0.2) + 
#   geom_point(inherit.aes = FALSE, data = true_burns_level1_full, aes(x = year, y = true_area), col = "red", size = 0.35) +
#   geom_line(inherit.aes = FALSE, data = true_burns_level1_full, aes(x = year, y = true_area), col = "red", linewidth = 0.35) +
#   facet_wrap(. ~ NA_L1NAME, scales = "free_y", nrow = 4) +
#   scale_y_log10() +
#   xlab("Year (1990-2020)") + 
#   ylab("Expected burn area (ha)") +
#   theme_classic() + theme(legend.position = "none")
# ggsave("full-model/figures/paper/burns_preds-vs-truth_both.pdf", width = 15, height = 12)
# 
# p <- burn_preds_both_winsor %>% 
#   ggplot(aes(x = year, y = winsor_total, group = interaction(year, model), color = model)) + 
#   geom_boxplot(outlier.size = 0.2) + 
#   geom_point(inherit.aes = FALSE, data = true_burns_level1_full, aes(x = year, y = true_area), col = "black", size = 0.55) +
#   geom_line(inherit.aes = FALSE, data = true_burns_level1_full, aes(x = year, y = true_area), col = "black", linewidth = 0.55) +
#   facet_wrap(. ~ NA_L1NAME, scales = "free_y", nrow = 4) +
#   scale_y_log10() +
#   xlab("Year (1990-2020)") + 
#   ylab("Expected burn area (ha)") +
#   theme_classic()
# ggsave("full-model/figures/paper/burns_preds-vs-truth_winsor_both.pdf", width = 15, height = 12)

## compare theta from both joint models
theta_both <- theta_gamma_ri %>% 
  mutate(model = 'ri') %>% 
  rbind(theta_gamma_cst %>% mutate(model = 'cst')) %>%
  mutate(model = as.factor(model)) %>%
  group_by(timepoint, model) %>%
  summarize(med_theta = median(theta)) %>%
  ungroup()

p <- theta_both %>% rename(time = timepoint) %>% left_join(time_df) %>%
  ggplot(aes(x = date, y = med_theta, group = model, color = model)) +
  geom_line() +
  scale_x_date(name = "Year (1990-2020)", date_breaks = "5 years", date_labels = "%Y") +
  theme_classic()
ggsave("full-model/figures/paper/theta_gamma-ri_vs_gamma-cst.pdf", width = 15)

## maps of xi, sigma, pi_prob, gamma ---------
gamma_map <- gamma %>% group_by(region) %>% summarize(gamma = median(gamma)) %>% ungroup()
breaks <- classIntervals(c(min(gamma_map$gamma) - .00001, gamma_map$gamma), style = 'fixed', 
                         fixedBreaks = seq(-1.5, 1.5, length.out = 9), intervalClosure = 'left')
eco_gamma <- ecoregions_geom %>%
  left_join(gamma_map %>% left_join(full_reg_key)) %>%
  mutate(gamma_cat = cut(gamma, unique(breaks$brks)))
p <- ecoregions_geom %>%
  ggplot() +
  geom_sf(size = .1, fill = 'white') +
  geom_sf(data = eco_gamma,
          aes(fill=gamma), alpha = 0.6, lwd = 0, inherit.aes = FALSE) +
  theme_void() + scale_fill_gradient2()
ggsave("full-model/figures/paper/gamma_map.pdf")

# sigma_map <- rand_int %>% select(-xi) %>% group_by(region) %>% summarize(sigma = median(sigma))
# breaks <- classIntervals(c(min(sigma_map$sigma) - .00001, sigma_map$sigma), style = 'quantile', intervalClosure = 'left')
# eco_sigma <- ecoregions_geom %>% 
#   left_join(sigma_map %>% left_join(full_reg_key)) %>% 
#   mutate(sigma_cat = cut(sigma, unique(breaks$brks)))
# p <- ecoregions_geom %>%
#   ggplot() +
#   geom_sf(size = .1, fill = 'white') +
#   geom_sf(data = eco_sigma,
#           aes(fill=sigma_cat), alpha = 0.6, lwd = 0, inherit.aes = FALSE) +
#   theme_void() + scale_fill_brewer(palette = "YlOrRd")
# ggsave("full-model/figures/paper/sigma_map.pdf", dpi = 320, bg ='white')
# 
xi_map_cst <- rand_int_cst %>% select(-sigma) %>% group_by(region) %>% summarize(xi = median(xi))
breaks <- classIntervals(c(min(xi_map_cst$xi) - .00001, xi_map_cst$xi), style = 'fixed', 
                         fixedBreaks = c(0, 0.2, 0.4, 0.6, 0.8, 2.0), intervalClosure = 'left')
eco_xi <- ecoregions_geom %>% 
  left_join(xi_map_cst %>% left_join(full_reg_key)) %>% 
  mutate(xi_cat = cut(xi, unique(breaks$brks)))
p1 <- ecoregions_geom %>%
  ggplot() +
  geom_sf(size = .1, fill = 'white') +
  geom_sf(data = eco_xi,
          aes(fill=xi_cat), alpha = 0.6, lwd = 0, inherit.aes = FALSE) +
  theme_void() + 
  scale_fill_brewer(palette = "YlOrRd", 
                    labels = c(expression(""<="0.2"), "(0.2, 0.4]", "(0.4, 0.6]", "(0.6, 0.8]", expression("">"0.8")),
                    expression(xi))
ggsave("full-model/figures/paper/xi_map_gamma-cst.pdf", dpi = 320, bg ='white')

xi_map_eda <- burn_eda
breaks <- classIntervals(c(min(xi_map_eda$shape_include) - .00001, xi_map_eda$shape_include), style = 'fixed', 
                         fixedBreaks = c(0, 0.2, 0.4, 0.6, 0.8, 2.0), intervalClosure = 'left')
eco_eda <- ecoregions_geom %>% left_join(xi_map_eda) %>%
  mutate(xi_cat = cut(shape_include, unique(breaks$brks)))

p2 <- ecoregions_geom %>%
  ggplot() +
  geom_sf(size = .1, fill = 'white') +
  geom_sf(data = eco_eda,
          aes(fill=xi_cat), alpha = 0.6, lwd = 0, inherit.aes = FALSE) +
  theme_void() + 
  scale_fill_brewer(palette = "YlOrRd") + theme(legend.position = "none")
p <- p1 + p2 + plot_layout(guides = "collect")
ggsave("full-model/figures/paper/xi_map_gamma-cst_witheda.pdf", dpi = 320, bg ='white')
# pi_prob_map <- pi_prob %>% 
#   mutate(pi_expit = exp(pi)/(1 + exp(pi))) %>%
#   group_by(region) %>% 
#   summarize(pi_prob = median(pi_expit))
# breaks <- classIntervals(c(min(pi_prob_map$pi_prob) - .00001, pi_prob_map$pi_prob), style = 'quantile', intervalClosure = 'left')
# eco_pi_prob <- ecoregions_geom %>% 
#   left_join(pi_prob_map %>% left_join(full_reg_key)) %>% 
#   mutate(pi_prob_cat = cut(pi_prob, unique(breaks$brks)))
# p <- ecoregions_geom %>%
#   ggplot() +
#   geom_sf(size = .1, fill = 'white') +
#   geom_sf(data = eco_pi_prob,
#           aes(fill=pi_prob_cat), alpha = 0.6, lwd = 0, inherit.aes = FALSE) +
#   theme_void() + scale_fill_brewer(palette = "YlOrRd")
# ggsave("full-model/figures/paper/pi_prob_map.pdf", dpi = 320, bg ='white')

## partial effects plots for lambda and kappa ---------
beta_count <- readRDS("~/research/egpd-fires/full-model/figures/paper/mcmc_draws/theta-time_gamma-cst_2000iter/beta_count.RDS")
r <- 84
t <- 12*21 

stan_data <- readRDS("full-model/data/stan_data_joint.RDS")
X <- stan_data$X_train_count
vars <- c('log_housing_density', 'vs',
          'pr', 'prev_12mo_precip', 'tmmx',
          'rmin')
X_covar <- c()
X_cols <- vector("list", length(vars))
start <- 2
for(i in seq_along(vars)) {
  X_covar[i] <- paste0("X_", vars[i])
  X_cols[[i]] <- c(1, start:(start+5))
  start = start + 6
}

lambda_betas <- beta_count %>%
  separate_wider_delim(cols = "name", delim = ",", names = c("coef", "region")) %>%
  mutate(coef = as.numeric(gsub("beta_count\\[", "", coef)),
         region = as.numeric(gsub("\\]", "", region))) %>%
  group_by(region, coef) %>% summarize(med_val = median(value)) %>% ungroup() %>%
  pivot_wider(names_from = "region", values_from = "med_val") %>% select(-coef) %>% as.matrix()

# rescale data
un_std_data <- readRDS("full-model/data/un_std_all.RDS")
X_unstd <- X
for(k in seq_along(vars)) {
  var_mean <- un_std_data %>% filter(variable == vars[[k]]) %>% select(mean) %>% as.numeric()
  var_sd <- un_std_data %>% filter(variable == vars[[k]]) %>% select(sd) %>% as.numeric()
  X_unstd[,,X_cols[[k]][2]] <- X[,,X_cols[[k]][2]] * var_sd + var_mean
}

reg_cols <- full_reg_key$region
covar_effect <- function(egpd_param_df, covar_term, linear_term) {
  return(
    egpd_param_df %>% as_tibble() %>% rename_with(., ~ as.character(reg_cols)) %>%
      mutate(time = c(1:t)) %>%
      pivot_longer(cols = c(1:all_of(r)), values_to = "effect", names_to = "region") %>%
      mutate(region = as.numeric(region), covar = covar_term, linear = linear_term)
  )
}

coef_df_list_lambda <- list()
for(k in seq_along(vars)) {
  var_mean <- un_std_data %>% filter(variable == vars[[k]]) %>% select(mean) %>% as.numeric()
  var_sd <- un_std_data %>% filter(variable == vars[[k]]) %>% select(sd) %>% as.numeric()
  stored_df_lambda <- matrix(NA, t, r)
  for(j in 1:r) {
    stored_df_lambda[, j] <- X[j, , X_cols[[k]][3:7]] %*% lambda_betas[X_cols[[k]][3:7], j] +
      X_unstd[j, , X_cols[[k]][2]] * lambda_betas[X_cols[[k]][2], j]/var_sd +
      lambda_betas[1, j] - (lambda_betas[X_cols[[k]][2], j] * var_mean)/var_sd
  }
  coef_df_list_lambda[[k]] <- covar_effect(stored_df_lambda, vars[k], c(X_unstd[,,X_cols[[k]][2]]))
}

lambda_effects <- bind_rows(coef_df_list_lambda) %>% as_tibble() %>% left_join(., full_reg_key)
p <- lambda_effects %>% 
  mutate(linear = case_when(covar == 'tmmx' ~ linear - 273.15,
                            TRUE ~ linear),
         covar = case_when(covar == 'log_housing_density' ~ "log(housing density (units/sq. km))",
                           covar == 'pr' ~ 'Precipitation (mm), same month',
                           covar == 'prev_12mo_precip' ~ 'Precipitation (mm), past 12months',
                           covar == 'rmin' ~ 'Min. relative humidity (%)',
                           covar == 'tmmx' ~ 'Max. air temperature (C)',
                           covar == 'vs' ~ 'Wind speed (m/s)',
                           TRUE ~ covar)) %>% 
  ggplot(aes(x = linear, y = effect, group = region)) + 
  geom_line(aes(color = NA_L2CODE)) +
  facet_wrap(. ~ covar, scales = "free_x") + theme_classic() + theme(legend.position = "none") +
  ylab("Partial effect") + xlab("")
file_name <- "full-model/figures/paper/partial_effects_lambda.pdf"
ggsave(file_name, p, dpi = 320, width = 14, height = 8, bg = "white")

# partial effects for kappa
X <- stan_data$X_train_burn
vars <- c('log_housing_density', 'erc')
X_covar <- c()
X_cols <- vector("list", length(vars))
start <- 2
for(i in seq_along(vars)) {
  X_covar[i] <- paste0("X_", vars[i])
  X_cols[[i]] <- c(1, start:(start+5))
  start = start + 6
}

kappa_betas <- beta_burn %>%
  separate_wider_delim(cols = "name", delim = ",", names = c("coef", "region")) %>%
  mutate(coef = as.numeric(gsub("beta_burn\\[", "", coef)),
         region = as.numeric(gsub("\\]", "", region))) %>%
  group_by(region, coef) %>% summarize(med_val = median(value)) %>% ungroup() %>%
  pivot_wider(names_from = "region", values_from = "med_val") %>% select(-coef) %>% as.matrix()

# rescale data
X_unstd <- X
for(k in seq_along(vars)) {
  var_mean <- un_std_data %>% filter(variable == vars[[k]]) %>% select(mean) %>% as.numeric()
  var_sd <- un_std_data %>% filter(variable == vars[[k]]) %>% select(sd) %>% as.numeric()
  X_unstd[,,X_cols[[k]][2]] <- X[,,X_cols[[k]][2]] * var_sd + var_mean
}

coef_df_list_kappa <- list()
for(k in seq_along(vars)) {
  var_mean <- un_std_data %>% filter(variable == vars[[k]]) %>% select(mean) %>% as.numeric()
  var_sd <- un_std_data %>% filter(variable == vars[[k]]) %>% select(sd) %>% as.numeric()
  stored_df_kappa <- matrix(NA, t, r)
  for(j in 1:r) {
    stored_df_kappa[, j] <- X[j, , X_cols[[k]][3:7]] %*% kappa_betas[X_cols[[k]][3:7], j] +
      X_unstd[j, , X_cols[[k]][2]] * kappa_betas[X_cols[[k]][2], j]/var_sd +
      kappa_betas[1, j] - (kappa_betas[X_cols[[k]][2], j] * var_mean)/var_sd
  }
  coef_df_list_kappa[[k]] <- covar_effect(stored_df_kappa, vars[k], c(X_unstd[,,X_cols[[k]][2]]))
}

kappa_effects <- bind_rows(coef_df_list_kappa) %>% as_tibble() %>% left_join(., full_reg_key)
p <- kappa_effects %>% 
  mutate(covar = case_when(covar == 'log_housing_density' ~ "log(housing density (units/sq. km))",
                           covar == 'erc' ~ 'Energy Release Component',
                           TRUE ~ covar)) %>% 
  ggplot(aes(x = linear, y = effect, group = region)) + 
  geom_line(aes(color = NA_L2CODE)) +
  facet_wrap(. ~ covar, scales = "free_x") + theme_classic() + theme(legend.position = "none") +
  ylab("Partial effect") + xlab("")
file_name <- "full-model/figures/paper/partial_effects_kappa.pdf"
ggsave(file_name, p, dpi = 320, width = 14, height = 8, bg = "white")

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
# 
