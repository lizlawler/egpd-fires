library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(tidyverse)
library(stringr)
library(posterior)
library(lubridate)
library(sf)
library(classInt)
library(RColorBrewer)
library(patchwork)
library(colorspace)

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
levels_l1 <- levels(as.factor(as.numeric(ecoregions$NA_L1CODE)))
levels_l2 <- levels(as.factor(as.numeric(ecoregions$NA_L2CODE)))
ecoregions_geom <- ecoregions %>% filter(!NA_L2NAME == "UPPER GILA MOUNTAINS (?)") %>% 
  mutate(NA_L2CODE = factor(NA_L2CODE, levels = levels_l2),
         NA_L1CODE = factor(NA_L1CODE, levels = levels_l1),
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

## ecoregion maps ----------
l1_only_map <- ecoregions_geom %>%
  ggplot() +
  geom_sf(size = .1, fill = "transparent") +
  geom_sf(data = ecoregions_geom %>% group_by(NA_L1CODE) %>% summarise(geometry = st_union(st_set_precision(geometry, 1e8))),
          aes(fill = NA_L1CODE), alpha = 0.8, lwd = 0.5, inherit.aes = FALSE, show.legend = FALSE) +
  geom_sf(data = ecoregions_geom %>% group_by(NA_L3CODE) %>% summarise(geometry = st_union(st_set_precision(geometry, 1e8))),
          fill = "transparent", lwd = 0.1, inherit.aes = FALSE) +
  theme_void() +
  coord_sf(ndiscr = FALSE)
ggsave(filename = "full-model/figures/paper/level1_map.png", plot = l1_only_map,
       dpi = 320, width = 9, height = 9)
knitr::plot_crop("full-model/figures/paper/level1_map.png")

er_map_l1 <- ecoregions_geom %>%
  ggplot() +
  geom_sf(size = .1, fill = "transparent") +
  geom_sf(data = ecoregions_geom %>% group_by(NA_L2CODE) %>% summarise(geometry = st_union(st_set_precision(geometry, 1e8))),
          aes(fill = NA_L2CODE), alpha = 0.8, lwd = 0.25, inherit.aes = FALSE, show.legend = FALSE) +
  geom_sf(data = ecoregions_geom %>% group_by(NA_L1CODE) %>% summarise(geometry = st_union(st_set_precision(geometry, 1e8))),
          fill = "transparent", lwd = 0.8, color = "gray20", inherit.aes = FALSE, show.legend = FALSE) +
  geom_sf(data = ecoregions_geom %>% group_by(NA_L3CODE) %>% summarise(geometry = st_union(st_set_precision(geometry, 1e8))),
          fill = "transparent", lwd = 0.1, inherit.aes = FALSE) +
  theme_void() +
  coord_sf(ndiscr = FALSE)
ggsave(filename = "full-model/figures/paper/er_map_l1.png", plot = er_map_l1,
       dpi = 320, width = 9, height = 9)
knitr::plot_crop("full-model/figures/paper/er_map_l1.png")

er_map_l2 <- ecoregions_geom %>%
  ggplot() +
  geom_sf(size = .1, fill = "transparent") +
  geom_sf(data = ecoregions_geom %>% group_by(NA_L2CODE) %>% summarise(geometry = st_union(st_set_precision(geometry, 1e8))),
          aes(fill = NA_L2CODE), alpha = 0.8, lwd = 0.45, inherit.aes = FALSE, show.legend = FALSE) +
  geom_sf(data = ecoregions_geom %>% group_by(NA_L3CODE) %>% summarise(geometry = st_union(st_set_precision(geometry, 1e8))),
          fill = "transparent", lwd = 0.1, inherit.aes = FALSE) +
  theme_void() +
  coord_sf(ndiscr = FALSE)
ggsave(filename = "full-model/figures/paper/er_map_l2.png", plot = er_map_l2,
       dpi = 320, width = 9, height = 9)
knitr::plot_crop("full-model/figures/paper/er_map_l2.png")

er_map_l3 <- ecoregions_geom %>% 
  ggplot() +
  geom_sf(size = .1, fill = "white") +
  geom_sf(data = ecoregions_geom %>% group_by(NA_L3CODE) %>% summarise(geometry = st_union(st_set_precision(geometry, 1e8))),
          fill = "transparent", lwd = 0.3, inherit.aes = FALSE) +
  theme_void() +
  coord_sf(ndiscr = FALSE)
ggsave(filename = "full-model/figures/paper/er_map_l3.png", plot = er_map_l3,
       dpi = 320, width = 9, height = 9)
knitr::plot_crop("full-model/figures/paper/er_map_l3.png")

# read in extracted mcmc draws ------
# files <- paste0("full-model/figures/paper/mcmc_draws/theta-time_gamma-ri_erc_fwi/",
#                 list.files(path = "full-model/figures/paper/mcmc_draws/theta-time_gamma-ri_erc_fwi/",
#                            pattern = ".RDS"))
# files <- files[!grepl("old", files)]
# files <- files[!grepl("beta", files)]
# files <- files[!grepl("phi", files)]
# files <- files[!grepl("kappa", files)]


# create return level plots ------
kappa <- readRDS("./full-model/figures/paper/mcmc_draws/theta-time_gamma-ri_erc_fwi/kappa.RDS")
rand_int <- readRDS("./full-model/figures/paper/mcmc_draws/theta-time_gamma-ri_erc_fwi/rand_int.RDS")
lambda <- readRDS("./full-model/figures/paper/mcmc_draws/theta-time_gamma-ri_erc_fwi/lambda.RDS")
pi_prob <- readRDS("./full-model/figures/paper/mcmc_draws/theta-time_gamma-ri_erc_fwi/pi_prob.RDS")

all_params <- kappa %>% 
  left_join(rand_int) %>% 
  left_join(lambda) %>% 
  left_join(pi_prob) %>%
  mutate(eta = exp_count(pi, lambda))

## calculate expected burn area for the 98th quantile given there are eta fires in that ecoregion
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
# saveRDS(returns_regional, "full-model/figures/paper/tibbles/returns_regional_gamma-ri.RDS")

# returns_regional <- readRDS("~/Desktop/research/egpd-fires/full-model/figures/paper/tibbles/returns_regional_gamma-ri.RDS")
returns_regional <- returns_regional %>% mutate(NA_L1CODE = factor(NA_L1CODE, levels = levels_l1))
p <- returns_regional %>% ggplot() + 
  geom_ribbon(aes(x=date, ymin=lower, ymax=upper, group = region, fill = NA_L1CODE, alpha = 0.5)) +
  geom_line(aes(x=date, y=med50, group = region, alpha = 0.5), linewidth = 0.5, color = 'darkgrey') + 
  scale_y_log10() +
  scale_x_date(name = "Year (1990-2020)", date_breaks = "5 years", date_labels = "%Y") + 
  ylab("Expected burn area (ha)") +
  facet_wrap(. ~ NA_L1NAME, nrow = 2) +
  theme_classic() + 
  theme(legend.position = "none") + 
  theme(strip.text.x = element_text(size = rel(1.4)))
file_name <- "full-model/figures/paper/50yr_returns_erc_fwi.pdf"
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
# saveRDS(level98_areas_regional, "full-model/figures/paper/tibbles/level98_areas_regional_gamma-ri.RDS")
# level98_areas_regional <- readRDS("./full-model/figures/paper/tibbles/level98_areas_regional_gamma-ri.RDS")
level98_areas_regional <- level98_areas_regional %>% mutate(NA_L1CODE = factor(NA_L1CODE, levels = levels_l1))
p <- level98_areas_regional %>% ggplot() + 
  geom_ribbon(aes(x=date, ymin=lower, ymax=upper, group = region, fill = NA_L1CODE, alpha = 0.5)) +
  geom_line(aes(x=date, y=med50, group = region), linewidth = 0.5, color = 'darkgrey') + scale_y_log10() +
  scale_x_date(name = "Year (1990-2020)", date_breaks = "5 years", date_labels = "%Y") + 
  ylab("Expected burn area (ha) given fire occurrence") +
  facet_wrap(. ~ NA_L1NAME, nrow = 2) +
  theme_classic() + 
  theme(legend.position = "none") + 
  theme(strip.text.x = element_text(size = rel(1.4)))
file_name <- "full-model/figures/paper/98th_quant_burn-areas_erc_fwi.pdf"
ggsave(file_name, p, dpi = 320, bg = "white", width = 17, height = 9)

## area-weighted average of burn area exceedances -------
eco_areas <- ecoregions_geom %>% as_tibble() %>% group_by(NA_L3CODE) %>% summarise(area = sum(Shape_Area))
returns_level1 <- returns %>% select(c(draw, time, region, yr50)) %>% 
  left_join(full_reg_key) %>% 
  left_join(eco_areas) %>% 
  filter(yr50 != Inf) %>%
  group_by(NA_L1NAME, time, draw, NA_L1CODE) %>% 
  summarize(wmean = weighted.mean(yr50, area)) %>% 
  ungroup() %>% 
  group_by(time, NA_L1NAME, NA_L1CODE) %>%
  summarize(med50 = median(wmean), 
            lower = quantile(wmean, probs = 0.025), 
            upper = quantile(wmean, probs = 0.975)) %>%
  ungroup()

returns_level1 <- returns_level1 %>% mutate(NA_L1CODE = factor(NA_L1CODE, levels = levels_l1))
p <- returns_level1 %>%
  left_join(time_df) %>%
  ggplot() + 
  geom_ribbon(aes(x=date, ymin=lower, ymax=upper, group = NA_L1CODE, fill = NA_L1CODE, alpha = 0.95)) +
  geom_line(aes(x=date, y=med50, group = NA_L1CODE), linewidth = 0.5, color = 'darkgrey') + 
  scale_y_log10() +
  scale_x_date(name = "Year (1990-2020)", date_breaks = "5 years", date_labels = "%Y") + 
  ylab("Expected burn area (ha)") +
  facet_wrap(. ~ NA_L1NAME, nrow = 2) +
  theme_classic() + 
  theme(legend.position = "none") + 
  theme(strip.text.x = element_text(size = rel(1.1)))
file_name <- "full-model/figures/paper/50yr_returns_level1.pdf"
ggsave(file_name, p, dpi = 320, width = 15, height = 8)

level98_areas_level1 <- level98_areas %>% select(c(draw, time, region, yr50)) %>% 
  left_join(full_reg_key) %>% 
  left_join(eco_areas) %>% 
  filter(yr50 != Inf) %>%
  group_by(NA_L1NAME, time, draw, NA_L1CODE) %>% 
  summarize(wmean = weighted.mean(yr50, area)) %>% 
  ungroup() %>% 
  group_by(time, NA_L1NAME, NA_L1CODE) %>%
  summarize(med50 = median(wmean), 
            lower = quantile(wmean, probs = 0.025), 
            upper = quantile(wmean, probs = 0.975)) %>%
  ungroup()

level98_areas_level1 <- level98_areas_level1 %>% mutate(NA_L1CODE = factor(NA_L1CODE, levels = levels_l1))
p <- level98_areas_level1 %>%
  left_join(time_df) %>%
  ggplot() + 
  geom_ribbon(aes(x=date, ymin=lower, ymax=upper, group = NA_L1CODE, fill = NA_L1CODE, alpha = 0.95)) +
  geom_line(aes(x=date, y=med50, group = NA_L1CODE), linewidth = 0.5, color = 'darkgrey') + 
  scale_y_log10() +
  scale_x_date(name = "Year (1990-2020)", date_breaks = "5 years", date_labels = "%Y") + 
  ylab("Expected burn area (ha) given fire occurrence") +
  facet_wrap(. ~ NA_L1NAME, nrow = 2) +
  theme_classic() + 
  theme(legend.position = "none") + 
  theme(strip.text.x = element_text(size = rel(1.1)))
file_name <- "full-model/figures/paper/98th_quant_burn-areas_level1.pdf"
ggsave(file_name, p, dpi = 320, width = 15, height = 8)



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
stan_data <- readRDS("full-model/data/stan_data_joint_erc_fwi.RDS")
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
ggsave("full-model/figures/paper/counts_preds-vs-truth_erc_fwi.pdf", width = 15)

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
ggsave("full-model/figures/paper/counts_preds-vs-truth_winsor_erc_fwi.pdf", width = 15)


## create boxplot of predicted burn area for entire US annually for all timepoints (holdout and training) and overlay truth ------
burn_preds <- readRDS("full-model/figures/paper/mcmc_draws/theta-time_gamma-ri_erc_fwi/burn_pred.RDS")

# create for 'gamma-ri' model
burn_preds_l1 <- burn_preds %>% 
  left_join(full_reg_key) %>% 
  left_join(time_df) %>% 
  mutate(year = year(date)) %>% rename(preds = value) %>%
  group_by(NA_L1NAME, year, draw) %>% 
  summarize(total_area = sum(preds[is.finite(preds)])*1000*0.405) %>% 
  ungroup() %>% 
  mutate(train = case_when(year %in% train_years ~ TRUE,
                           year %in% test_years ~ FALSE))

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
# g1_burn_preds_winsor_limits <- g1_burns_full_annual %>% group_by(NA_L1NAME, year) %>%
#   summarize(lower = quantile(total_area, 0.05, na.rm = TRUE),
#             upper = quantile(total_area, 0.95, na.rm = TRUE)) %>%
#   ungroup()
# g1_burn_preds_winsor <- g1_burns_full_annual %>%
#   left_join(g1_burn_preds_winsor_limits) %>%
#   mutate(winsor_total = case_when(total_area <= lower ~ lower,
#                                   total_area >= upper ~ upper,
#                                   .default = total_area))

# G3 burn preds ------
# g3_burns_full <- rbind(g3_burn_preds_1_greaterthan500, g3_burn_preds_1_lessthan500) %>% 
#   rbind(rbind(g3_burn_preds_3_greaterthan500, g3_burn_preds_3_lessthan500) %>% mutate(draw = 1000 + draw))
# rm(list = ls(pattern = "than"))
# 
# g3_burns_full <- g3_burns_full %>% left_join(full_reg_key) %>% left_join(time_df) %>% mutate(year = year(date))
# g3_burns_full_annual <- g3_burns_full %>%
#   group_by(NA_L1NAME, year, draw) %>%
#   summarize(total_area = sum(preds[is.finite(preds)])*1000*0.405) %>%
#   ungroup() %>%
#   mutate(train = case_when(year %in% train_years ~ TRUE,
#                            year %in% test_years ~ FALSE))
# 
# g3_burn_preds_winsor_limits <- g3_burns_full_annual %>% group_by(NA_L1NAME, year) %>%
#   summarize(lower = quantile(total_area, 0.05, na.rm = TRUE),
#             upper = quantile(total_area, 0.95, na.rm = TRUE)) %>%
#   ungroup()
# g3_burn_preds_winsor <- g3_burns_full_annual %>%
#   left_join(g3_burn_preds_winsor_limits) %>%
#   mutate(winsor_total = case_when(total_area <= lower ~ lower,
#                                   total_area >= upper ~ upper,
#                                   .default = total_area))

# read in true burn areas -------
true_burns <- readRDS("full-model/data/burn_df_agg.RDS") %>% 
  filter(NA_L2NAME != "UPPER GILA MOUNTAINS (?)") %>%
  mutate(NA_L1NAME = as.factor(str_to_title(NA_L1NAME)))
true_burns_level1 <- true_burns %>% group_by(NA_L1NAME, fire_yr) %>%
  summarize(true_area = sum(total_burns)*0.405) %>% ungroup()

true_burns_level1_full <- true_burns_level1 %>% complete(NA_L1NAME, fire_yr) %>% rename(year = fire_yr)
true_burns_full <- true_burns %>% 
  group_by(NA_L3NAME, fire_yr, NA_L1NAME) %>% 
  summarize(burn_er = sum(total_burns)) %>% 
  ungroup() %>% 
  full_join(full_reg_key)%>% 
  complete(NA_L3NAME, fire_yr) %>% 
  rename(year = fire_yr) %>% 
  filter(!is.na(year)) %>% 
  select(-c(NA_L1NAME, NA_L3CODE, NA_L2CODE, NA_L1CODE, region, NA_L1NAME)) %>% 
  left_join(full_reg_key)

# burn_preds_l1_annual <- burn_preds_l1 %>%
#   group_by(NA_L1NAME, year, draw, train) %>% 
#   summarize(total_area = sum(preds[is.finite(preds)])) %>% 
#   ungroup()

burn_preds_winsor_limits <- burn_preds_l1 %>% group_by(NA_L1NAME, year) %>%
  summarize(lower = quantile(total_area, 0.025, na.rm = TRUE),
            upper = quantile(total_area, 0.975, na.rm = TRUE)) %>%
  ungroup()
burn_preds_winsor <- burn_preds_l1 %>%
  left_join(burn_preds_winsor_limits) %>%
  mutate(winsor_total = case_when(total_area <= lower ~ lower,
                                  total_area >= upper ~ upper,
                                  .default = total_area))

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

# p <- g3_burns_full_annual %>%
#   ggplot(aes(x = year, y = total_area, group = year, color = train)) +
#   geom_boxplot(outlier.size = 0.2) + scale_color_grey(start = 0.4, end = 0.6) +
#   geom_point(inherit.aes = FALSE, data = true_burns_level1_full, aes(x = year, y = true_area), col = "red", size = 0.35) +
#   geom_line(inherit.aes = FALSE, data = true_burns_level1_full, aes(x = year, y = true_area), col = "red", linewidth = 0.35) +
#   facet_wrap(. ~ NA_L1NAME, scales = "free_y", nrow = 2) +
#   scale_y_log10() +
#   xlab("Year (1990-2020)") +
#   ylab("Expected burn area (ha)") +
#   theme_classic() + theme(legend.position = "none")
# ggsave("full-model/figures/paper/burns_preds-vs-truth_g3.pdf", width = 15)
# 
# p <- g3_burn_preds_winsor %>%
#   ggplot(aes(x = year, y = winsor_total, group = year, color = train)) +
#   geom_boxplot(outlier.size = 0.2) + scale_color_grey(start = 0.4, end = 0.6) +
#   geom_point(inherit.aes = FALSE, data = true_burns_level1_full, aes(x = year, y = true_area), col = "red", size = 0.35) +
#   geom_line(inherit.aes = FALSE, data = true_burns_level1_full, aes(x = year, y = true_area), col = "red", linewidth = 0.35) +
#   facet_wrap(. ~ NA_L1NAME, scales = "free_y", nrow = 2) +
#   scale_y_log10() +
#   xlab("Year (1990-2020)") +
#   ylab("Expected burn area (ha)") +
#   theme_classic() + theme(legend.position = "none")
# ggsave("full-model/figures/paper/burns_preds-vs-truth_winsor_g3.pdf", width = 15)

p <- burn_preds_l1 %>% 
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
ggsave("full-model/figures/paper/burns_preds-vs-truth_erc_fwi.pdf", width = 15)

p <- burn_preds_winsor %>% 
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
ggsave("full-model/figures/paper/burns_preds-vs-truth_winsor_erc_fwi.pdf", width = 15)

## maps of xi, sigma, pi_prob, gamma ---------
gamma <- readRDS("full-model/figures/paper/mcmc_draws/theta-time_gamma-ri_erc_fwi/gamma.RDS")
gamma_map <- gamma %>% group_by(region) %>% summarize(gamma = median(gamma)) %>% ungroup()
# breaks <- classIntervals(c(min(gamma_map$gamma) - .00001, gamma_map$gamma), style = 'fixed', 
#                          fixedBreaks = seq(-1.5, 1.5, length.out = 9), intervalClosure = 'left')
eco_gamma <- ecoregions_geom %>%
  left_join(gamma_map %>% left_join(full_reg_key)) 
# %>%
#   mutate(gamma_cat = cut(gamma, unique(breaks$brks)))
p <- ecoregions_geom %>%
  ggplot() +
  geom_sf(size = .1, fill = 'transparent') +
  geom_sf(data = eco_gamma %>% group_by(NA_L3CODE, gamma) %>% summarise(geometry = st_union(st_set_precision(geometry, 1e8))),
          aes(fill=gamma), alpha = 1, lwd = 0.1, inherit.aes = FALSE) +
  theme_void() + scale_fill_gradient2(name = bquote(gamma~value))
p_legend <- p + theme(legend.position = c(0.95, 0.4), 
                      legend.key.size = unit(1.3, "cm"),
                      legend.text = element_text(size = 12),
                      legend.title = element_text(size = 15))
ggsave(filename = "full-model/figures/paper/gamma_map.png", plot = p_legend, 
       dpi = 320, 
       width = 10, height = 10)
knitr::plot_crop("full-model/figures/paper/gamma_map.png")
# parse(text = paste("tau[", 1:2, "]")))

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
burn_eda <- readRDS("./full-model/figures/paper/burn_eda.RDS")
eco_eda <- ecoregions_geom %>% left_join(burn_eda)

eda_90_plot <- ecoregions_geom %>%
  ggplot() +
  geom_sf(size = .1, fill = 'transparent') +
  geom_sf(data = eco_eda %>% group_by(NA_L2NAME, shape_inc_90) %>% summarise(geometry = st_union(st_set_precision(geometry, 1e8))), 
          aes(fill=shape_inc_90), alpha = 0.6, lwd = 0.1, inherit.aes = FALSE) +
  theme_void()
eda_90_hcl <- eda_90_plot + scale_fill_continuous_sequential(palette = 'Mint',
                                                             na.value = "transparent",
                                                             name = bquote(xi~value)) + 
  theme(legend.position = c(0.95, 0.4), 
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15))
ggsave(filename = "full-model/figures/paper/xi_eda_90th-quant.png", 
       plot = eda_90_hcl, 
       width = 10, height = 10,
       dpi = 320)
knitr::plot_crop("full-model/figures/paper/xi_eda_90th-quant.png")

# eda_90_gradient <- eda_90_plot + scale_fill_gradient2(na.value = "transparent",
#                                     name = bquote(xi~value)) + 
#   theme(legend.position = c(0.95, 0.4), 
#         legend.key.size = unit(1, "cm"),
#         legend.text = element_text(size = 12),
#         legend.title = element_text(size = 15))
# ggsave("full-model/figures/paper/xi_eda_90th-quant.pdf", width = 10, height = 8)

rand_int <- readRDS("./full-model/figures/paper/mcmc_draws/theta-time_gamma-ri_erc_fwi/rand_int.RDS")
xi_map <- rand_int %>% select(-sigma) %>% group_by(region) %>% summarize(xi = median(xi))
eco_xi <- ecoregions_geom %>% 
  left_join(xi_map %>% left_join(full_reg_key))

xi_model_map <- ecoregions_geom %>%
  ggplot() +
  geom_sf(size = .1, fill = 'transparent') +
  geom_sf(data = eco_xi %>% group_by(NA_L3CODE, xi) %>% summarise(geometry = st_union(st_set_precision(geometry, 1e8))),
          aes(fill=xi), alpha = 0.6, lwd = 0.1, inherit.aes = FALSE) +
  theme_void()
xi_model_hcl <- xi_model_map + scale_fill_continuous_sequential(palette = 'Mint',
                                                                na.value = "transparent",
                                                                name = bquote(xi~value),
                                                                limits = c(0,1.2),
                                                                breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)) + 
  theme(legend.position = c(0.9, 0.4), 
        legend.key.size = unit(1.2, "cm"),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 13))
# ggsave("full-model/figures/paper/xi_model_map.pdf", width = 10)

eda_90_combo_plot <- eda_90_plot + scale_fill_continuous_sequential(palette = 'Mint',
                                                                    na.value = "transparent",
                                                                    limits = c(0,1.2)) +
  theme(legend.position = "none")
# eda_95_combo_plot <- eda_95_plot + scale_fill_continuous_sequential(palette = 'Mint',
#                                                                     na.value = "transparent",
#                                                                     limits = c(0,1.2)) +
#   theme(legend.position = "none")

p90 <- eda_90_combo_plot + xi_model_hcl + plot_layout(guides = 'collect')
ggsave("full-model/figures/paper/xi_map_with-eda90.png", dpi = 320, width = 20, height = 20, bg = 'transparent')
knitr::plot_crop("full-model/figures/paper/xi_map_with-eda90.png")

sigma_map <- rand_int %>% select(-xi) %>% group_by(region) %>% summarize(sigma = median(sigma))
eco_sigma <- ecoregions_geom %>% 
  left_join(sigma_map %>% left_join(full_reg_key))

sigma_model_map <- ecoregions_geom %>%
  ggplot() +
  geom_sf(size = .1, fill = 'transparent') +
  geom_sf(data = eco_sigma %>% group_by(NA_L3CODE, sigma) %>% summarise(geometry = st_union(st_set_precision(geometry, 1e8))),
          aes(fill=sigma), alpha = 0.6, lwd = 0.1, inherit.aes = FALSE) +
  theme_void()

sigma_model_hcl <- sigma_model_map + scale_fill_continuous_sequential(palette = 'Mint',
                                                                na.value = "transparent",
                                                                name = bquote(sigma~value),
                                                                limits = c(0.25,4.5)) + 
  theme(legend.position = c(0.9, 0.4), 
        legend.key.size = unit(1.2, "cm"),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 13))
ggsave("full-model/figures/paper/sigma_map.png", dpi = 320, width = 20, height = 20, bg = 'transparent')
knitr::plot_crop("full-model/figures/paper/sigma_map.png")


# p95 <- eda_95_combo_plot + xi_model_hcl + plot_layout(guides = 'collect')
# ggsave("full-model/figures/paper/xi_map_with-eda95.pdf")


## partial effects plots for lambda and kappa ---------
beta_count <- readRDS("./full-model/figures/paper/mcmc_draws/theta-time_gamma-ri_erc_fwi/beta_count.RDS")
r <- 84
t <- 12*21 

stan_data <- readRDS("full-model/data/stan_data_joint_erc_fwi.RDS")
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

lambda_effects <- bind_rows(coef_df_list_lambda) %>% as_tibble() %>% left_join(., full_reg_key) %>% mutate(NA_L2CODE = factor(NA_L2CODE, levels = levels_l2))
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
file_name <- "full-model/figures/paper/partial_effects_lambda_erc_fwi.pdf"
ggsave(file_name, p, width = 15, height = 9)

# partial effects for kappa
beta_burn <- readRDS("./full-model/figures/paper/mcmc_draws/theta-time_gamma-ri_erc_fwi/beta_burn.RDS")
X <- stan_data$X_train_burn
vars <- c('log_housing_density', 'erc', 'fwi')
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

kappa_effects <- bind_rows(coef_df_list_kappa) %>% as_tibble() %>% left_join(., full_reg_key) %>% mutate(NA_L2CODE = factor(NA_L2CODE, levels = levels_l2))
p <- kappa_effects %>% 
  mutate(covar = case_when(covar == 'log_housing_density' ~ "log(housing density (units/sq. km))",
                           covar == 'erc' ~ 'Energy Release Component',
                           covar == 'fwi' ~ 'Fire Weather Index',
                           TRUE ~ covar)) %>% 
  ggplot(aes(x = linear, y = effect, group = region)) + 
  geom_line(aes(color = NA_L2CODE)) +
  facet_wrap(. ~ covar, scales = "free_x") + theme_classic() + theme(legend.position = "none") +
  ylab("Partial effect") + xlab("")
file_name <- "full-model/figures/paper/partial_effects_kappa_erc_fwi.pdf"
ggsave(file_name, p, width = 15, height = 4.5)

## look at rho values ----------
# 1=lambda, 2=kappa, 3=pi, 4=delta, 5=sigma, 6=xi, 7 = gamma

rho1 <- readRDS("~/research/egpd-fires/full-model/figures/paper/mcmc_draws/theta-time_gamma-ri_erc_fwi/rho1.RDS")
rho2 <- readRDS("~/research/egpd-fires/full-model/figures/paper/mcmc_draws/theta-time_gamma-ri_erc_fwi/rho2.RDS")

rho1 <- rho1 %>% mutate(rho1 = as.numeric(str_remove(str_remove(name, "rho1\\["), "\\]"))) %>%
  mutate(rho1 = case_when(rho1 == 1 ~ 'lambda',
                          rho1 == 2 ~ 'kappa',
                          rho1 == 3 ~ 'pi',
                          rho1 == 4 ~ 'delta',
                          rho1 == 5 ~ 'sigma',
                          rho1 == 6 ~ 'xi',
                          rho1 == 7 ~ 'gamma')) %>% 
  select(-name)

rho2 <- rho2 %>% mutate(rho2 = as.numeric(str_remove(str_remove(name, "rho2\\["), "\\]"))) %>%
  mutate(rho2 = case_when(rho2 == 1 ~ 'lambda',
                          rho2 == 2 ~ 'kappa',
                          rho2 == 3 ~ 'pi',
                          rho2 == 4 ~ 'delta',
                          rho2 == 5 ~ 'sigma',
                          rho2 == 6 ~ 'xi',
                          rho2 == 7 ~ 'gamma')) %>% 
  select(-name)

rho1_meds <- rho1 %>% group_by(rho1) %>% summarize(med_val = median(value)) %>% ungroup()
rho2_meds <- rho2 %>% group_by(rho2) %>% summarize(med_val = median(value)) %>% ungroup()
rho_full <- rho1_meds %>% 
  rename(param = rho1, rho1 = med_val) %>% 
  left_join(rho2_meds %>% rename(param = rho2, rho2 = med_val)) %>%
  mutate(rho_sum = rho1 + rho2,
         model_part = case_when(param == 'delta' ~ 'ZINB',
                                param == 'gamma' ~ 'joint',
                                param == 'kappa' ~ 'EGPD',
                                param == 'lambda' ~ 'ZINB', 
                                param == 'pi' ~ 'ZINB',
                                param == 'sigma' ~ 'EGPD', 
                                param == 'xi' ~ 'EGPD')) %>%
  relocate(model_part, .before = param) %>% arrange(model_part)


