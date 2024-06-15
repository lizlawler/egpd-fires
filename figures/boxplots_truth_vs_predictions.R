##############################################################################
##  This script produces two sets of boxplots comparing:                    ##
##  (1) total expected wildfire counts to total observed counts             ##
##  (2) total expected burned area to total observed burned area            ##
##############################################################################

library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(tidyverse)
library(stringr)
library(posterior)
library(lubridate)
source("./data/load_eco.R")

# load ecoregion shapes -----------------------------
ecoregion_shp <- load_ecoregions()
ecoregion_df <- ecoregion_shp %>% sf::st_as_sf()
levels_l1 <- levels(as.factor(as.numeric(ecoregion_df$NA_L1CODE)))
levels_l2 <- levels(as.factor(as.numeric(ecoregion_df$NA_L2CODE)))
## read in ecoregion key ---------------------------------------------------
region_key <- readRDS(file = "./data/processed/region_key.rds") %>%
  mutate(NA_L2CODE = factor(NA_L2CODE, levels = levels_l2),
         NA_L1CODE = factor(NA_L1CODE, levels = levels_l1),
         NA_L3CODE = as.factor(NA_L3CODE),
         NA_L1NAME = as.factor(str_to_title(NA_L1NAME)))

## create dataframe of timepoints -------------------------------------------
all_years <- 1990:2020
first_five <- 1990:1994
last_five <- 2020:2016
test_years <- sort(c(first_five, last_five))
train_years <- setdiff(all_years, test_years)

date_seq <- seq(as.Date("1990-01-01"), by = "1 month", length.out = 372) %>% as_tibble() %>% rename(date = value)
time_df <- date_seq %>% mutate(time = 1:372)

## expected counts from ZINB posterior parameters -----------------------------
exp_count <- function(pi, lambda) {
  pi_prob <- exp(pi)/(1+exp(pi))
  return((1-pi_prob) * exp(lambda))
}

## read in extracted MCMC draws -----------------------------------------------
lambda <- readRDS("./figures/mcmc_draws/lambda.RDS")
pi_prob <- readRDS("./figures/mcmc_draws/pi_prob.RDS")
size_preds <- readRDS("./figures/mcmc_draws/size_pred.RDS")

## calculations for plot (1): wildfire count predictions vs truth -------------
count_preds <- lambda %>% 
  left_join(pi_prob) %>% 
  mutate(preds = exp_count(pi, lambda)) %>% 
  select(c(draw, time, region, preds)) %>%
  left_join(region_key) %>% 
  left_join(time_df) %>% 
  mutate(year = year(date)) %>% 
  group_by(NA_L1NAME, year, draw) %>% 
  summarize(total_fires = sum(preds)) %>% 
  ungroup() %>% 
  mutate(train = case_when(year %in% train_years ~ TRUE,
                           year %in% test_years ~ FALSE))

# calculate 95% winsorization for count predictions 
preds_winsor_limits <- count_preds %>% 
  group_by(NA_L1NAME, year) %>% 
  summarize(lower = quantile(total_fires, 0.025),
            upper = quantile(total_fires, 0.975)) %>%
  ungroup()
count_preds_winsor <- count_preds %>% 
  left_join(preds_winsor_limits) %>%
  mutate(winsor_total = case_when(total_fires <= lower ~ lower,
                                  total_fires >= upper ~ upper,
                                  .default = total_fires))

## read in observed wildfire counts ------------------------------------------
stan_data <- readRDS("./data/stan_lists/data_joint.RDS")
y_train_count <- stan_data$y_train_count %>% 
  as_tibble() %>% 
  rowid_to_column(var = "time") %>% 
  pivot_longer(!time, names_to = "region", values_to = "value") %>%
  mutate(region = as.numeric(gsub("V", "", region)),
         time = time + 60) %>%
  left_join(region_key) %>%
  left_join(time_df) %>% 
  mutate(year = year(date)) %>%
  group_by(NA_L1NAME, year) %>%
  summarize(true_count = sum(value)) %>%
  ungroup()
y_hold_count <- stan_data$y_hold_count %>% 
  as_tibble() %>% 
  rowid_to_column(var = "time") %>% 
  pivot_longer(!time, names_to = "region", values_to = "value") %>%
  mutate(region = as.numeric(gsub("V", "", region)),
         time = case_when(time > 60 ~ time + 252,
                          TRUE ~ time)) %>%
  left_join(region_key) %>%
  left_join(time_df) %>% 
  mutate(year = year(date)) %>%
  group_by(NA_L1NAME, year) %>%
  summarize(true_count = sum(value)) %>%
  ungroup()
true_counts <- bind_rows(y_train_count, y_hold_count) %>% 
  arrange(year, NA_L1NAME, locale = ".en")

# counts_boxplot <- count_preds %>% 
#   ggplot(aes(x = year, y = total_fires, group = year, color = train)) + 
#   geom_boxplot(outlier.size = 0.2) + 
#   scale_color_grey(start = 0.4, end = 0.6) +
#   geom_point(inherit.aes = FALSE, data = true_counts, aes(x = year, y = true_count), col = "red", size = 0.35) +
#   geom_line(inherit.aes = FALSE, data = true_counts, aes(x = year, y = true_count), col = "red", linewidth = 0.35) +
#   facet_wrap(. ~ NA_L1NAME, scales = "free_y", nrow = 2) + 
#   xlab("Year (1990-2020)") + 
#   ylab("Expected number of fires") +
#   theme_classic() + 
#   theme(legend.position = "none")
# ggsave("full-model/figures/paper/counts_preds_vs_truth.pdf", counts_boxplot, 
#        dpi = 320, width = 8.5, height = 4)

counts_boxplot_winsor <- count_preds_winsor %>% 
  ggplot(aes(x = year, y = winsor_total, group = year, color = train)) + 
  geom_boxplot(outlier.size = 0.2) + 
  scale_color_grey(start = 0.4, end = 0.6) +
  geom_point(inherit.aes = FALSE, data = true_counts, aes(x = year, y = true_count), col = "red", size = 0.35) +
  geom_line(inherit.aes = FALSE, data = true_counts, aes(x = year, y = true_count), col = "red", linewidth = 0.35) +
  facet_wrap(. ~ NA_L1NAME, scales = "free_y", ncol = 2) + 
  xlab("Year (1990-2020)") + 
  ylab("Expected number of fires") +
  theme_classic() + 
  theme(legend.position = "none") 
ggsave("./figures/paper_figures/counts_preds_vs_truth.pdf", 
       counts_boxplot_winsor,
       bg = 'transparent',
       dpi = 320)
knitr::plot_crop("./figures/paper_figures/counts_preds_vs_truth.pdf")

## calculations for plot (2): wildfire burned area predictions vs truth -------
size_preds <- size_preds %>% 
  left_join(region_key) %>% 
  left_join(time_df) %>% 
  mutate(year = year(date)) %>% rename(preds = value) %>%
  group_by(NA_L1NAME, year, draw) %>% 
  summarize(total_area = sum(preds[is.finite(preds)])*1000*0.405) %>% 
  ungroup() %>% 
  mutate(train = case_when(year %in% train_years ~ TRUE,
                           year %in% test_years ~ FALSE))

size_preds_winsor_limits <- size_preds %>% 
  group_by(NA_L1NAME, year) %>%
  summarize(lower = quantile(total_area, 0.025, na.rm = TRUE),
            upper = quantile(total_area, 0.975, na.rm = TRUE)) %>%
  ungroup()
size_preds_winsor <- size_preds %>%
  left_join(size_preds_winsor_limits) %>%
  mutate(winsor_total = case_when(total_area <= lower ~ lower,
                                  total_area >= upper ~ upper,
                                  .default = total_area))

## read in observed wildfire burned area --------------------------------------
true_burns <- readRDS("./data/processed/obs_burned_areas.RDS") %>%
  mutate(NA_L1NAME = as.factor(str_to_title(NA_L1NAME)))
true_burns <- true_burns %>% group_by(NA_L1NAME, fire_yr) %>%
  summarize(true_area = sum(total_burns)*0.405) %>% 
  ungroup() %>% 
  complete(NA_L1NAME, fire_yr) %>% 
  rename(year = fire_yr)

# areas_boxplot <- size_preds %>% 
#   ggplot(aes(x = year, y = total_area, group = year, color = train)) + 
#   geom_boxplot(outlier.size = 0.2) + 
#   scale_color_grey(start = 0.4, end = 0.6) +
#   geom_point(inherit.aes = FALSE, data = true_burns_level1_full, aes(x = year, y = true_area), col = "red", size = 0.35) +
#   geom_line(inherit.aes = FALSE, data = true_burns_level1_full, aes(x = year, y = true_area), col = "red", linewidth = 0.35) +
#   facet_wrap(. ~ NA_L1NAME, scales = "free_y", nrow = 2) +
#   scale_y_log10() +
#   xlab("Year (1990-2020)") + 
#   ylab("Expected burn area (ha)") +
#   theme_classic() + theme(legend.position = "none",
# strip.text.x = element_text(size = rel(1.3)), 
# axis.text.x = element_text(size = rel(1.5)),
# axis.title.x = element_text(size = rel(1.4)),
# axis.text.y = element_text(size = rel(1.3)),
# axis.title.y = element_text(size = rel(1.4)))
# ggsave("full-model/figures/paper/areas_preds_vs_truth.pdf", width = 15)

areas_boxplot_winsor <- size_preds_winsor %>% 
  ggplot(aes(x = year, y = winsor_total, group = year, color = train)) + 
  geom_boxplot(outlier.size = 0.2) + 
  scale_color_grey(start = 0.4, end = 0.6) +
  geom_point(inherit.aes = FALSE, data = true_burns, aes(x = year, y = true_area), col = "red", size = 0.35) +
  geom_line(inherit.aes = FALSE, data = true_burns, aes(x = year, y = true_area), col = "red", linewidth = 0.35) +
  facet_wrap(. ~ NA_L1NAME, scales = "free_y", ncol = 2) +
  scale_y_log10() +
  xlab("Year (1990-2020)") + 
  ylab("Expected burn area (ha)") +
  theme_classic() + 
  theme(legend.position = "none")
ggsave("./figures/paper_figures/areas_preds_vs_truth.pdf", 
       areas_boxplot_winsor,
       bg = 'transparent',
       dpi = 320)
knitr::plot_crop("./figures/paper_figures/areas_preds_vs_truth.pdf")
