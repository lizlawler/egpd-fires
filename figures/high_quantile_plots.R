##############################################################################
##  This script produces two plots:                                         ##
##  (1) 50-year return levels of wildfire burned areas                      ##
##  (2) 98th quantile of the wildfire burned area distributions             ## 
##     marginalized over the distribution of wildfire counts                ##
##############################################################################

library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(tidyverse)
library(stringr)
library(lubridate)
library(sf)
library(RColorBrewer)
library(patchwork)
library(colorspace)
source("./data/load_eco.R")

## load in ecoregion shapes ----------------
ecoregion_shp <- load_ecoregions()
ecoregion_df <- ecoregion_shp %>% st_as_sf()
levels_l1 <- levels(as.factor(as.numeric(ecoregion_df$NA_L1CODE)))
levels_l2 <- levels(as.factor(as.numeric(ecoregion_df$NA_L2CODE)))
ecoregion_df <- ecoregion_df %>% 
  mutate(NA_L2CODE = factor(NA_L2CODE, levels = levels_l2),
         NA_L1CODE = factor(NA_L1CODE, levels = levels_l1),
         NA_L3CODE = as.factor(NA_L3CODE),
         NA_L1NAME = as.factor(str_to_title(NA_L1NAME)))

## read in L3 ecoregion key
region_key <- readRDS(file = "./data/processed/region_key.rds") %>%
  mutate(NA_L2CODE = factor(NA_L2CODE, levels = levels_l2),
         NA_L1CODE = factor(NA_L1CODE, levels = levels_l1),
         NA_L3CODE = as.factor(NA_L3CODE),
         NA_L1NAME = as.factor(str_to_title(NA_L1NAME)))

## create dataframe of timepoints -------------------------------------------
date_seq <- seq(as.Date("1990-01-01"), by = "1 month", length.out = 372) %>% as_tibble() %>% rename(date = value)
time_df <- date_seq %>% mutate(time = 1:372)

## G1 EGPD distributional functions -------------------------------------------
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

## expected counts from ZINB posterior parameters -----------------------------
exp_count <- function(pi, lambda) {
  pi_prob <- exp(pi)/(1+exp(pi))
  return((1-pi_prob) * exp(lambda))
}

## exceedance probabilities function ------------------------------------------
rlevel <- function(N, kappa, sigma, xi, eta) {
  p <- ((N-1)/N)^(1/(12 * eta)) * (1-pegpd(1.001, kappa, sigma, xi)) + pegpd(1.001, kappa, sigma, xi)
  return(qegpd(p, kappa, sigma, xi))
}

## high quantiles function ----------------------------------------------------
high_quant <- function(N, kappa, sigma, xi) {
  p <- ((N-1)/N) * (1-pegpd(1.001, kappa, sigma, xi)) + pegpd(1.001, kappa, sigma, xi)
  return(qegpd(p, kappa, sigma, xi))
}

## read in extracted MCMC draws -----------------------------------------------
kappa <- readRDS("./figures/mcmc_draws/kappa.RDS")
rand_int <- readRDS("./figures/mcmc_draws/rand_int.RDS")
lambda <- readRDS("./figures/mcmc_draws/lambda.RDS")
pi_prob <- readRDS("./figures/mcmc_draws/pi_prob.RDS")

all_params <- kappa %>% 
  left_join(rand_int) %>% 
  left_join(lambda) %>% 
  left_join(pi_prob) %>%
  mutate(eta = exp_count(pi, lambda))

## calculation for plot (1): 50-year return levels for each MCMC iteration --------------------
# rescale back to 1000s of acres; convert to hectares
returns <- all_params %>%
  mutate(yr50 = rlevel(50, kappa, sigma, xi, eta)*1000*0.405) 

# calculate median and 95% credible interval for every month and L3 ecoregion 
returns_summary <- returns %>% group_by(time, region) %>%
  summarize(med50 = median(yr50[is.finite(yr50)]),
            lower = quantile(yr50[is.finite(yr50)], probs = 0.025),
            upper = quantile(yr50[is.finite(yr50)], probs = 0.975)) %>%
  ungroup()
returns_regional <- returns_summary %>%
  left_join(region_key) %>%
  left_join(time_df) %>% 
  mutate(NA_L1CODE = factor(NA_L1CODE, levels = levels_l1))

returns_plot <- returns_regional %>% 
  ggplot() + 
  geom_ribbon(aes(x=date, ymin=lower, ymax=upper, group = region, fill = NA_L1CODE, alpha = 0.5)) +
  geom_line(aes(x=date, y=med50, group = region, alpha = 0.5), linewidth = 0.5, color = 'darkgrey') + 
  scale_y_log10() +
  scale_x_date(name = "Year (1990-2020)", date_breaks = "5 years", date_labels = "%Y") + 
  ylab("Expected burn area (ha)") +
  facet_wrap(. ~ NA_L1NAME, ncol = 2) +
  theme_classic() + 
  theme(legend.position = "none",
        strip.text.x = element_text(size = rel(1.3)),
        axis.text.x = element_text(size = rel(1.3)),
        axis.title.x = element_text(size = rel(1.3)),
        axis.text.y = element_text(size = rel(1.3)),
        axis.title.y = element_text(size = rel(1.3)))

ggsave("./figures/paper_figures/50yr_returns.pdf", 
       returns_plot, 
       dpi = 320, 
       bg = "white", 
       width = 8.5, 
       height = 8.5)
knitr::plot_crop("./figures/paper_figures/50yr_returns.pdf")

## calculation for plot (2): 98th quantile of expected burn areas -------------
quant98 <- all_params %>%
  mutate(yr50 = high_quant(50, kappa, sigma, xi)*1000*0.405) 
quant98_summary <- quant98 %>% group_by(time, region) %>%
  summarize(med50 = median(yr50[is.finite(yr50)]),
            lower = quantile(yr50[is.finite(yr50)], probs = 0.025),
            upper = quantile(yr50[is.finite(yr50)], probs = 0.975)) %>%
  ungroup()
quant98_summary_regional <- quant98_summary %>%
  left_join(region_key) %>%
  left_join(time_df) %>% 
  mutate(NA_L1CODE = factor(NA_L1CODE, levels = levels_l1))

quant98_plot <- quant98_summary_regional %>% 
  ggplot() + 
  geom_ribbon(aes(x=date, ymin=lower, ymax=upper, group = region, fill = NA_L1CODE, alpha = 0.5)) +
  geom_line(aes(x=date, y=med50, group = region), linewidth = 0.5, color = 'darkgrey') + scale_y_log10() +
  scale_x_date(name = "Year (1990-2020)", date_breaks = "5 years", date_labels = "%Y",  ) + 
  ylab("Expected burn area (ha) given fire occurrence") +
  facet_wrap(. ~ NA_L1NAME, ncol = 2) +
  theme_classic() + 
  theme(legend.position = "none",
        strip.text.x = element_text(size = rel(1.3)),
        axis.text.x = element_text(size = rel(1.3)),
        axis.title.x = element_text(size = rel(1.3)),
        axis.text.y = element_text(size = rel(1.3)),
        axis.title.y = element_text(size = rel(1.3)))
ggsave("./figures/paper_figures/98th_quant_sizes.pdf", 
       quant98_plot, 
       dpi = 320, 
       bg = "white", 
       width = 8.5, 
       height = 8.5)
knitr::plot_crop("./figures/paper_figures/98th_quant_sizes.pdf")

## calculate area-weighted average of the above values by L1 ecoregion -------
## these plots were not included in the paper
# eco_areas <- ecoregion_df %>% 
#   as_tibble() %>% 
#   group_by(NA_L3CODE) %>%
#   summarise(area = sum(Shape_Area)) 
# returns_level1 <- returns %>% select(c(draw, time, region, yr50)) %>% 
#   left_join(region_key) %>% 
#   left_join(eco_areas) %>% 
#   filter(yr50 != Inf) %>%
#   group_by(NA_L1NAME, time, draw, NA_L1CODE) %>% 
#   summarize(wmean = weighted.mean(yr50, area)) %>% 
#   ungroup() %>% 
#   group_by(time, NA_L1NAME, NA_L1CODE) %>%
#   summarize(med50 = median(wmean), 
#             lower = quantile(wmean, probs = 0.025), 
#             upper = quantile(wmean, probs = 0.975)) %>%
#   ungroup() %>% mutate(NA_L1CODE = factor(NA_L1CODE, levels = levels_l1))
# returns_level1_plot <- returns_level1 %>%
#   left_join(time_df) %>%
#   ggplot() + 
#   geom_ribbon(aes(x=date, ymin=lower, ymax=upper, group = NA_L1CODE, fill = NA_L1CODE, alpha = 0.95)) +
#   geom_line(aes(x=date, y=med50, group = NA_L1CODE), linewidth = 0.5, color = 'darkgrey') + 
#   scale_y_log10() +
#   scale_x_date(name = "Year (1990-2020)", date_breaks = "5 years", date_labels = "%Y") + 
#   ylab("Expected burn area (ha)") +
#   facet_wrap(. ~ NA_L1NAME, nrow = 2) +
#   theme_classic() + 
#   theme(legend.position = "none")
# ggsave("figures/paper/50yr_returns_level1.pdf", returns_level1_plot, 
#        dpi = 320, width = 15, height = 8)
# 
# quant98_level1 <- quant98 %>% select(c(draw, time, region, yr50)) %>% 
#   left_join(region_key) %>% 
#   left_join(eco_areas) %>% 
#   filter(yr50 != Inf) %>%
#   group_by(NA_L1NAME, time, draw, NA_L1CODE) %>% 
#   summarize(wmean = weighted.mean(yr50, area)) %>% 
#   ungroup() %>% 
#   group_by(time, NA_L1NAME, NA_L1CODE) %>%
#   summarize(med50 = median(wmean), 
#             lower = quantile(wmean, probs = 0.025), 
#             upper = quantile(wmean, probs = 0.975)) %>%
#   ungroup() %>% mutate(NA_L1CODE = factor(NA_L1CODE, levels = levels_l1))
# quant98_level1_plot <- quant98_level1 %>%
#   left_join(time_df) %>%
#   ggplot() + 
#   geom_ribbon(aes(x=date, ymin=lower, ymax=upper, group = NA_L1CODE, fill = NA_L1CODE, alpha = 0.95)) +
#   geom_line(aes(x=date, y=med50, group = NA_L1CODE), linewidth = 0.5, color = 'darkgrey') + 
#   scale_y_log10() +
#   scale_x_date(name = "Year (1990-2020)", date_breaks = "5 years", date_labels = "%Y") + 
#   ylab("Expected burn area (ha) given fire occurrence") +
#   facet_wrap(. ~ NA_L1NAME, nrow = 2) +
#   theme_classic() + 
#   theme(legend.position = "none")
# ggsave("figures/paper/98th_quant_burn_areas_level1.pdf", quant98_level1_plot, 
#        dpi = 320, width = 15, height = 8)
