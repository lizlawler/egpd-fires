##############################################################################
##  This script produces three sets of maps:                                ##
##  (1) MLE estimates of GPD shape parameter for burned areas               ##
##  (2) MLE estimates vs posterior medians of EGPD shape parameter          ##
##  (3) posterior medians of modulating parameter (gamma)                   ##
##############################################################################

library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(tidyverse)
library(stringr)
library(posterior)
# library(lubridate)
source("./data/load_eco.R")

ecoregion_shp <- load_ecoregions()

# read in posterior values of xi and gamma ---------
gamma <- readRDS("./figures/mcmc_draws/gamma.RDS")
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
ggsave(filename = "full-model/figures/paper/gamma_map.png", plot = p, 
       dpi = 320, 
       width = 6, height = 6)
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
  geom_sf(data = eco_eda %>% group_by(NA_L2CODE, shape_inc_90) %>% summarise(geometry = st_union(st_set_precision(geometry, 1e8))), 
          aes(fill=shape_inc_90), alpha = 0.6, lwd = 0.1, inherit.aes = FALSE) +
  theme_void()
eda_90_hcl <- eda_90_plot + scale_fill_continuous_sequential(palette = 'Mint',
                                                             na.value = "transparent",
                                                             name = bquote(xi~value))
# theme(legend.position = c(0.95, 0.4), 
#       legend.key.size = unit(1, "cm"),
#       legend.text = element_text(size = 12),
#       legend.title = element_text(size = 15))
ggsave(filename = "full-model/figures/paper/xi_eda_90th-quant.pdf", 
       plot = eda_90_hcl, 
       width = 6, height = 6,
       dpi = 320)
knitr::plot_crop("full-model/figures/paper/xi_eda_90th-quant.pdf")

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
                                                                breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))  
# theme(legend.position = c(0.9, 0.4), 
#       legend.key.size = unit(1.2, "cm"),
#       legend.text = element_text(size = 11),
#       legend.title = element_text(size = 13))
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
ggsave("full-model/figures/paper/xi_map_with-eda90.png", dpi = 320, width = 8, height = 8, bg = 'transparent')
knitr::plot_crop("full-model/figures/paper/xi_map_with-eda90.png")