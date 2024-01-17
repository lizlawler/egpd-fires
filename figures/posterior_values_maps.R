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
library(sf)
library(colorspace)
library(patchwork)
source("./data/load_eco.R")
source("./figures/sizes_gpd_shape_mle.R")

# load ecoregion shapes
ecoregion_shp <- load_ecoregions()
ecoregion_df <- ecoregion_shp %>% st_as_sf()
levels_l1 <- levels(as.factor(as.numeric(ecoregion_df$NA_L1CODE)))
levels_l2 <- levels(as.factor(as.numeric(ecoregion_df$NA_L2CODE)))
ecoregion_df <- ecoregion_df %>% 
  mutate(NA_L2CODE = factor(NA_L2CODE, levels = levels_l2),
         NA_L1CODE = factor(NA_L1CODE, levels = levels_l1),
         NA_L3CODE = as.factor(NA_L3CODE),
         NA_L1NAME = as.factor(str_to_title(NA_L1NAME)))

# read in L3 ecoregion key
region_key <- readRDS(file = "./data/processed/region_key.rds") %>%
  mutate(NA_L2CODE = factor(NA_L2CODE, levels = levels_l2),
         NA_L1CODE = factor(NA_L1CODE, levels = levels_l1),
         NA_L3CODE = as.factor(NA_L3CODE),
         NA_L1NAME = as.factor(str_to_title(NA_L1NAME)))

# read in posterior values of xi and gamma ---------
gamma <- readRDS("./figures/mcmc_draws/gamma.RDS")
rand_int <- readRDS("./figures/mcmc_draws/rand_int.RDS")

# create map of CONUS by L3 ecoregion filled by median gamma value
gamma_map <- gamma %>% group_by(region) %>% summarize(gamma = median(gamma)) %>% ungroup()
ecoreg_gamma <- ecoregion_df %>%
  left_join(gamma_map %>% left_join(region_key)) 
gamma_post_map <- ecoregion_df %>%
  ggplot() +
  geom_sf(size = .1, fill = 'transparent') +
  geom_sf(data = ecoreg_gamma %>% group_by(NA_L3CODE, gamma) %>% summarise(geometry = st_union(st_set_precision(geometry, 1e8))),
          aes(fill=gamma), alpha = 1, lwd = 0.1, inherit.aes = FALSE) +
  theme_void() + scale_fill_gradient2(name = bquote(gamma~value))
ggsave(filename = "./figures/paper_figures/gamma_map.pdf", 
       plot = gamma_post_map, 
       dpi = 320, 
       width = 6, height = 6)
knitr::plot_crop("./figures/paper_figures/gamma_map.pdf")

# create map of MLE estimates of GPD shape param by L2 ecoregion
shape_mle <- readRDS("./figures/sizes_gpd_shape_mle_df.RDS")
ecoreg_shape_mle <- ecoregion_df %>% left_join(shape_mle)

shape_mle_plot <- ecoregion_df %>% 
  ggplot() +
  geom_sf(data = ecoreg_shape_mle %>% group_by(NA_L2CODE, shape_inc_90) %>% summarise(geometry = st_union(st_set_precision(geometry, 1e8))), 
          aes(fill=shape_inc_90), alpha = 0.6, lwd = 0.2, inherit.aes = FALSE) +
  theme_void()
shape_mle_plot_hcl <- shape_mle_plot + scale_fill_continuous_sequential(palette = 'Mint',
                                                             na.value = "transparent",
                                                             name = bquote(xi~value))
ggsave(filename = "./figures/paper_figures/shape_mle_map.pdf", 
       plot = shape_mle_plot_hcl, 
       width = 6, height = 6,
       dpi = 320)
knitr::plot_crop("./figures/paper_figures/shape_mle_map.pdf")

# create map of posterior medians of GPD shape param, with MLE estimates for comparison
xi_map <- rand_int %>% select(-sigma) %>% group_by(region) %>% summarize(xi = median(xi))
ecoreg_xi <- ecoregion_df %>% 
  left_join(xi_map %>% left_join(region_key))

shape_post_map <- ecoregion_df %>%
  ggplot() +
  geom_sf(data = ecoreg_xi %>% group_by(NA_L3CODE, xi) %>% summarise(geometry = st_union(st_set_precision(geometry, 1e8))),
          aes(fill=xi), alpha = 0.6, lwd = 0.2, inherit.aes = FALSE) +
  theme_void()
shape_post_map_hcl <- shape_post_map + scale_fill_continuous_sequential(palette = 'Mint',
                                                                na.value = "transparent",
                                                                name = bquote(xi~value),
                                                                limits = c(0,1.2),
                                                                breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0))  
# make limits of MLE map the same as the posterior map
shape_mle_combo_plot <- shape_mle_plot + scale_fill_continuous_sequential(palette = 'Mint',
                                                                    na.value = "transparent",
                                                                    limits = c(0,1.2)) +
  theme(legend.position = "none")
# combine the maps into one plot
shape_combo_plot <- shape_mle_combo_plot + shape_post_map_hcl + plot_layout(guides = 'collect')
ggsave("./figures/paper_figures/xi_post_with_mle_map.pdf", 
       plot = shape_combo_plot,
       dpi = 320, 
       width = 8, 
       height = 8, 
       bg = 'transparent')
knitr::plot_crop("./figures/paper_figures/xi_post_with_mle_map.pdf")
