library(tidyverse)
library(lubridate)
#library(rstan)
library(splines)
library(spdep)
library(extRemes)
library(sf)
library(RColorBrewer)
library(patchwork)
library(classInt)
library(assertthat)
library(spatialreg)


# nb <- read_rds('./sim-study/shared-data/nb.rds')
ecoregions <- read_rds(file = "ecoregions.RDS")

ecoregions_geom <- ecoregions %>% filter(!NA_L2NAME == "UPPER GILA MOUNTAINS (?)")

library(posterior)
xi_vals <- bestmod_xivals$post_warmup_draws %>% as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>% 
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  mutate(region = readr::parse_number(name)) %>%
  select(-name) %>% 
  mutate(exp_val = exp(value))

xi_vals_med <- xi_vals %>% group_by(region) %>% summarize(med_val = median(exp_val)) %>% ungroup()
full_reg_key <- as_tibble(region_key) %>%
  mutate(region = c(1:84),
         NA_L2CODE = as.factor(NA_L2CODE),
         NA_L1CODE = as.factor(NA_L1CODE),
         NA_L3CODE = as.factor(NA_L3CODE))
xi_vals_med_regions <- xi_vals_med %>% left_join(full_reg_key)

eco_xi <- ecoregions_geom %>% 
  mutate(NA_L2CODE = as.factor(NA_L2CODE), 
         NA_L1CODE = as.factor(NA_L1CODE), 
         NA_L3CODE = as.factor(NA_L3CODE)) %>% 
  left_join(xi_vals_med_regions)

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
ggsave("full-model/figures/plots_eva2023/xi_vals.png", dpi = 320, bg ='white')

xi_vals <- bestmod_xivals_rerun$post_warmup_draws %>% as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>% 
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  mutate(region = readr::parse_number(name)) %>%
  select(-name) %>% 
  mutate(exp_val = exp(value))

xi_vals_med <- xi_vals %>% group_by(region) %>% summarize(med_val = median(exp_val)) %>% ungroup()

full_reg_key <- as_tibble(region_key_new) %>%
  mutate(region = c(1:84),
         NA_L2CODE = as.factor(NA_L2CODE),
         NA_L1CODE = as.factor(NA_L1CODE),
         NA_L3CODE = as.factor(NA_L3CODE))
xi_vals_med_regions <- xi_vals_med %>% left_join(full_reg_key)

eco_xi <- ecoregions_geom %>% 
  mutate(NA_L2CODE = as.factor(NA_L2CODE), 
         NA_L1CODE = as.factor(NA_L1CODE), 
         NA_L3CODE = as.factor(NA_L3CODE)) %>% 
  left_join(xi_vals_med_regions)

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
ggsave("full-model/figures/plots_eva2023/xi_vals_newdata.png", dpi = 320, bg ='white')


er_map_l1 <- ecoregions_geom %>%
  ggplot() +
  geom_sf(size = .2, fill = "white") +
  geom_sf(data = ecoregions_geom,
          aes(fill = NA_L2CODE), alpha = 0.6, lwd = 0, inherit.aes = FALSE, show.legend = FALSE) +
  geom_sf(data = ecoregions_geom %>% group_by(NA_L1CODE) %>% summarise(),
          fill = "transparent", lwd = 1, color = "gray20", inherit.aes = FALSE, show.legend = FALSE) +
  theme_void() +
  coord_sf(ndiscr = FALSE)
ggsave("er_map_l1.png", dpi = 320)

er_map_l3 <- ecoregions_geom %>%
  ggplot() +
  geom_sf(size = .2, fill = "white", lwd = .3) +
  theme_void() +
  coord_sf(ndiscr = FALSE)
ggsave("er_map_l3.png", dpi = 320)

er_map_l2 <- ecoregions_geom %>%
  ggplot() +
  geom_sf(size = .2, fill = "white") +
  geom_sf(data = ecoregions_geom,
          aes(fill = NA_L2CODE), alpha = 0.6, lwd = .3, inherit.aes = FALSE, show.legend = FALSE) +
  theme_void() +
  coord_sf(ndiscr = FALSE)
ggsave("er_map_l2.png", dpi = 320)

library(maps)
states <- map_data("state")
state_map <- ggplot(data = states) +
  geom_polygon(aes(x = long, y = lat, fill = region, group = group), color = "gray") +
  theme_void() +
  coord_fixed(1.3) +
  guides(fill = FALSE)
ggsave("usa.png", dpi = 320)

world <- map_data("world")
world_map <- ggplot(data = world, mapping = aes(x = long, y = lat, group = group)) + 
  geom_polygon(color = "gray20", fill = "white") + 
  theme_void() +
  coord_fixed(1.3)
ggsave("world.png", dpi = 320)

ggplot(data = ca_df, mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(1.3) + 
  geom_polygon(color = "black", fill = "gray")


# 
# # histograms #### ---------
# hist(params_l2_95$shape)
# hist(params_l2_95$scale)
# 
# hist(params_l2_97$shape)
# hist(params_l2_97$scale)
# 
# hist(params_l3_95$shape)
# hist(params_l3_95$scale)
# 
# hist(params_l3_97$shape)
# hist(params_l3_97$scale)
# 
# plot(params_l3_97$scale)

## make maps ####
## LEVEL 3 ## --------
pre_phi_times <- phi_mat %>% as_tibble() %>% 
  rename_with(., ~ reg_cols) %>%
  mutate(timepoint = 1:all_of(t)) %>%
  pivot_longer(cols = c(1:all_of(r)), values_to = "phi_val", names_to = "region") %>%
  left_join(., mod_reg_key) %>% mutate(type = "truth")

joined_l3_phi_pre_times <- right_join(ecoregions_geom, pre_phi_times)
breaks <- classIntervals(c(min(joined_l3_phi_pre_times$phi_val) - .00001, joined_l3_phi_pre_times$phi_val), n = 5, style = "quantile")
joined_l3_phi_pre_times <- joined_l3_phi_pre_times %>% mutate(phi_cat = cut(phi_val, unique(breaks$brks)))

pre_phi_l3_onetime <- ecoregions_geom %>%
  ggplot() +
  geom_sf(size = .1, fill = 'white') +
  geom_sf(data = joined_l3_phi_pre_times[joined_l3_phi_pre_times$timepoint == 1, ],
          aes(fill=phi_cat), alpha = 0.6, lwd = 0, inherit.aes = FALSE) +
  theme_minimal() +
  theme(panel.grid.major = element_line(colour = "lightgrey"))
pre_phi_l3_onetime

post_phi <- rstan::extract(egpd_fit, pars = "phi")
median_phi <- apply(post_phi$phi, c(2,3), median)
post_phi_times <- median_phi %>% as_tibble() %>%
  rename_with(., ~ reg_cols) %>%
  mutate(timepoint = 1:all_of(t)) %>%
  pivot_longer(cols = c(1:all_of(r)), values_to = "phi_val", names_to = "region") %>%
  left_join(., mod_reg_key) %>% mutate(type = "sim")

joined_l3_phi_post_times <- right_join(ecoregions_geom, post_phi_times)
joined_l3_phi_post_times <- joined_l3_phi_post_times %>% mutate(phi_cat = cut(phi_val, unique(breaks$brks)))

# post_phi_l3_onetime <- ecoregions_geom %>%
#   ggplot() +
#   geom_sf(size = .1, fill = 'white') + 
#   geom_sf(data = joined_l3_phi_post_times[joined_l3_phi_post_times$timepoint == 1, ], aes(fill=phi_cat), alpha = 0.6, lwd = 0, inherit.aes = FALSE) + 
#   theme_minimal() + 
#   theme(panel.grid.major = element_line(colour = "lightgrey"))
# post_phi_l3_onetime

phi_full <- rbind(pre_phi_times, post_phi_times) %>% mutate(type = factor(type, levels = c("truth", "sim")))
joined_l3_phi_full <- right_join(ecoregions_geom, phi_full)
breaks <- classIntervals(c(min(joined_l3_phi_full$phi_val) - .00001, joined_l3_phi_full$phi_val), n=5, style = "quantile")
joined_l3_phi_full <- joined_l3_phi_full %>% mutate(phi_cat = cut(phi_val, unique(breaks$brks)))

full_phi_l3_onetime <- ecoregions_geom %>%
  ggplot() +
  geom_sf(size = .1, fill = 'white') + 
  geom_sf(data = joined_l3_phi_full[joined_l3_phi_full$timepoint == 99, ], 
          aes(fill=phi_cat), alpha = 0.6, lwd = 0, inherit.aes = FALSE) + 
  facet_grid(type ~ .) +
  theme_minimal() + 
  theme(panel.grid.major = element_line(colour = "lightgrey"))
full_phi_l3_onetime

