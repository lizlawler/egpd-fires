##############################################################################
##  This script produces four sets of maps:                                 ##
##  (1) CONUS split by L3 ecoregion                                         ##
##  (2) CONUS split by L3 ecoregion, colored by L2 ecoregion                ##
##  (3) map from (2), with dark grey borders by L1 ecoregion                ##
##  (4) CONUS colored by L1 ecoregion, with light grey borders by L3        ##
##############################################################################

library(tidyverse)
library(stringr)
library(sf)
library(colorspace)
library(patchwork)
source("./data/load_eco.R")

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

# map (1)
er_map_l3 <- ecoregion_df %>% 
  ggplot() +
  geom_sf(data = ecoregion_df %>% group_by(NA_L3CODE) %>% summarise(geometry = st_union(st_set_precision(geometry, 1e8))),
          fill = "transparent", lwd = 0.3, inherit.aes = FALSE) +
  theme_void() +
  coord_sf(ndiscr = FALSE)
ggsave(filename = "./figures/paper_figures/er_l3_map.pdf", 
       plot = er_map_l3,
       dpi = 320, 
       width = 6, 
       height = 6)
knitr::plot_crop("./figures/paper_figures/er_l3_map.pdf")

# map (2)
er_map_l2 <- ecoregion_df %>%
  ggplot() +
  geom_sf(data = ecoregion_df %>% group_by(NA_L2CODE) %>% summarise(geometry = st_union(st_set_precision(geometry, 1e8))),
          aes(fill = NA_L2CODE), alpha = 0.8, lwd = 0.45, inherit.aes = FALSE, show.legend = FALSE) +
  geom_sf(data = ecoregion_df %>% group_by(NA_L3CODE) %>% summarise(geometry = st_union(st_set_precision(geometry, 1e8))),
          fill = "transparent", lwd = 0.1, inherit.aes = FALSE) +
  theme_void() +
  coord_sf(ndiscr = FALSE)
ggsave(filename = "./figures/paper_figures/er_l2_map.pdf", 
       plot = er_map_l2,
       dpi = 320, 
       width = 6, 
       height = 6)
knitr::plot_crop("./figures/paper_figures/er_l2_map.pdf")

# map (3)
er_map_l1 <- ecoregion_df %>%
  ggplot() +
  geom_sf(data = ecoregion_df %>% group_by(NA_L2CODE) %>% summarise(geometry = st_union(st_set_precision(geometry, 1e8))),
          aes(fill = NA_L2CODE), alpha = 0.8, lwd = 0.25, inherit.aes = FALSE, show.legend = FALSE) +
  geom_sf(data = ecoregion_df %>% group_by(NA_L1CODE) %>% summarise(geometry = st_union(st_set_precision(geometry, 1e8))),
          fill = "transparent", lwd = 0.8, color = "gray20", inherit.aes = FALSE, show.legend = FALSE) +
  geom_sf(data = ecoregion_df %>% group_by(NA_L3CODE) %>% summarise(geometry = st_union(st_set_precision(geometry, 1e8))),
          fill = "transparent", lwd = 0.1, inherit.aes = FALSE) +
  theme_void() +
  coord_sf(ndiscr = FALSE)
ggsave(filename = "./figures/paper_figures/er_l1_map.pdf", 
       plot = er_map_l1,
       dpi = 320, 
       width = 6, 
       height = 6)
knitr::plot_crop("./figures/paper_figures/er_l1_map.pdf")

# map (4)
l1_only_map <- ecoregion_df %>%
  ggplot() +
  geom_sf(data = ecoregion_df %>% group_by(NA_L1CODE) %>% summarise(geometry = st_union(st_set_precision(geometry, 1e8))),
          aes(fill = NA_L1CODE), alpha = 0.8, lwd = 0.5, inherit.aes = FALSE, show.legend = FALSE) +
  geom_sf(data = ecoregion_df %>% group_by(NA_L3CODE) %>% summarise(geometry = st_union(st_set_precision(geometry, 1e8))),
          fill = "transparent", lwd = 0.1, inherit.aes = FALSE) +
  theme_void() +
  coord_sf(ndiscr = FALSE)
ggsave(filename = "./figures/paper_figures/level1_map.pdf", 
       plot = l1_only_map,
       dpi = 320, 
       width = 6, 
       height = 6)
knitr::plot_crop("./figures/paper_figures/level1_map.pdf")
