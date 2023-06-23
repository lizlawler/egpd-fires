# Summarizing housing density at the ecoregion level ----------------------
library(purrr)
library(tidyverse)
library(sf)
library(stars)
library(terra)
library(assertthat)

#####################################################################################

# before running this script, need to extract housing density by ecoregion in QGIS:
# 1) use "join attribute by location (summary), with input layer as the ecoregion file and 
# comparison layer as the housing file;
# 2) keep predicate for features as "intersect"
# 3) select appropriate housing density field as "field to summarize"
# 4) select "mean" as aggregation method
# 5) save as .shp file
# Note: to avoid any errors, choose "do not filter" as method of dealing with "invalid geometry"
# code to run from the python console within QGIS is as follows for 'HUDEN1990':
# processing.run("qgis:joinbylocationsummary", 
  # {'INPUT':QgsProcessingFeatureSourceDefinition('/Users/lizlawler/Desktop/research/egpd-fires/full-model/data/raw/us_eco_l3/us_eco_l3.shp', 
  # selectedFeaturesOnly=False, featureLimit=-1, flags=QgsProcessingFeatureSourceDefinition.FlagOverrideDefaultGeometryCheck, 
  # geometryCheck=QgsFeatureRequest.GeometryNoCheck),'PREDICATE':[0],
  # 'JOIN':QgsProcessingFeatureSourceDefinition('/Users/lizlawler/Desktop/research/egpd-fires/full-model/data/raw/conus_wui_blk20.gdb|layername=CONUS_WUI_block_1990_2020_change', 
  # selectedFeaturesOnly=False, featureLimit=-1, flags=QgsProcessingFeatureSourceDefinition.FlagOverrideDefaultGeometryCheck, 
  # geometryCheck=QgsFeatureRequest.GeometryNoCheck),'JOIN_FIELDS':['HUDEN1990'],'SUMMARIES':[6],'DISCARD_NONMATCHING':False,
  # 'OUTPUT':'/Users/lizlawler/Desktop/research/egpd-fires/full-model/data/processed/housing-density/huden1990_byeco.shp'})

dens_files <- paste0("full-model/data/processed/housing-density/", c("1990", "2000", "2010", "2020"), "/")
dens_layer <- paste0("huden", c("1990", "2000", "2010", "2020"))
extractions <- mapply(function(z, y) vect(x = z, layer = y), z = dens_files, y = dens_layer)

extraction_df <- lapply(
  extractions, function(x) pivot_longer(
    (data.frame(x) %>% as_tibble() %>% 
       select(c("NA_L3NAME", "Shape_Area"), contains("HUDEN"))), 
    !c("NA_L3NAME", "Shape_Area"), names_to = "year", values_to = "value")) %>% 
  bind_rows %>%
  filter(!is.na(value)) %>%
  mutate(year = case_when(
    year == 'HUDEN1990_' ~ 1990, 
    year == 'HUDEN2000_' ~ 2000, 
    year == 'HUDEN2010_' ~ 2010, 
    year == 'HUDEN2020_' ~ 2021 # changing to 2021 so that all of 2020 is reflected in interpolation below
    ),
    NA_L3NAME = ifelse(NA_L3NAME == 'Chihuahuan Desert','Chihuahuan Deserts', NA_L3NAME)) %>%
  group_by(NA_L3NAME, year) %>%
  summarize(wmean = weighted.mean(value, Shape_Area)) %>%
  ungroup()
# 
# ecoregions <- read_rds(file = "./sim-study/shared-data/ecoregions.RDS")
# 
# ecoregions_geom <- ecoregions %>% filter(!NA_L2NAME == "UPPER GILA MOUNTAINS (?)") %>%
#   mutate(NA_L3NAME = ifelse(NA_L3NAME == 'Chihuahuan Desert','Chihuahuan Deserts', NA_L3NAME))

# joined_hdens <- right_join(ecoregions_geom, extraction_df)
# breaks <- classIntervals(c(min(joined_hdens$wmean) - .00001, joined_hdens$wmean), n = 7)
# joined_hdens <- joined_hdens %>% mutate(wmean_cat = cut(wmean, unique(breaks$brks)))
# 
# dens_map_year <- ecoregions_geom %>%
#   ggplot() +
#   geom_sf(size = .1, fill = 'white') +
#   geom_sf(data = joined_hdens, aes(fill=wmean_cat), alpha = 0.6, lwd = 0, inherit.aes = FALSE) +
#   facet_wrap(year ~ .) +
#   theme(panel.grid.major = element_line(colour = "lightgrey"))
# ggsave("dens_er_map.png", dpi = 300, type = "cairo")

# Then interpolate for each month and year from 1990 - 2020 (so the last full year is 2019)
# using a simple linear sequence
impute_density <- function(df) {
  # year_seq <- min(df$year):max(df$year)
  # predict_seq <- seq(min(df$year),
  #                    max(df$year),
  #                    length.out = (length(year_seq) - 1) * 12)
  year_seq <- sort(rep.int(1990:2020, times = 12))
  time_int <- rep(seq(0.1,0.99,length.out = 12), time = length(1990:2020))
  predict_seq <- year_seq + time_int
  preds <- approx(x = df$year,
         y = df$wmean,
         xout = predict_seq)
  res <- as_tibble(preds) %>%
    rename(t = x, wmean = y) %>%
    mutate(year = floor(t),
           month = rep(1:12, times = length(1990:2020))) %>%
    filter(year < 2021)
  res$NA_L3NAME <- unique(df$NA_L3NAME)
  res
}

monthly_dens <- extraction_df %>%
  split(.$NA_L3NAME) %>%
  map(~impute_density(.)) %>%
  bind_rows %>%
  rename(housing_density = wmean)

out_file <- 'full-model/data/processed/housing_density.csv'

monthly_dens %>%
  write_csv(out_file)

if (file.exists(out_file)) {
  print(paste(out_file, 'successfully written'))
}
