library(tidyverse)
library(sf)
library(zoo)
library(assertthat)
library(lubridate)
library(terra)

# Albers equal area (AEA) conic projection of North America
# aea_proj <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

source('./full-model/data/process-data/helpers.R')
# Read ecoregion data
ecoregion_shp <- load_ecoregions()

# fix names for Chihuahuan Deserts (L3 ecoregion) 
ecoregion_shp$NA_L3NAME <- as.character(ecoregion_shp$NA_L3NAME)
ecoregion_shp$NA_L3NAME <- ifelse(ecoregion_shp$NA_L3NAME == 'Chihuahuan Desert',
                                  'Chihuahuan Deserts',
                                  ecoregion_shp$NA_L3NAME)

# Read fire data ----------------------
if (!dir.exists("./full-model/data/raw/mtbs_fod_pts_data/")) {
  download.file("https://edcintl.cr.usgs.gov/downloads/sciweb1/shared/MTBS_Fire/data/composite_data/fod_pt_shapefile/mtbs_fod_pts_data.zip",
                destfile = "./full-model/data/raw/mtbs_fod_pts_data.zip")
  unzip("./full-model/data/raw/mtbs_fod_pts_data.zip",
        exdir = "./full-model/data/raw/mtbs_fod_pts_data/")
}

# longitude and latitude bounds are for the lower 48 US, including DC; can be found using lower48_bounds.sh
mtbs <- terra::vect(x = './full-model/data/raw/mtbs_fod_pts_data/', layer = 'mtbs_FODpoints_DD')
# will lose spat vector aspects if convert to tibble and then manipulate
# using base R to manipulate
mtbs$BurnBndLat <- as.numeric(mtbs$BurnBndLat)
mtbs$BurnBndLon <- as.numeric(mtbs$BurnBndLon)
# only keep lower 48 of US (including DC)
idx <- which(mtbs$BurnBndLat < 49.39 & mtbs$BurnBndLat > 24.39 & 
               mtbs$BurnBndLon > -124.848 & mtbs$BurnBndLon < - 66.89 &
               mtbs$BurnBndAc > 1e3 & mtbs$Incid_Type == 'Wildfire')
mtbs <- mtbs[idx,]
mtbs$fire_yr <- year(mtbs$Ig_Date)
mtbs$fire_mon <- month(mtbs$Ig_Date)
mtbs$ym <- as.yearmon(paste0(mtbs$fire_yr, "-", mtbs$fire_mon))
mtbs <- mtbs[which(mtbs$fire_yr >= 1990 & mtbs$fire_yr <= 2020),]

# make ecoregions the same crs as MTBS 
ecoregion_shp <- terra::project(ecoregion_shp, crs(mtbs))

# match each ignition to an ecoregion
mtbs_er <- terra::intersect(mtbs, ecoregion_shp)

unique_er_yms <- expand.grid(
  NA_L3NAME = unique(ecoregion_shp$NA_L3NAME),
  fire_yr = unique(mtbs$fire_yr),
  fire_mon = unique(mtbs$fire_mon)) %>%
  as_tibble()

# count the number of fires in each ecoregion in each month
count_df <- mtbs_er %>%
  as_tibble() %>%
  group_by(NA_L3NAME, fire_yr, fire_mon) %>%
  summarize(n_fire = n()) %>%
  ungroup() %>%
  full_join(unique_er_yms) %>%
  mutate(n_fire = ifelse(is.na(n_fire), 0, n_fire),
         ym = as.yearmon(paste0(fire_yr, "-", fire_mon))) %>%
  arrange(ym)

assert_that(0 == sum(is.na(count_df$NA_L3NAME)))
assert_that(sum(count_df$n_fire) == nrow(mtbs_er))
assert_that(all(ecoregion_shp$NA_L3NAME %in% count_df$NA_L3NAME))

# load climate covariate data and link to count data frame ----
ecoregion_summaries <- read_csv('./full-model/data/processed/ecoregion_summaries.csv') %>%
  mutate(month = match(month, month.abb),
         ym = as.yearmon(paste0(year, "-", month))) %>%
  pivot_wider(names_from = variable, values_from = value)

# Compute previous 12 months total precip
if (!file.exists('./full-model/data/processed/lagged_precip.rds')) {
  lagged_precip <- ecoregion_summaries %>% group_by(NA_L3NAME) %>%
    dplyr::mutate(prev_12mo_precip = rollsumr(pr, k = 12, fill = NA)) %>%
    ungroup() %>%
    dplyr::select(NA_L3NAME, ym, prev_12mo_precip) 
  
  lagged_precip %>%
    write_rds('./full-model/data/processed/lagged_precip.rds')
} else {
  lagged_precip <- read_rds('./full-model/data/processed/lagged_precip.rds')
}

housing_df <- read_csv('./full-model/data/processed/housing_density.csv') %>%
  mutate(month = as.integer(month),
         ym = as.yearmon(paste0(year, "-", month))) %>%
  arrange(NA_L3NAME, year, month, .locale = "en")

ecoregion_summaries <- ecoregion_summaries %>%
  left_join(lagged_precip) %>%
  left_join(housing_df) %>%
  filter(year >= 1990, year <= 2020)

count_df_climate <- left_join(ecoregion_summaries, count_df) %>%
  mutate(er_ym = paste(NA_L3NAME, ym, sep = "_"))

write_csv(count_df_climate, file = './full-model/data/processed/count_df.csv')
print('Count, climate, and housing data integrated and saved in count_df.csv')

# load ERC data and link to count data frame
ecoregion_summaries_erc <- read_csv('./full-model/data/processed/ecoregion_summaries_erc.csv') %>%
  mutate(month = match(month, month.abb),
         ym = as.yearmon(paste0(year, "-", month))) %>%
  pivot_wider(names_from = variable, values_from = value)

ecoregion_summaries_erc <- ecoregion_summaries_erc %>%
  left_join(housing_df) %>%
  filter(year >= 1990, year <= 2020)

count_df_erc <- left_join(ecoregion_summaries_erc, count_df) %>%
  mutate(er_ym = paste(NA_L3NAME, ym, sep = "_"))

write_csv(count_df_erc, file = './full-model/data/processed/count_df_erc.csv')
print('Count, ERC, and housing data integrated and saved in count_df_erc.csv')

# load FWI data and link to count data frame
ecoregion_summaries_fwi <- read_csv('./full-model/data/processed/ecoregion_summaries_fwi.csv') %>%
  mutate(month = match(month, month.abb),
         ym = as.yearmon(paste0(year, "-", month))) %>%
  pivot_wider(names_from = variable, values_from = value)

ecoregion_summaries_fwi <- ecoregion_summaries_fwi %>%
  left_join(housing_df) %>%
  filter(year >= 1990, year <= 2020)

count_df_fwi <- left_join(ecoregion_summaries_fwi, count_df) %>%
  mutate(er_ym = paste(NA_L3NAME, ym, sep = "_"))

write_csv(count_df_fwi, file = './full-model/data/processed/count_df_fwi.csv')
print('Count, FWI, and housing data integrated and saved in count_df_fwi.csv')
