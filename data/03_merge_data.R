##############################################################################
##  This script merges climate, housing, ecoregion, and fire datasets for   ##
##  use in the Stan models.                                                 ##
##############################################################################

## If the two R scripts below have not already been run, please uncomment and run.
## The csvs produced by these two R scripts are required to run this script.
# source('data/01_get_ecoregion_summaries.R')
# source('data/02_summarize_housing_density.R')

library(tidyverse)
library(sf)
library(zoo)
library(assertthat)
library(lubridate)
library(terra)
source('data/load_eco.R')

ecoregion_shp <- load_ecoregions()

##  Read MTBS fire data -------------------------------------------------------
if (!dir.exists("./data/raw/mtbs_fod_pts_data/")) {
  download.file("https://edcintl.cr.usgs.gov/downloads/sciweb1/shared/MTBS_Fire/data/composite_data/fod_pt_shapefile/mtbs_fod_pts_data.zip",
                destfile = "./data/raw/mtbs_fod_pts_data.zip")
  unzip("./data/raw/mtbs_fod_pts_data.zip",
        exdir = "./data/raw/mtbs_fod_pts_data/")
}

mtbs <- terra::vect(x = './data/raw/mtbs_fod_pts_data/', layer = 'mtbs_FODpoints_DD')
# will lose spat vector aspects if convert to tibble and then manipulate
# using base R to manipulate
mtbs$BurnBndLat <- as.numeric(mtbs$BurnBndLat)
mtbs$BurnBndLon <- as.numeric(mtbs$BurnBndLon)

# only keep lower 48 of US (including DC)
# longitude and latitude bounds can be can be found using conus_bounds.sh
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

# load climate covariate and fire indices data and link to count data frame ----
ecoregion_summaries <- read_csv('./data/processed/ecoregion_summaries.csv') %>%
  mutate(month = match(month, month.abb),
         ym = as.yearmon(paste0(year, "-", month))) %>%
  pivot_wider(names_from = variable, values_from = value)

# compute previous 12 months total precipitation
if (!file.exists('./data/processed/lagged_precip.rds')) {
  lagged_precip <- ecoregion_summaries %>% group_by(NA_L3NAME) %>%
    dplyr::mutate(prev_12mo_precip = rollsumr(pr, k = 12, fill = NA)) %>%
    ungroup() %>%
    dplyr::select(NA_L3NAME, ym, prev_12mo_precip) 
  
  lagged_precip %>%
    write_rds('./data/processed/lagged_precip.rds')
} else {
  lagged_precip <- read_rds('./data/processed/lagged_precip.rds')
}

housing_df <- read_csv('./data/processed/housing_density.csv') %>%
  mutate(month = as.integer(month),
         ym = as.yearmon(paste0(year, "-", month))) %>%
  arrange(NA_L3NAME, year, month, .locale = "en")

ecoregion_summaries <- ecoregion_summaries %>%
  left_join(lagged_precip) %>%
  left_join(housing_df) %>%
  filter(year >= 1990, year <= 2020)

count_df <- left_join(ecoregion_summaries, count_df) %>%
  mutate(er_ym = paste(NA_L3NAME, ym, sep = "_"))

write_csv(count_df, file = './data/processed/count_df.csv')
print('Count, climate, fire indice, and housing data integrated and saved in count_df.csv')
