library(tidyverse)
library(sf)
library(terra)
library(zoo)
library(assertthat)
library(lubridate)

# Albers equal area (AEA) conic projection of North America
# aea_proj <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

# Read fire data ----------------------
if (!dir.exists("./full-model/data/raw/mtbs_fod_pts_data/")) {
  download.file("https://edcintl.cr.usgs.gov/downloads/sciweb1/shared/MTBS_Fire/data/composite_data/fod_pt_shapefile/mtbs_fod_pts_data.zip",
                destfile = "./full-model/data/raw/mtbs_fod_pts_data.zip")
  unzip("./full-model/data/raw/mtbs_fod_pts_data.zip",
        exdir = "./full-model/data/raw/mtbs_fod_pts_data/")
}

# longitude and latitude bounds are for the lower 48 US, including DC; can be found using lower48_bounds.sh
mtbs <- st_read('./full-model/data/raw/mtbs_fod_pts_data/mtbs_FODpoints_DD.shp')
test_eco <- terra::project(vect(ecoregions), crs(mtbs))
mtbs <- st_read('./full-model/data/raw/mtbs_fod_pts_data/mtbs_FODpoints_DD.shp') %>%
  mutate(BurnBndLat = as.numeric(BurnBndLat),
         BurnBndLon = as.numeric(BurnBndLon)) %>%
  filter(BurnBndLat < 49.39, BurnBndLat > 24.39, 
         BurnBndLon > -124.848, BurnBndLon < -66.89, 
         BurnBndAc > 1e3, 
         Incid_Type == 'Wildfire') %>%
  st_transform(st_crs(ecoregions)) %>%
  mutate(ym = as.yearmon(Ig_Date), 
         FIRE_YEAR = year(Ig_Date), 
         FIRE_MON = month(Ig_Date))

# Read ecoregion data -----
ecoregions <- vect('./full-model/data/raw/us_eco_l3/us_eco_l3.shp') %>% project(., crs(mtbs))

# fix names for chihuahuan desert
ecoregions <- ecoregions %>%
  mutate(NA_L3NAME = as.character(NA_L3NAME),
         NA_L3NAME = ifelse(NA_L3NAME == 'Chihuahuan Desert',
                            'Chihuahuan Deserts',
                            NA_L3NAME))


# match each ignition to an ecoregion
if (!file.exists("./full-model/data/processed/ov.rds")) {
  st_over <- function(x, y) {
    sapply(st_intersects(x,y), function(z) if (length(z)==0) NA_integer_ else z[1])
  }
  ov <- st_over(mtbs, ecoregions)
  write_rds(ov, "./full-model/data/processed/ov.rds")
}
ov <- read_rds("./full-model/data/processed/ov.rds")

mtbs <- mtbs %>%
  mutate(NA_L3NAME = ecoregions$NA_L3NAME[ov]) %>%
  filter(!is.na(NA_L3NAME))

unique_er_yms <- expand.grid(
  NA_L3NAME = unique(ecoregions$NA_L3NAME),
  FIRE_YEAR = unique(mtbs$FIRE_YEAR),
  FIRE_MON = unique(mtbs$FIRE_MON)) %>%
  as_tibble

# count the number of fires in each ecoregion in each month
count_df <- mtbs %>%
  as_tibble() %>%
  dplyr::select(-geometry) %>%
  group_by(NA_L3NAME, FIRE_YEAR, FIRE_MON) %>%
  summarize(n_fire = n()) %>%
  ungroup %>%
  full_join(unique_er_yms) %>%
  mutate(n_fire = ifelse(is.na(n_fire), 0, n_fire),
         ym = as.yearmon(paste(FIRE_YEAR, sprintf("%02d", FIRE_MON), sep = "-"))) %>%
  arrange(ym) %>%
  filter(ym >= 'Jan 1984') # first record is feb 1984 in mtbs data

assert_that(0 == sum(is.na(count_df$NA_L3NAME)))
assert_that(sum(count_df$n_fire) == nrow(mtbs))
assert_that(all(ecoregions$NA_L3NAME %in% count_df$NA_L3NAME))

# load covariate data and link to count data frame
ecoregion_summaries <- read_csv('./full-model/data/processed/ecoregion_summaries.csv') %>%
  mutate(ym = as.yearmon(paste(year,
                 sprintf("%02d", month),
                 sep = "-"))) %>%
  pivot_wider(names_from = variable, values_from = value)

# Compute previous 12 months total precip
if (!file.exists('./full-model/data/processed/lagged_precip.rds')) {
  ecoregion_summaries$prev_12mo_precip <- NA
  pb <- txtProgressBar(max = nrow(ecoregion_summaries), style = 3)
  for (i in 1:nrow(ecoregion_summaries)) {
    setTxtProgressBar(pb, i)
    if (ecoregion_summaries$year[i] > 1983) {
      start_ym <- ecoregion_summaries$ym[i] - 1 # minus one year

      ecoregion_summaries$prev_12mo_precip[i] <- ecoregion_summaries %>%
        filter(NA_L3NAME == ecoregion_summaries$NA_L3NAME[i],
               ym >= start_ym,
               ym <= ecoregion_summaries$ym[i]) %>%
        summarize(twelve_month_total = sum(pr)) %>%
        c %>%
        unlist
    }
  }

  ecoregion_summaries %>%
    dplyr::select(NA_L3NAME, ym, prev_12mo_precip) %>%
    write_rds('./full-model/data/processed/lagged_precip.rds')
}

housing_df <- read_csv('./full-model/data/processed/housing_density.csv') %>%
  mutate(ym = as.yearmon(paste(year, sprintf("%02d", month), sep = "-"))) %>%
  arrange(NA_L3NAME, year, month)

ecoregion_summaries <- ecoregion_summaries %>%
  left_join(read_rds('./full-model/data/processed/lagged_precip.rds')) %>%
  left_join(housing_df) %>%
  filter(ym >= 'Jan 1984')

count_df <- left_join(count_df, ecoregion_summaries) %>%
  mutate(er_ym = paste(NA_L3NAME, ym, sep = "_"))

write_csv(count_df, file = './full-model/data/processed/count_df.csv')

print('Count, climate, and housing data integrated and saved in count_df.csv')
