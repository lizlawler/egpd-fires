# Summarizing housing density at the ecoregion level ----------------------
library(purrr)
library(tidyverse)
library(sf)
library(stars)
library(terra)
library(assertthat)

# before running this script, need to extract housing density by ecoregion in QGIS:
# 1) use "join attribute by location (summary), with input layer as the ecoregion file and 
# comparison layer as the housing file;
# 2) keep predicate for features as "intersect"
# 3) select appropriate housing density field as "field to summarize"
# 4) select "mean" as aggregation method
# 5) save as .shp file
# Note: to avoid any errors, choose "do not filter" as method of dealing with "invalid geometry"

dens_files <- paste0("full-model/data/processed/qgis_extraction/", c("1990", "2000", "2010", "2020"), "/")
dens_layer <- paste0("eco_conus_", c("1990", "2000", "2010", "2020"))
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
    year == 'HUDEN2020_' ~ 2020
    ),
    NA_L3NAME = ifelse(NA_L3NAME == 'Chihuahuan Desert','Chihuahuan Deserts', NA_L3NAME)) %>%
  group_by(NA_L3NAME, year) %>%
  summarize(wmean = weighted.mean(value, Shape_Area)) %>%
  ungroup

geom(extractions[[1]]) %>% group_by(NA_L3NAME)

# Then interpolate for each month and year from 1990 - 2020
# using a simple linear sequence
impute_density <- function(df) {
  year_seq <- min(df$year):max(df$year)
  predict_seq <- seq(min(df$year),
                     max(df$year),
                     length.out = (length(year_seq) - 1) * 12)
  preds <- approx(x = df$year,
         y = df$wmean,
         xout = predict_seq)
  res <- as_tibble(preds) %>%
    rename(t = x, wmean = y) %>%
    mutate(year = floor(t),
           month = rep(1:12, times = length(year_seq) - 1)) %>%
    filter(year < 2030)
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
