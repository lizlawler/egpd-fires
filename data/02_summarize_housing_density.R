##############################################################################
##  This script summarizes housing density at the L3 ecoregion level.       ##
##############################################################################

## If '01_extract_housing_dens.txt' has not been executed, please start there
## first. The tifs produced by that code are required to proceed with the below.

library(purrr)
library(tidyverse)
library(sf)
library(terra)

dens_files <- paste0("data/processed/housing_density/", c("1990", "2000", "2010", "2020"), "/")
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
    year == 'HUDEN2020_' ~ 2021 # use 2021 instead of 2020 so that all of 2020 is reflected in interpolation below
    ),
    NA_L3NAME = ifelse(NA_L3NAME == 'Chihuahuan Desert','Chihuahuan Deserts', NA_L3NAME)) %>%
  group_by(NA_L3NAME, year) %>%
  summarize(wmean = weighted.mean(value, Shape_Area)) %>%
  ungroup()

## Interpolate for each month and year from 1990 - 2021 (so the last full year is 2020),
## using a simple linear sequence
impute_density <- function(df) {
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

out_file <- 'data/processed/housing_density.csv'

monthly_dens %>%
  write_csv(out_file)

if (file.exists(out_file)) {
  print(paste(out_file, 'successfully written'))
}
