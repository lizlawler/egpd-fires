##############################################################################
##  This script downloads four climate covariates one year at a time,       ##
##  then calculates the daily FWI as a raster, then aggregates to a         ##
##  monthly timepoint.                                                      ##
##############################################################################

library(terra)
library(tidyverse)
library(lubridate)
library(assertthat)
library(progress)
source("data/fwiraster_terra.R")
options(timeout = max(600, getOption("timeout")))

climate_data_urls <- read.csv('data/raw/climate_data.csv',
                              stringsAsFactors = FALSE)

climate_data_df <- climate_data_urls |> as_tibble() |> 
  mutate(baseurl = basename(url),
         baseurl = str_remove(baseurl, ".nc")) |>
  separate_wider_delim(cols = baseurl, delim = "_", names = c("var", "file_year")) |>
  mutate(file_year = as.numeric(file_year)) |> pivot_wider(names_from = "var", values_from = "url")

start_year <- 1983
end_year <- 2021
total_years <- end_year - start_year + 1
pb <- progress_bar$new(
  format = "  [:bar] :percent Estimated Time Remaining: :eta",
  total = total_years,
  clear = FALSE
)

for(year in start_year:end_year) {
  pb$tick()
  # download 4 netCDF files for each year (one for each variable)
  year_urls <- climate_data_df |> filter(file_year == year) |> select(-file_year) |> as.character()
  temp_files <- basename(year_urls)
  file_year <- as.numeric(str_extract(temp_files[1], pattern = "\\d{4}"))
  assert_that(year == file_year)
  print(paste0("Downloading these files to your local machine now: ", temp_files))
  download.file(year_urls, destfile = temp_files, method = "libcurl")
  ndays <- nlyr(rast(temp_files[1]))
  
  # month sequence is used in fwi function; dates are used in monthly aggregation after for loop
  start_date <- as.Date(paste(year, "01", "01", sep = "-"))
  end_date <- as.Date(paste(year, "12", "31", sep = "-"))
  date_seq <- seq(start_date, end_date, by = "1 day")
  date_seq <- date_seq[1:ndays]
  month_seq <- lubridate::month(date_seq)
  
  # need to initialize with the first day of 1983
  if(year == 1983) {
    init_day <- c(terra::rast(temp_files[1], lyrs = 1), terra::rast(temp_files[2], lyrs = 1), 
                  terra::rast(temp_files[3], lyrs = 1), terra::rast(temp_files[4], lyrs = 1))
    names(init_day) <- c("prec", "temp", "ws", "rh")
    fire_indices <- fwiRaster_terra(init_day, mon = month_seq[1], uppercase = FALSE)
    fwi_rast <- subset(fire_indices, "fwi")
    start <- 2
  } else {
    fwi_rast <- rast()
    start <- 1
  }
  
  for(i in start:ndays) {
    prev_day <- fire_indices
    day_all <- c(terra::rast(temp_files[1], lyrs = i), terra::rast(temp_files[2], lyrs = i), 
                 terra::rast(temp_files[3], lyrs = i), terra::rast(temp_files[4], lyrs = i))
    names(day_all) <- c("prec", "temp", "ws", "rh")
    fire_indices <- fwiRaster_terra(day_all, init = prev_day, mon = month_seq[i], uppercase = FALSE)
    terra::add(fwi_rast) <- subset(fire_indices, "fwi")
    print(paste0("Completed day ", i, " of ", year))
  }
  
  print(paste0("Completed daily FWI for ", year, "; now aggregating to monthly means"))
  res <- terra::tapp(fwi_rast, month_seq, fun = mean, na.rm = TRUE)
  names(res) <- paste("fwi", year,
                      unique(lubridate::month(date_seq, label = TRUE)),
                      sep = "_")
  # mask to usa_shp
  rast_crs <- terra::crs(res)
  usa_shp <- terra::vect("data/raw/cb_2021_us_nation_20m",
                         layer = "cb_2021_us_nation_20m")
  mask_shp <- terra::project(usa_shp, rast_crs)
  masked_res <- terra::mask(res, mask_shp)
  out_name <- paste0("data/processed/climate_data/fixed_fwi/", paste("monthly", "fwi", year, sep = "_"), ".tif")
  terra::writeRaster(masked_res, out_name, filetype = "GTiff")
  print(paste0("Monthly FWI means for ", year, " have been written to disk"))
  
  # proceed with next year
  unlink(temp_files)
  rm(fwi_rast)
  gc()
}
