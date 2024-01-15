##############################################################################
##  This script downloads individual daily climate covariates by year from  ##
##  gridMET and aggregates them to a monthly timestep.                      ##
##############################################################################

library(lubridate)
library(tidyverse)
library(parallel)
library(pbapply)
library(terra)

if (!dir.exists("data/raw/cb_2021_us_nation_20m/")) {
  download.file("https://www2.census.gov/geo/tiger/GENZ2021/shp/cb_2021_us_nation_20m.zip",
                destfile = "data/raw/cb_2021_us_nation_20m.zip")
  unzip("data/raw/cb_2021_us_nation_20m.zip",
        exdir = "data/raw/cb_2021_us_nation_20m/")
}

## This function takes a climate data url (or erc data url) and a masking 
## spatial polygon as input and makes a monthly summary of the climate data, 
## masked by the polygon
summarize_by_month <- function(url) {
  usa_shp <- terra::vect("data/raw/cb_2021_us_nation_20m",
                         layer = "cb_2021_us_nation_20m")
  file <- basename(url)

  # generate an output filename
  nc_name <- gsub(file,
                  pattern = basename(file),
                  replacement = paste0("monthly_",
                                       basename(file)))
  out_name <- gsub(x = nc_name, pattern = ".nc", replacement = ".tif")

  # determine which climate variable and which year we have
  file_split <- unlist(strsplit(basename(url), split = "_"))
  var <- file_split[1]
  year <- substr(file_split[2], start = 1, stop = 4)

  # pass over file if necessary
  if (out_name %in% list.files(pattern = basename(out_name), recursive = TRUE)) {
    return(paste("File", out_name, "already exists"))
  }
  if (as.numeric(year) < 1983) {
    return("Year outside range of consideration")
  }

  # download the original data file
  download.file(url = url, destfile = file, timeout = 300)

  # determine which function to use to aggregate
  if (var == "pr") {
    fun <- sum
  } else {
    fun <- mean
  }

  # generate monthly summary
  raster <- terra::rast(file)
  start_date <- as.Date(paste(year, "01", "01", sep = "-"))
  end_date <- as.Date(paste(year, "12", "31", sep = "-"))
  date_seq <- seq(start_date, end_date, by = "1 day")
  date_seq <- date_seq[1:terra::nlyr(raster)]
  month_seq <- lubridate::month(date_seq)

  res <- terra::tapp(raster, month_seq, fun = fun, na.rm = TRUE)
  names(res) <- paste(var, year,
                      unique(lubridate::month(date_seq, label = TRUE)),
                      sep = "_")
  new_crs <- terra::crs(res)
  mask_shp <- terra::project(usa_shp, new_crs)
  masked_res <- terra::mask(res, mask_shp)
  terra::writeRaster(masked_res, out_name, filetype = "GTiff")
  unlink(file)
  print(paste0("File ", out_name, " written"))
}

## Summarize daily climate data by month ---------------------------
climate_data_urls <- read.csv('data/raw/climate_data.csv',
                              stringsAsFactors = FALSE)
print('Aggregating daily weather data to monthly means. May take a while...')
pboptions(type = 'timer', use_lb = TRUE)
cl <- makeCluster(getOption("cl.cores", detectCores() / 2), outfile="")
clusterEvalQ(cl, c(library("terra"), "summarize_by_month", options(timeout = 300)))
pblapply(X = climate_data_urls$url,
         FUN = summarize_by_month,
         cl = cl)
stopCluster(cl)

## Summarize daily ERC data by month ---------------------------
erc_data_urls <- read.csv('data/raw/erc_data.csv',
                          stringsAsFactors = FALSE)
print('Aggregating daily ERC data to monthly means. May take a while...')
pboptions(type = 'timer', use_lb = TRUE)
cl <- makeCluster(getOption("cl.cores", detectCores() / 2), outfile="")
clusterEvalQ(cl, c(library("terra"), "summarize_by_month", options(timeout = 300)))
pblapply(X = erc_data_urls$url,
         FUN = summarize_by_month,
         cl = cl)
stopCluster(cl)

## Move files to unique directories, by climate variable ----------------------
tifs_in_home_dir <- list.files(pattern = ".tif") %>%
  tibble(filename = .) %>%
  separate(filename, into = c("time_interval", "variable", "year"), sep = "_") %>%
  mutate(current_name = paste(time_interval, variable, year, sep = "_"))

variables <- distinct(tifs_in_home_dir, variable) %>%
  unlist()

dirs_to_make <- file.path("data", "processed", "climate_data", variables)
sapply(dirs_to_make, dir.create, recursive = TRUE)

# move each tif file to the proper location
tifs_in_home_dir %>%
  mutate(dest_dir = file.path("data", "processed", 'climate_data', variable),
         dest_file = file.path(dest_dir, current_name)) %>%
  do(file.rename(.$current_name, .$dest_file) %>% tibble)

print('Daily climate data aggregated to monthly summaries.')

