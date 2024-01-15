##############################################################################
##  This script extracts and summarizes the climate covariates and fire     ##
##  indices by L3 ecoregion.                                                ##
##############################################################################

## If the two R scripts below have not been run, please uncomment and run.
## The files produced by these two R scripts  are required to run this script.
# source('data/00_aggregate_climate_data.R')
# source('data/00_aggregate_fwi.R')

library(tidyverse)
library(parallel)
library(pbapply)
library(terra)
library(assertthat)
source('data/load_eco.R')

ecoregion_shp <- load_ecoregions() 
# fix names for Chihuahuan Deserts (L3 ecoregion) 
ecoregion_shp$NA_L3NAME <- ifelse(ecoregion_shp$NA_L3NAME == 'Chihuahuan Desert',
                                  'Chihuahuan Deserts',
                                  ecoregion_shp$NA_L3NAME)
# fix names for Upper Gila Mountains (L2 ecoregion) 
ecoregion_shp$NA_L2NAME <- ifelse(ecoregion_shp$NA_L2NAME == 'UPPER GILA MOUNTAINS (?)',
                                  'UPPER GILA MOUNTAINS',
                                  ecoregion_shp$NA_L2NAME)

tifs <- list.files("./data/processed/climate_data/",
                   pattern = ".tif",
                   recursive = TRUE,
                   full.names = TRUE)

# remove any housing density geotiffs that matched the file listing
tifs <- tifs[!grepl('den[0-9]{2}\\.tif', tifs)]
tifs <- tifs[!grepl('daily', tifs)]

# grab ERC and FWI tifs
erc_tifs <- tifs[grepl('erc_', tifs)]
fwi_tifs <- tifs[grepl('fwi_', tifs)]

# grab remaining climate variable tifs
tifs <- setdiff(tifs, c(erc_tifs, fwi_tifs))

# confirm crs is the same across tifs, erc_tifs, and fwi_tifs
same.crs(terra::rast(tifs[1]), terra::rast(fwi_tifs[1]))
same.crs(terra::rast(tifs[1]), terra::rast(erc_tifs[1]))

##  Generate indices from polygons for raster extraction ---------------------
r <- terra::rast(tifs[1])
ecoregion_shp <- terra::project(ecoregion_shp, terra::crs(r))
shp_raster_idx <- terra::cells(r, ecoregion_shp)
shp_raster_idx_list <- split(shp_raster_idx, shp_raster_idx[,1])
names(shp_raster_idx_list) <- ecoregion_shp$NA_L3NAME

# the above list of indices has one element per polygon, but we want one per region
ecoregion_raster_idx <- vector(mode = 'list',
                               length = length(unique(ecoregion_shp$NA_L3NAME)))
ecoregion_names <- sort(unique(ecoregion_shp$NA_L3NAME))
names(ecoregion_raster_idx) <- ecoregion_names

for (i in seq_along(ecoregion_names)) {
  list_elements <- names(shp_raster_idx_list) == names(ecoregion_raster_idx)[i]
  assert_that(any(list_elements))
  ecoregion_raster_idx[[i]] <- shp_raster_idx_list[list_elements] %>%
    unlist() %>% unique()
}

# verify that no cells are duplicated
assert_that(ecoregion_raster_idx %>%
              sapply(FUN = function(x) any(duplicated(x))) %>%
              sum == 0)

# verify that all ecoregions have some cells
assert_that(ecoregion_raster_idx %>%
              sapply(FUN = function(x) length(x)) %>%
              min > 0)

# define an efficient extraction function to get mean values by polygon
fast_extract <- function(rasterfile, index_list) {
  r <- terra::rast(rasterfile)

  polygon_means <- lapply(index_list, function(x) {
    terra::extract(r, x) %>% colMeans(., na.rm = TRUE)})

  list_of_dfs <- lapply(polygon_means, function(x) {
    as.data.frame(x) %>% tibble::rownames_to_column() %>%
      dplyr::rename(var_ym = "rowname")})

  merged_dfs <- dplyr::bind_rows(list_of_dfs, .id = 'NA_L3NAME')
  wide_df <- tidyr::pivot_wider(merged_dfs, names_from = var_ym, values_from = x)
  return(tibble::as_tibble(wide_df))
}

##  Extract monthly climate data by ecoregion ---------------------------------
print('Aggregating monthly climate data to ecoregion means. May take a while...')
pboptions(type = 'timer', use_lb = TRUE)
cl <- makeCluster(getOption("cl.cores", detectCores() / 2), outfile="")
clusterEvalQ(cl, c(library("terra"), "fast_extract", library("tidyverse")))
extractions <- pblapply(tifs, 
                        fast_extract, 
                        index_list = ecoregion_raster_idx, 
                        cl = cl)
stopCluster(cl)

##  Process extracted climate data values into a usable data frame ------------
ecoregion_climate <- lapply(extractions, function(x) pivot_longer(x, !NA_L3NAME, names_to = "variable", values_to = "value")) %>% 
  bind_rows() %>%
  separate(variable, into = c("variable", "year", "month"), sep = "_") %>%
  mutate(year = parse_number(year),
         month = parse_factor(month)) %>%
  arrange(year, month, variable, NA_L3NAME, .locale = "en")
# including '.locale' in `arrange` is necessary, otherwise a few ecoregions are 
# ordered differently due to different treatment of the word 'and' 
# in base R vs tidyverse

##  Extract monthly ERC data by ecoregion -------------------------------------
print('Aggregating monthly ERC data to ecoregion means. May take a while...')
pboptions(type = 'timer', use_lb = TRUE)
cl <- makeCluster(getOption("cl.cores", detectCores() / 2), outfile="")
clusterEvalQ(cl, c(library("terra"), "fast_extract", library("tidyverse")))
extractions <- pblapply(erc_tifs, 
                        fast_extract, 
                        index_list = ecoregion_raster_idx, 
                        cl = cl)
stopCluster(cl)

##  Process extracted ERC values into a usable data frame ---------------------
ecoregion_erc <- lapply(extractions, function(x) pivot_longer(x, !NA_L3NAME, names_to = "variable", values_to = "value")) %>% 
  bind_rows() %>%
  separate(variable, into = c("variable", "year", "month"), sep = "_") %>%
  mutate(year = parse_number(year),
         month = parse_factor(month)) %>%
  arrange(year, month, variable, NA_L3NAME, .locale = "en")

##  Extract monthly FWI data by ecoregion -------------------------------------
print('Aggregating monthly FWI data to ecoregion means. May take a while...')
pboptions(type = 'timer', use_lb = TRUE)
cl <- makeCluster(getOption("cl.cores", detectCores() / 2), outfile="")
clusterEvalQ(cl, c(library("terra"), "fast_extract", library("tidyverse")))
extractions <- pblapply(fwi_tifs, 
                        fast_extract, 
                        index_list = ecoregion_raster_idx, 
                        cl = cl)
stopCluster(cl)

##  Process extracted FWI values into a usable data frame ---------------------
ecoregion_fwi <- lapply(extractions, function(x) pivot_longer(x, !NA_L3NAME, names_to = "variable", values_to = "value")) %>% 
  bind_rows() %>%
  separate(variable, into = c("variable", "year", "month"), sep = "_") %>%
  mutate(year = parse_number(year),
         month = parse_factor(month)) %>%
  arrange(year, month, variable, NA_L3NAME, .locale = "en")

## Combine all three tibbles to create one summary of ecoregions --------------
ecoregion_summaries <- ecoregion_climate %>% rbind(ecoregion_erc) %>% rbind(ecoregion_fwi)
destfile <- "./data/processed/ecoregion_summaries.csv"
write_csv(ecoregion_summaries, destfile)
print(paste('Ecoregion summaries written to', destfile))
