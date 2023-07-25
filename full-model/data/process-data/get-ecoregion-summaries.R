library(tidyverse)
library(parallel)
library(pbapply)
library(terra)
library(assertthat)
source('./full-model/data/process-data/helpers.R')

# Extracting monthly climate summaries for ecoregions ---------------------
ecoregion_shp <- load_ecoregions() 
ecoregion_shp$NA_L3NAME <- as.character(ecoregion_shp$NA_L3NAME)
ecoregion_shp$NA_L3NAME <- ifelse(ecoregion_shp$NA_L3NAME == 'Chihuahuan Desert',
                                  'Chihuahuan Deserts',
                                  ecoregion_shp$NA_L3NAME)

tifs <- list.files("./full-model/data/processed/climate-data/",
                   pattern = ".tif",
                   recursive = TRUE,
                   full.names = TRUE)

# remove any housing density geotiffs that matched the file listing
tifs <- tifs[!grepl('den[0-9]{2}\\.tif', tifs)]
tifs <- tifs[!grepl('daily', tifs)] # remove daily weather tifs

# grab ERC and FWI tifs only
erc_tifs <- tifs[grepl('erc_', tifs)]
fwi_tifs <- tifs[grepl('fwi_', tifs)]

# grab remaining climate tifs
tifs <- setdiff(tifs, c(erc_tifs, fwi_tifs))

# confirm crs is the same across tifs, erc_tifs, and fwi_tifs
same.crs(terra::rast(tifs[1]), terra::rast(fwi_tifs[1]))
same.crs(terra::rast(tifs[1]), terra::rast(erc_tifs[1]))

# Generate indices from polygons for raster extraction --------------------
r <- terra::rast(tifs[1])
ecoregion_shp <- terra::project(ecoregion_shp, terra::crs(r))
shp_raster_idx <- terra::cells(r, ecoregion_shp)
shp_raster_idx_list <- split(shp_raster_idx, shp_raster_idx[,1])
names(shp_raster_idx_list) <- ecoregion_shp$NA_L3NAME

# this list of indices has one element per polygon, but we want one per region
ecoregion_raster_idx <- vector(mode = 'list',
                               length = length(unique(ecoregion_shp$NA_L3NAME)))
ecoregion_names <- sort(unique(ecoregion_shp$NA_L3NAME))
names(ecoregion_raster_idx) <- ecoregion_names
# the above sorting places "Southern and Baja ..." after "Southeastern Wisconsin..."
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

# Extract climate data ---------------------------------------
print('Aggregating monthly climate data to ecoregion means. May take a while...')
pboptions(type = 'timer', use_lb = TRUE)
cl <- makeCluster(getOption("cl.cores", detectCores() / 2), outfile="")
clusterEvalQ(cl, c(library("terra"), "fast_extract", library("tidyverse")))
extractions <- pblapply(tifs, 
                        fast_extract, 
                        index_list = ecoregion_raster_idx, 
                        cl = cl)
stopCluster(cl)

# Process extracted climate data values into a usable data frame -----------------------
ecoregion_summaries <- lapply(extractions, function(x) pivot_longer(x, !NA_L3NAME, names_to = "variable", values_to = "value")) %>% 
  bind_rows() %>%
  separate(variable, into = c("variable", "year", "month"), sep = "_") %>%
  mutate(year = parse_number(year),
         month = parse_factor(month)) %>%
  arrange(year, month, variable, NA_L3NAME, .locale = "en")

destfile <- "./full-model/data/processed/ecoregion_summaries.csv"
write_csv(ecoregion_summaries, destfile)
print(paste('Ecoregion climate summaries written to', destfile))

# Extract ERC data ---------------------------------------
print('Aggregating monthly ERC data to ecoregion means. May take a while...')
pboptions(type = 'timer', use_lb = TRUE)
cl <- makeCluster(getOption("cl.cores", detectCores() / 2), outfile="")
clusterEvalQ(cl, c(library("terra"), "fast_extract", library("tidyverse")))
extractions <- pblapply(erc_tifs, 
                        fast_extract, 
                        index_list = ecoregion_raster_idx, 
                        cl = cl)
stopCluster(cl)

# Process extracted ERC values into a usable data frame -----------------------
ecoregion_summaries <- lapply(extractions, function(x) pivot_longer(x, !NA_L3NAME, names_to = "variable", values_to = "value")) %>% 
  bind_rows() %>%
  separate(variable, into = c("variable", "year", "month"), sep = "_") %>%
  mutate(year = parse_number(year),
         month = parse_factor(month)) %>%
  arrange(year, month, variable, NA_L3NAME, .locale = "en")

destfile <- "./full-model/data/processed/ecoregion_summaries_erc.csv"
write_csv(ecoregion_summaries, destfile)
print(paste('Ecoregion ERC summaries written to', destfile))

# Extract FWI data ---------------------------------------
print('Aggregating monthly FWI data to ecoregion means. May take a while...')
pboptions(type = 'timer', use_lb = TRUE)
cl <- makeCluster(getOption("cl.cores", detectCores() / 2), outfile="")
clusterEvalQ(cl, c(library("terra"), "fast_extract", library("tidyverse")))
extractions <- pblapply(fwi_tifs, 
                        fast_extract, 
                        index_list = ecoregion_raster_idx, 
                        cl = cl)
stopCluster(cl)

# Process extracted values into a usable data frame -----------------------
ecoregion_summaries <- lapply(extractions, function(x) pivot_longer(x, !NA_L3NAME, names_to = "variable", values_to = "value")) %>% 
  bind_rows() %>%
  separate(variable, into = c("variable", "year", "month"), sep = "_") %>%
  mutate(year = parse_number(year),
         month = parse_factor(month)) %>%
  arrange(year, month, variable, NA_L3NAME, .locale = "en")

destfile <- "./full-model/data/processed/ecoregion_summaries_fwi.csv"
write_csv(ecoregion_summaries, destfile)
print(paste('Ecoregion FWI summaries written to', destfile))
