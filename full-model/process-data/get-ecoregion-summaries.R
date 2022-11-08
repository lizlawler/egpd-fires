library(tidyverse)
library(raster)
library(parallel)
library(pbapply)
library(rgdal)
library(assertthat)
source('./full-model/process-data/helpers.R')


# Extracting monthly climate summaries for ecoregions ---------------------
ecoregion_shp <- rgdal::readOGR(dsn = './full-model/data/raw/us_eco_l3', layer = 'us_eco_l3') %>%
  spTransform(CRSobj = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))

ecoregion_shp$NA_L3NAME <- as.character(ecoregion_shp$NA_L3NAME)

ecoregion_shp$NA_L3NAME <- ifelse(ecoregion_shp$NA_L3NAME == 'Chihuahuan Desert',
                                  'Chihuahuan Deserts',
                                  ecoregion_shp$NA_L3NAME)

tifs <- list.files("./full-model/data/processed/climate-data",
                   pattern = ".tif",
                   recursive = TRUE,
                   full.names = TRUE)

# remove any housing density geotiffs that matched the file listing
tifs <- tifs[!grepl('den[0-9]{2}\\.tif', tifs)]

# Generate indices from polygons for raster extraction --------------------
r <- raster::brick(tifs[1])
shp_raster_idx <- cellFromPolygon(r, ecoregion_shp)
names(shp_raster_idx) <- ecoregion_shp$NA_L3NAME

# this list of indices has one element per polygon, but we want one per region
ecoregion_raster_idx <- vector(mode = 'list',
                               length = length(unique(ecoregion_shp$NA_L3NAME)))
ecoregion_names <- sort(unique(ecoregion_shp$NA_L3NAME))
names(ecoregion_raster_idx) <- ecoregion_names
for (i in seq_along(ecoregion_names)) {
  list_elements <- names(shp_raster_idx) == ecoregion_names[i]
  assert_that(any(list_elements))
  ecoregion_raster_idx[[i]] <- shp_raster_idx[list_elements] %>%
    unlist
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
  r <- raster::brick(rasterfile)
  var_yr <- stringr::str_replace(stringr::str_extract(rasterfile, "\\w+.\\d{4}"), "monthly_", "")

  polygon_means <- lapply(index_list, function(x) {
    extracts <- raster::extract(r, x)
    colMeans(extracts, na.rm = TRUE)})

  list_of_dfs <- lapply(polygon_means, function(x) {
    df <- as.data.frame(x)
    df <- tibble::rownames_to_column(df)
    df <- dplyr::mutate(df, rowname = gsub("layer.", paste0(var_yr, "_"), rowname, fixed = TRUE))
    dplyr::rename(df, var_ym = rowname)})

  merged_dfs <- dplyr::bind_rows(list_of_dfs, .id = 'NA_L3NAME')
  wide_df <- tidyr::pivot_wider(merged_dfs, names_from = var_ym, values_from = x)
  return(tibble::as_tibble(wide_df))
}

# Extract climate data ---------------------------------------
print('Aggregating monthly climate data to ecoregion means. May take a while...')
pboptions(type = 'timer', use_lb = TRUE)
cl <- makeCluster(getOption("cl.cores", detectCores() / 2))
extractions <- pblapply(tifs, 
                        fast_extract, 
                        index_list = ecoregion_raster_idx, 
                        cl = cl)
stopCluster(cl)


# Process extracted values into a usable data frame -----------------------
ecoregion_summaries <- lapply(extractions, function(x) pivot_longer(x, !NA_L3NAME, names_to = "variable", values_to = "value")) %>% 
  bind_rows() %>%
  separate(variable, into = c("variable", "year", "month"), sep = "_") %>%
  mutate(year = parse_number(year),
         month = parse_number(month)) %>%
  arrange(year, month, variable, NA_L3NAME)

destfile <- "./full-model/data/processed/ecoregion_summaries.csv"
write_csv(ecoregion_summaries, destfile)

print(paste('Ecoregion climate summaries written to', destfile))
