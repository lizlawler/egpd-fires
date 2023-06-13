library(tidyverse)
library(parallel)
library(pbapply)
library(terra)
library(assertthat)
source('./full-model/data/process-data/helpers.R')

# Extracting monthly climate summaries for ecoregions ---------------------
ecoregion_shp <- load_ecoregions()
# %>%
#   project("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

ecoregion_shp$NA_L3NAME <- as.character(ecoregion_shp$NA_L3NAME)
ecoregion_shp$NA_L3NAME <- ifelse(ecoregion_shp$NA_L3NAME == 'Chihuahuan Desert',
                                  'Chihuahuan Deserts',
                                  ecoregion_shp$NA_L3NAME)

# tifs <- list.files("./full-model/data/processed/climate-data",
#                    pattern = ".tif",
#                    recursive = TRUE,
#                    full.names = TRUE)

daily_tifs <- list.files("./full-model/data/processed/climate-data/daily",
                   pattern = ".tif",
                   recursive = TRUE,
                   full.names = TRUE)

# remove any housing density geotiffs that matched the file listing
# tifs <- tifs[!grepl('den[0-9]{2}\\.tif', tifs)]
# 
# # grab ERC tifs only
# erc_tifs <- tifs[grepl('erc_', tifs)]

# Generate indices from polygons for raster extraction --------------------
r <- terra::rast(daily_tifs[1])
tif_crs <- terra::crs(r)
ecoregion_shp <- terra::project(ecoregion_shp, tif_crs)
eco_names <- ecoregion_shp$NA_L3NAME %>% as_tibble() %>% rowid_to_column(var = "id") 
shp_raster_idx_df <- terra::cells(r, ecoregion_shp) %>% as_tibble() %>%
  rename(id = ID) %>% left_join(eco_names)

# the above tibble of cell indices has one cell per polygon, but we want one per region (ie no duplicates)
ecoregion_raster_idx <- vector(mode = 'list',
                               length = length(unique(ecoregion_shp$NA_L3NAME)))
ecoregion_names <- sort(unique(ecoregion_shp$NA_L3NAME))
names(ecoregion_raster_idx) <- ecoregion_names
for (i in seq_along(ecoregion_names)) {
  df_idx <- shp_raster_idx_df$value == ecoregion_names[i]
  assert_that(any(df_idx))
  ecoregion_raster_idx[[i]] <- unique(shp_raster_idx_df[df_idx, 2])
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
# fast_extract_monthly <- function(rasterfile, index_list) {
#   r <- terra::rast(rasterfile)
# 
#   polygon_means <- lapply(index_list, function(x) {
#     extracts <- terra::extract(r, x$cell)
#     colMeans(extracts, na.rm = TRUE)})
# 
#   list_of_dfs <- lapply(polygon_means, function(x) {
#     df <- as.data.frame(x)
#     colnames(df) <- "mean_by_er"
#     df <- tibble::rownames_to_column(df) %>% 
#       dplyr::rename(var_ymd = rowname)})
# 
#   merged_dfs <- dplyr::bind_rows(list_of_dfs, .id = 'NA_L3NAME')
#   wide_df <- tidyr::pivot_wider(merged_dfs, names_from = var_ymd, values_from = mean_by_er)
#   return(wide_df)
# }
# the below function creates csvs of the extracted raster values by ecoregion and by variable and year
# takes about 2.5 hours on Macbook pro with M1 chip, 32 GB memory
fast_extract_daily <- function(rasterfile, index_list) {
  r <- terra::rast(rasterfile)
  temp_name <- gsub(x = rasterfile,
                  pattern = basename(rasterfile),
                  replacement = paste0("by_er_", basename(rasterfile)))
  out_name <- gsub(x = temp_name, pattern = ".tif", replacement = ".csv")
  
  if (out_name %in% list.files(pattern = basename(out_name), recursive = TRUE)) {
    return(paste("File", out_name, "already exists"))
  }

  extracts <- lapply(index_list, function(x) terra::extract(r, x$cell))
  
  list_of_dfs <- lapply(extracts, function(x) {
    df <- tibble::as_tibble(x)
    df <- tibble::rowid_to_column(df) %>% 
      dplyr::rename(id = rowid) %>%
      tidyr::pivot_longer(!id, names_to = "var_ymd", values_to = "value") %>%
      dplyr::filter(!is.na(value))})
  
  merged_dfs <- dplyr::bind_rows(list_of_dfs, .id = 'NA_L3NAME')
  wide_df <- tidyr::pivot_wider(merged_dfs, names_from = var_ymd, values_from = value)
  write_csv(wide_df, file = out_name)
  print(paste0("File ", basename(out_name), " written"))
  gc()
}

print('Aggregating daily climate data to ecoregion means. May take a while...')
pboptions(type = 'timer', use_lb = TRUE)
cl <- makeCluster(getOption("cl.cores", detectCores() / 2), outfile="")
clusterEvalQ(cl, c(library("terra"), library("tidyverse"), "ecoregion_raster_idx", "fast_extract_daily"))
extractions <- pblapply(daily_tifs, 
                        fast_extract_daily, 
                        index_list = ecoregion_raster_idx, 
                        cl = cl)
stopCluster(cl)

# now create ecoregion specific csvs
longer_daily <- by_er_daily_pr_2015 %>% pivot_longer(cols = !c("NA_L3NAME", "id"), names_to = "combo", values_to = "value")
daily_rast_er <- longer_daily %>% mutate(rast_er = paste0())

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
daily_climate_eco <- lapply(extractions, function(x) pivot_longer(x, !NA_L3NAME, names_to = "var_ymd", values_to = "value")) %>% 
  bind_rows() %>%
  separate(var_ymd, into = c("variable", "date"), sep = "_") %>%
  separate(date, into = c("year", "month", "day"), sep = "-") %>%
  mutate(year = parse_number(year),
         month = parse_number(month),
         day = parse_number(day)) %>%
  arrange(year, month, day, variable, NA_L3NAME) %>% 
  pivot_wider(., names_from = variable, values_from = value)
eco_geom <- geom(ecoregion_shp) %>% as_tibble() %>% rename(id = geom) %>% left_join(eco_names) %>% dplyr::select(-c("id", "part", "hole"))
eco_lat_long <- eco_geom %>% group_by(value) %>% summarise(long = mean(x), lat = mean(y)) %>% rename(NA_L3NAME = value)
daily_climate_list <- daily_climate_eco %>% rename(yr = year, mon = month, temp = tmmx, rh = rmin, prec = pr, ws = vs) %>%
  left_join(eco_lat_long) %>% split(., .$NA_L3NAME)
daily_fwi_eco <- lapply(daily_climate_list, function(x) fwi(x, uppercase = FALSE, batch = TRUE, lat.adjust = TRUE))

destfile <- "./full-model/data/processed/ecoregion_summaries.csv"
write_csv(ecoregion_summaries, destfile)

print(paste('Ecoregion climate summaries written to', destfile))

# Process extracted values into a usable data frame -----------------------
ecoregion_summaries <- extractions %>%
  bind_cols %>%
  gather(variable, value, -NA_L3NAME) %>%
  filter(!grepl(pattern = 'NA_L3NAME', x = variable)) %>%
  separate(variable,
           into = c("interval", "variable", "timestep"),
           sep = "_") %>%
  separate(timestep, into = c("year", "month"), sep = "\\.") %>%
  dplyr::select(-interval) %>%
  mutate(year = parse_number(year),
         month = parse_number(month)) %>%
  arrange(year, month, variable, NA_L3NAME)

destfile <- "data/processed/ecoregion_summaries.csv"
write_csv(ecoregion_summaries, destfile)

print(paste('Ecoregion climate summaries written to', destfile))

