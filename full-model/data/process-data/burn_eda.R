library(tidyverse)
library(lubridate)
library(cmdstanr)
library(assertthat)
library(raster)
library(sp)
library(spdep)
library(spatialreg)
library(splines)
library(extRemes)

source("./full-model/data/process-data/merge-data.R")
ecoregions <- read_rds(file = "ecoregions.RDS")
ecoregions_geom <- ecoregions %>% filter(!NA_L2NAME == "UPPER GILA MOUNTAINS (?)") %>% 
  mutate(NA_L2CODE = as.factor(NA_L2CODE),
         NA_L1CODE = as.factor(NA_L1CODE),
         NA_L3CODE = as.factor(NA_L3CODE),
         NA_L1NAME = as.factor(str_to_title(NA_L1NAME)))

burn_df <- mtbs_er %>% as_tibble() %>% 
  dplyr::select(BurnBndAc, fire_yr, fire_mon, ym, 
                NA_L3CODE, NA_L3NAME, NA_L2CODE, NA_L2NAME, NA_L1CODE, NA_L1NAME, 
                Shape_Leng, Shape_Area) %>%
  arrange(NA_L3NAME, ym, .locale = "en") %>%
  filter(NA_L2NAME != "UPPER GILA MOUNTAINS (?)")

burn_df_l2 <- split(burn_df, burn_df$NA_L2NAME)

# function to build GPD model ####
evd_fit <- function(df, t) {
  fevd(x = df$BurnBndAc, threshold = t, type = "GP",
                                time.units = "months", na.action = na.omit)
}
# function that will allow mapply to still work even if there is an error
trycatch_fit <- function(df, t) {
  return(tryCatch(evd_fit(df, t), error=function(e) NULL))
}

# get empirical quantile (95) to find threshold -------
thres_95 <- function(df) as.numeric(quantile(df$BurnBndAc, probs = 0.95)[1])
thres_l2_95 <- lapply(burn_df_l2, thres_95)

# find regions with less than 20 observations
exc_20 <- function(df) ifelse(nrow(df)*.05 < 20, 0, 1)
exc_20_l2_95 <- lapply(burn_df_l2, exc_20)

### fit model and obtain estimators - level 2 ------
fit_l2_95 <- mapply(df = burn_df_l2, t = thres_l2_95, trycatch_fit)

get_estimators <- function(model) return(t(as.matrix(model$results$par)))
catch_estimators <- function(model) return(tryCatch(get_estimators(model), error=function(e) NULL))
estimators_l2_95 <- lapply(fit_l2_95, catch_estimators)
params_l2_95 <- do.call(rbind.data.frame, estimators_l2_95) # this does NOT keep the regions that did not fit a model
params_l2_95$NA_L2NAME <- rownames(params_l2_95)
exc_20_df <- as.data.frame(unlist(exc_20_l2_95)) %>% rename(include = "unlist(exc_20_l2_95)") %>%
  rownames_to_column(var = "NA_L2NAME")

# add column to determine which regions to exclude from map
params_l2_95 <- params_l2_95 %>% left_join(exc_20_df)
# add column of threshold values
thres_l2_95_col <- do.call(cbind.data.frame, thres_l2_95) %>% 
  t() %>% 
  as.data.frame() %>%
  rename(threshold = "V1") %>%
  rownames_to_column(var = "NA_L2NAME")
params_l2_95 <- params_l2_95 %>% mutate(shape_include = shape*include)
saveRDS(params_l2_95, file = "full-model/figures/paper/burn_eda.RDS")
