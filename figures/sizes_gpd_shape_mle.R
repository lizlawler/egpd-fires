##############################################################################
##  This script produces a dataframe of MLE estimates of the GPD shape      ##
##  parameter for burned areas by L2 ecoregion.                             ##
##############################################################################

library(tidyverse)
library(lubridate)
library(extRemes)

source("./data/03_merge_data.R")
ecoregions <- ecoregion_shp %>% as_tibble() %>%
  mutate(NA_L2CODE = as.factor(NA_L2CODE),
         NA_L1CODE = as.factor(NA_L1CODE),
         NA_L3CODE = as.factor(NA_L3CODE),
         NA_L1NAME = as.factor(str_to_title(NA_L1NAME)))

size_df <- mtbs_er %>% as_tibble() %>% 
  dplyr::select(BurnBndAc, fire_yr, fire_mon, ym, 
                NA_L3CODE, NA_L3NAME, NA_L2CODE, NA_L2NAME, NA_L1CODE, NA_L1NAME, 
                Shape_Leng, Shape_Area) %>%
  dplyr::rename(year = fire_yr, month = fire_mon) %>%
  arrange(NA_L2NAME, ym, .locale = "en")

size_df_l2 <- split(size_df, size_df$NA_L2NAME)

# function to fit GPD
evd_fit <- function(df, t) {
  fevd(x = df$BurnBndAc, data = df, threshold = t, type = "GP", time.units = "months")
}
# function that will allow mapply to still work even if there is an error
trycatch_fit <- function(df, t) {
  return(tryCatch(evd_fit(df, t), error=function(e) NULL))
}

# get empirical quantile (95) to find threshold -------
thres_90 <- function(df) as.numeric(quantile(df$BurnBndAc, probs = 0.90)[1])
thres_l2_90 <- lapply(size_df_l2, thres_90)

# find regions with less than 20 observations
n_20_90 <- function(df) floor(nrow(df) * 0.10)
n_20_l2_90 <- lapply(size_df_l2, n_20_90)

### fit model and obtain estimators - level 2 ------
fit_l2_90 <- mapply(df = size_df_l2, t = thres_l2_90, trycatch_fit)

get_estimators <- function(model) return(t(as.matrix(model$results$par)))
trycatch_estimators <- function(model) return(tryCatch(get_estimators(model), error=function(e) NULL))
estimators_l2_90 <- lapply(fit_l2_90, trycatch_estimators)

params_l2_90 <- do.call(rbind.data.frame, estimators_l2_90) %>% rownames_to_column(var = "NA_L2NAME") %>% as_tibble() %>% dplyr::rename(scale90 = scale, shape90 = shape)

n_20_90_df <- as.data.frame(unlist(n_20_l2_90)) %>% rownames_to_column(var = "NA_L2NAME") %>% 
  as_tibble() %>% rename(n_over_90 = "unlist(n_20_l2_90)") %>% mutate(include90 = ifelse(n_over_90 < 20, NA, 1))

thres_90_df <- as.data.frame(unlist(thres_l2_90)) %>% rownames_to_column(var = "NA_L2NAME") %>% 
  as_tibble() %>% rename(thres_90 = "unlist(thres_l2_90)")

# add column to determine which regions to exclude from map
params_l2_90_df <- params_l2_90 %>% right_join(n_20_90_df) %>% mutate(shape_inc_90 = shape90*include90) %>% left_join(thres_90_df)
saveRDS(params_l2_90_df, file = "./figures/sizes_shape_mle.RDS")
