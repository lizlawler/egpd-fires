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
library(classInt)

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
  dplyr::rename(year = fire_yr, month = fire_mon) %>%
  arrange(NA_L2NAME, ym, .locale = "en") %>%
  filter(NA_L2NAME != "UPPER GILA MOUNTAINS (?)")

burn_df_l2 <- split(burn_df, burn_df$NA_L2NAME)

# function to build GPD model ####
evd_fit <- function(df, t) {
  fevd(x = df$BurnBndAc, data = df, threshold = t, type = "GP", time.units = "months")
}
# function that will allow mapply to still work even if there is an error
trycatch_fit <- function(df, t) {
  return(tryCatch(evd_fit(df, t), error=function(e) NULL))
}

# get empirical quantile (95) to find threshold -------
thres_95 <- function(df) as.numeric(quantile(df$BurnBndAc, probs = 0.95)[1])
thres_l2_95 <- lapply(burn_df_l2, thres_95)

thres_90 <- function(df) as.numeric(quantile(df$BurnBndAc, probs = 0.90)[1])
thres_l2_90 <- lapply(burn_df_l2, thres_90)

# find regions with less than 20 observations
n_20_95 <- function(df) floor(nrow(df) * 0.05)
n_20_l2_95 <- lapply(burn_df_l2, n_20_95)

n_20_90 <- function(df) floor(nrow(df) * 0.10)
n_20_l2_90 <- lapply(burn_df_l2, n_20_90)

### fit model and obtain estimators - level 2 ------
fit_l2_95 <- mapply(df = burn_df_l2, t = thres_l2_95, trycatch_fit)
fit_l2_90 <- mapply(df = burn_df_l2, t = thres_l2_90, trycatch_fit)

get_estimators <- function(model) return(t(as.matrix(model$results$par)))
trycatch_estimators <- function(model) return(tryCatch(get_estimators(model), error=function(e) NULL))
estimators_l2_95 <- lapply(fit_l2_95, trycatch_estimators)
estimators_l2_90 <- lapply(fit_l2_90, trycatch_estimators)

params_l2_95 <- do.call(rbind.data.frame, estimators_l2_95) %>% rownames_to_column(var = "NA_L2NAME") %>% as_tibble() %>% dplyr::rename(scale95 = scale, shape95 = shape)
params_l2_90 <- do.call(rbind.data.frame, estimators_l2_90) %>% rownames_to_column(var = "NA_L2NAME") %>% as_tibble() %>% dplyr::rename(scale90 = scale, shape90 = shape)

n_20_90_df <- as.data.frame(unlist(n_20_l2_90)) %>% rownames_to_column(var = "NA_L2NAME") %>% 
  as_tibble() %>% rename(n_over_90 = "unlist(n_20_l2_90)") %>% mutate(include90 = ifelse(n_over_90 < 20, NA, 1))
n_20_95_df <- as.data.frame(unlist(n_20_l2_95)) %>% rownames_to_column(var = "NA_L2NAME") %>% 
  as_tibble() %>% rename(n_over_95 = "unlist(n_20_l2_95)") %>% mutate(include95 = ifelse(n_over_95 < 20, NA, 1))

thres_90_df <- as.data.frame(unlist(thres_l2_90)) %>% rownames_to_column(var = "NA_L2NAME") %>% 
  as_tibble() %>% rename(thres_90 = "unlist(thres_l2_90)")
thres_95_df <- as.data.frame(unlist(thres_l2_95)) %>% rownames_to_column(var = "NA_L2NAME") %>% 
  as_tibble() %>% rename(thres_95 = "unlist(thres_l2_95)")

# add column to determine which regions to exclude from map
params_l2_95_df <- params_l2_95 %>% right_join(n_20_95_df) %>% mutate(shape_inc_95 = shape95*include95) %>% left_join(thres_95_df)
params_l2_90_df <- params_l2_90 %>% right_join(n_20_90_df) %>% mutate(shape_inc_90 = shape90*include90) %>% left_join(thres_90_df)

params_l2_df <- params_l2_90_df %>% left_join(params_l2_95_df) %>% dplyr::select(c(NA_L2NAME, shape90, n_over_90, shape_inc_90, shape95, n_over_95, shape_inc_95, thres_90, thres_95)) %>%
  relocate(c(n_over_90, n_over_95, thres_90, thres_95), .before = "shape90") %>% relocate(shape95, .before = "shape_inc_90")
saveRDS(params_l2_df, file = "full-model/figures/paper/burn_eda.RDS")

# breaks <- classIntervals(c(min(params_l2_df$shape_inc_95) - .00001, params_l2_df$shape_inc_95), style = 'fixed', 
#                          fixedBreaks = c(0, 0.2, 0.4, 0.6, 0.8, 2.0), intervalClosure = 'left')
breaks <- classIntervals(c(min(params_l2_df$shape_inc_95) - .00001, params_l2_df$shape_inc_95), style = 'pretty', intervalClosure = 'left')
eco_eda <- ecoregions_geom %>% left_join(params_l2_95_df) %>%
  mutate(xi_cat = cut(shape_inc_95, unique(breaks$brks)))

p <- ecoregions_geom %>% group_by(NA_L2NAME) %>% summarise() %>%
  ggplot() +
  geom_sf(size = .1, fill = 'transparent', inherit.aes = FALSE) +
  geom_sf(data = eco_eda %>% group_by(NA_L2NAME, shape_inc_95) %>% summarise(), 
          aes(fill=shape_inc_95), alpha = 0.6, lwd = 0, inherit.aes = FALSE) +
  theme_void()

p_color <- p + scale_fill_gradient2(low = "yellow", mid = "orange", high = "red", na.value = "transparent")
p_color <- p + scale_fill_gradient2(na.value = "transparent")
# + scale_fill_gradient2(na.value = "transparent")
ggsave("eda_map.png", bg ='white')
