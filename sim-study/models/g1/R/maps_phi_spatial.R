library(tidyverse)
library(lubridate)
#library(rstan)
library(splines)
library(spdep)
library(extRemes)
library(sf)
library(RColorBrewer)
library(patchwork)
library(classInt)
library(assertthat)
library(spatialreg)


nb <- read_rds('~/Desktop/csu/research/josephs_paper/data/processed/nb.rds')
ecoregions <- read_rds(file = "ecoregions.RDS")

ecoregions_geom <- ecoregions %>% filter(!NA_L2NAME == "UPPER GILA MOUNTAINS (?)")

ecoregions_geom %>%
  ggplot() +
  geom_sf(size = .1, fill = 'white') +
  geom_sf(data = ecoregions_geom,
          aes(fill=NA_L1CODE), alpha = 0.6, lwd = 0, inherit.aes = FALSE) +
  geom_sf_label(aes(label = NA_L2CODE)) +
  theme_minimal() +
  theme(panel.grid.major = element_line(colour = "lightgrey"))

# # split into list of dataframes, one for each region at level 2 ------
# sep_regions_l2 <- split(full_burns_complete_nogeo, full_burns_complete_nogeo$NA_L2NAME)
# 
# # function to build GPD model ####
# evd_fit <- function(df, t) { 
#   fevd(x = BurnBndAc, data = df, threshold = t, type = "GP",
#                                 time.units = "months", na.action = na.omit)
# }
# # function that will allow mapply to still work even if there is an error
# trycatch_fit <- function(df, t) {
#   return(tryCatch(evd_fit(df, t), error=function(e) NULL))
# }
# 
# # get empirical quantile (95) to find threshold -------
# thres_95 <- function(df) as.numeric(quantile(df$BurnBndAc, probs = 0.95)[1])
# thres_l2_95 <- lapply(sep_regions_l2, thres_95)
# 
# # find regions with less than 20 observations
# exc_20 <- function(df) ifelse(nrow(df)*.05 < 20, 0, 1)
# exc_20_l2_95 <- lapply(sep_regions_l2, exc_20)
# 
# ### fit model and obtain estimators - level 2 ------
# fit_l2_95 <- mapply(df = sep_regions_l2, t = thres_l2_95, trycatch_fit)
# 
# get_estimators <- function(model) return(t(as.matrix(model$results$par)))
# catch_estimators <- function(model) return(tryCatch(get_estimators(model), error=function(e) NULL))
# 
# estimators_l2_95 <- lapply(fit_l2_95, catch_estimators)
# params_l2_95 <- do.call(rbind.data.frame, estimators_l2_95) # this does NOT keep the regions that did not fit a model
# params_l2_95$NA_L2NAME <- rownames(params_l2_95)
# exc_20_df <- as.data.frame(unlist(exc_20_l2_95))
# colnames(exc_20_df) <- "include"
# exc_20_df$NA_L2NAME <- rownames(exc_20_df)
# 
# # add column to determine which regions to exclude from map
# params_l2_95 <- params_l2_95 %>% left_join(exc_20_df)
# # add column of threshold values
# thres_l2_95_col <- do.call(cbind.data.frame, thres_l2_95) %>% t() %>% as.data.frame()
# thres_l2_95_col$NA_L2NAME <- rownames(thres_l2_95_col)
# class(thres_l2_95)
# class(estimators_l2_95)
# 
# t(thres_l2_95_col)
# 
# ## at 97 quantile ####
# thres_97 <- function(df) as.numeric(quantile(df$Acres, probs = 0.97)[1])
# thres_l2_97 <- lapply(sep_regions_l2, thres_97)
# fit_l2_97 <- mapply(df = sep_regions_l2, t = thres_l2_97, trycatch_fit)
# estimators_l2_97 <- lapply(fit_l2_97, catch_estimators)
# params_l2_97 <- do.call(rbind.data.frame, estimators_l2_97) # this does NOT keep the regions that did not fit a model
# params_l2_97$NA_L2NAME <- rownames(params_l2_97)
# 
# # exclude anything with less than 20 obs
# exc_20_97 <- function(df) ifelse(nrow(df)*.03 < 20, 0, 1)
# exc_20_l2_97 <- lapply(sep_regions_l2, exc_20_97)
# 
# exc_20_df_97 <- as.data.frame(unlist(exc_20_l2_97))
# colnames(exc_20_df_97) <- "include"
# exc_20_df_97$NA_L2NAME <- rownames(exc_20_df_97)
# 
# params_l2_97 <- params_l2_97 %>% left_join(exc_20_df_97)
# 
# ## level 3 ####
# sep_regions_l3 <- split(full_burns_complete_nogeo, full_burns_complete_nogeo$NA_L3NAME)
# thres_l3_95 <- lapply(sep_regions_l3, thres_95)
# fit_l3_95 <- mapply(df = sep_regions_l3, t = thres_l3_95, trycatch_fit)
# estimators_l3_95 <- lapply(fit_l3_95, catch_estimators)
# params_l3_95 <- do.call(rbind.data.frame, estimators_l3_95) # this does NOT keep the regions that did not fit a model
# params_l3_95$NA_L3NAME <- rownames(params_l3_95)
# # exclude anything with < 20 observations
# exc_20_l3_95 <- lapply(sep_regions_l3, exc_20) %>% unlist() %>% as.data.frame()
# colnames(exc_20_l3_95) <- "include"
# exc_20_l3_95$NA_L3NAME <- rownames(exc_20_l3_95)
# 
# params_l3_95 <- params_l3_95 %>% left_join(exc_20_l3_95)
# 
# ## at 97 quantile, level 3 ####
# thres_l3_97 <- lapply(sep_regions_l3, thres_97)
# fit_l3_97 <- mapply(df = sep_regions_l3, t = thres_l3_97, trycatch_fit)
# estimators_l3_97 <- lapply(fit_l3_97, catch_estimators)
# params_l3_97 <- do.call(rbind.data.frame, estimators_l3_97) # this does NOT keep the regions that did not fit a model
# params_l3_97$NA_L3NAME <- rownames(params_l3_97)
# 
# # exclude anything with < 20 observations
# exc_20_l3_97 <- lapply(sep_regions_l3, exc_20_97) %>% unlist() %>% as.data.frame()
# colnames(exc_20_l3_97) <- "include"
# exc_20_l3_97$NA_L3NAME <- rownames(exc_20_l3_97)
# 
# params_l3_97 <- params_l3_97 %>% left_join(exc_20_l3_97)
# 
# ## level 1 #### -----
# sep_regions_l1 <- split(full_burns_complete_nogeo, full_burns_complete_nogeo$NA_L1NAME)
# thres_l1_95 <- lapply(sep_regions_l1, thres_95)
# fit_l1_95 <- mapply(df = sep_regions_l1, t = thres_l1_95, trycatch_fit, SIMPLIFY = FALSE)
# estimators_l1_95 <- lapply(fit_l1_95, catch_estimators)
# params_l1_95 <- do.call(rbind.data.frame, estimators_l1_95) # this does NOT keep the regions that did not fit a model
# params_l1_95$NA_L1NAME <- rownames(params_l1_95)
# 
# # exclude anything with < 20 observations
# exc_20_l1_95 <- lapply(sep_regions_l1, exc_20) %>% unlist() %>% as.data.frame()
# colnames(exc_20_l1_95) <- "include"
# exc_20_l1_95$NA_L1NAME <- rownames(exc_20_l1_95)
# 
# params_l1_95 <- params_l1_95 %>% left_join(exc_20_l1_95)
# 
# ## at 97 quantile, level 1 ####
# thres_l1_97 <- lapply(sep_regions_l1, thres_97)
# fit_l1_97 <- mapply(df = sep_regions_l1, t = thres_l1_97, trycatch_fit)
# estimators_l1_97 <- lapply(fit_l1_97, catch_estimators)
# params_l1_97 <- do.call(rbind.data.frame, estimators_l1_97) # this does NOT keep the regions that did not fit a model
# params_l1_97$NA_L1NAME <- rownames(params_l1_97)
# 
# # exclude anything with < 20 observations
# exc_20_l1_97 <- lapply(sep_regions_l1, exc_20_97) %>% unlist() %>% as.data.frame()
# colnames(exc_20_l1_97) <- "include"
# exc_20_l1_97$NA_L1NAME <- rownames(exc_20_l1_97)
# 
# params_l1_97 <- params_l1_97 %>% left_join(exc_20_l1_97)
# 
# # histograms #### ---------
# hist(params_l2_95$shape)
# hist(params_l2_95$scale)
# 
# hist(params_l2_97$shape)
# hist(params_l2_97$scale)
# 
# hist(params_l3_95$shape)
# hist(params_l3_95$scale)
# 
# hist(params_l3_97$shape)
# hist(params_l3_97$scale)
# 
# plot(params_l3_97$scale)

## make maps ####
## LEVEL 3 ## --------
pre_phi_times <- phi_mat %>% as_tibble() %>% 
  rename_with(., ~ reg_cols) %>%
  mutate(timepoint = 1:all_of(t)) %>%
  pivot_longer(cols = c(1:all_of(r)), values_to = "phi_val", names_to = "region") %>%
  left_join(., mod_reg_key) %>% mutate(type = "truth")

joined_l3_phi_pre_times <- right_join(ecoregions_geom, pre_phi_times)
breaks <- classIntervals(c(min(joined_l3_phi_pre_times$phi_val) - .00001, joined_l3_phi_pre_times$phi_val), n = 5, style = "quantile")
joined_l3_phi_pre_times <- joined_l3_phi_pre_times %>% mutate(phi_cat = cut(phi_val, unique(breaks$brks)))

pre_phi_l3_onetime <- ecoregions_geom %>%
  ggplot() +
  geom_sf(size = .1, fill = 'white') +
  geom_sf(data = joined_l3_phi_pre_times[joined_l3_phi_pre_times$timepoint == 1, ],
          aes(fill=phi_cat), alpha = 0.6, lwd = 0, inherit.aes = FALSE) +
  theme_minimal() +
  theme(panel.grid.major = element_line(colour = "lightgrey"))
pre_phi_l3_onetime

post_phi <- rstan::extract(egpd_fit, pars = "phi")
median_phi <- apply(post_phi$phi, c(2,3), median)
post_phi_times <- median_phi %>% as_tibble() %>%
  rename_with(., ~ reg_cols) %>%
  mutate(timepoint = 1:all_of(t)) %>%
  pivot_longer(cols = c(1:all_of(r)), values_to = "phi_val", names_to = "region") %>%
  left_join(., mod_reg_key) %>% mutate(type = "sim")

joined_l3_phi_post_times <- right_join(ecoregions_geom, post_phi_times)
joined_l3_phi_post_times <- joined_l3_phi_post_times %>% mutate(phi_cat = cut(phi_val, unique(breaks$brks)))

# post_phi_l3_onetime <- ecoregions_geom %>%
#   ggplot() +
#   geom_sf(size = .1, fill = 'white') + 
#   geom_sf(data = joined_l3_phi_post_times[joined_l3_phi_post_times$timepoint == 1, ], aes(fill=phi_cat), alpha = 0.6, lwd = 0, inherit.aes = FALSE) + 
#   theme_minimal() + 
#   theme(panel.grid.major = element_line(colour = "lightgrey"))
# post_phi_l3_onetime

phi_full <- rbind(pre_phi_times, post_phi_times) %>% mutate(type = factor(type, levels = c("truth", "sim")))
joined_l3_phi_full <- right_join(ecoregions_geom, phi_full)
breaks <- classIntervals(c(min(joined_l3_phi_full$phi_val) - .00001, joined_l3_phi_full$phi_val), n=5, style = "quantile")
joined_l3_phi_full <- joined_l3_phi_full %>% mutate(phi_cat = cut(phi_val, unique(breaks$brks)))

full_phi_l3_onetime <- ecoregions_geom %>%
  ggplot() +
  geom_sf(size = .1, fill = 'white') + 
  geom_sf(data = joined_l3_phi_full[joined_l3_phi_full$timepoint == 99, ], 
          aes(fill=phi_cat), alpha = 0.6, lwd = 0, inherit.aes = FALSE) + 
  facet_grid(type ~ .) +
  theme_minimal() + 
  theme(panel.grid.major = element_line(colour = "lightgrey"))
full_phi_l3_onetime

