##############################################################################
##  This script produces tables of scores in for the burned areas           ##
##  and occurrences submodels, according to the approach outlined in the    ##
##  paper.                                                                  ##
##                                                                          ##
##  The tables are saved for use in figures/scores_tables.Rmd to create     ##
##  LaTeX tables.                                                           ##
##############################################################################

library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(tidyverse)
library(stringr)
library(posterior)

# function to read each score file and appropriately rename
read_scores <- function(file, model_name) {
  temp <- readRDS(file)
  assign(model_name, temp, parent.frame())
  rm(temp)
  gc()
}
# function to extract scores from stan csvs
extraction <- function(file_group, model) {
  holdout_loglik <- read_cmdstan_csv(file_group, variables = "holdout_loglik")$post_warmup_draws
  holdout_twcrps <- read_cmdstan_csv(file_group, variables = "holdout_twcrps")$post_warmup_draws
  temp <- list(holdout_loglik = holdout_loglik, holdout_twcrps = holdout_twcrps)
  assign(model, temp, parent.frame())
  gc()
}

# Size submodel comparisons and pipeline ----------------------
### Phase one: pick best parameter permutation within G1 carrier family of EGPD ------
score_files <- paste0("./scores_traceplots/extracted_scores/",
                      list.files("./scores_traceplots/extracted_scores/", 
                                 "g1"))
# score_files <- score_files[grepl("climate", score_files)]
g1_model_names <- basename(score_files) %>% 
  str_remove(., "_climate_\\d{2}\\w{3}2023_\\d{4}_scores.RDS") %>%
  str_remove(., "g1_")

for(i in seq_along(g1_model_names)) {
  read_scores(score_files[i], g1_model_names[i])
}

## aggregate scores across the four models
# not looking at training scores here, but this can easily be 
# modified to accommodate aggregating those as well
nfits <- length(g1_model_names)
holdout_loglik <- vector("list", nfits)
holdout_twcrps <- vector("list", nfits)
for(i in seq_along(g1_model_names)) {
  holdout_loglik[[i]] <- get(g1_model_names[i])[["holdout_loglik"]] %>%
    as_draws_df() %>%
    select(-c(".iteration", ".chain")) %>%
    pivot_longer(cols = !".draw") %>%
    rename(draw = ".draw") %>%
    group_by(draw) %>%
    summarize(loglik = sum(value)) %>%
    mutate(model = g1_model_names[i])
  holdout_twcrps[[i]] <- get(g1_model_names[i])[["holdout_twcrps"]] %>%
    as_draws_df() %>%
    select(-c(".iteration", ".chain")) %>%
    pivot_longer(cols = !".draw") %>%
    rename(draw = ".draw") %>%
    group_by(draw) %>%
    summarize(twcrps = mean(value)) %>%
    mutate(model = g1_model_names[i])
}

holdout_loglik <- bind_rows(holdout_loglik)
holdout_twcrps <- bind_rows(holdout_twcrps)

## determine best model based on log-likelihood scores 
ll_ranked <- holdout_loglik %>%
  group_by(model) %>% 
  summarize(med_train_ll = median(loglik)) %>% 
  arrange(-med_train_ll)
top_mod_ll <- as.character(ll_ranked$model[1])
ll_comp <- holdout_loglik %>%
  select(c(draw, model, loglik)) %>% 
  pivot_wider(names_from = model, values_from = loglik, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_ll))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(ll_med_diff = median(value[is.finite(value)]), ll_sd_diff = sd(value[is.finite(value)])) %>% 
  arrange(-ll_med_diff)
ll_comp

## determine best model based on twCRPS
twcrps_ranked <- holdout_twcrps %>% 
  group_by(model) %>% 
  summarize(mean_twcrps = mean(twcrps, na.rm = TRUE)) %>% 
  arrange(mean_twcrps)
top_mod_twcrps <- as.character(twcrps_ranked$model[1])
twcrps_comp <- holdout_twcrps %>% mutate(twcrps = twcrps * 10^2) %>%
  select(c(draw, model, twcrps)) %>% 
  pivot_wider(names_from = model, values_from = twcrps, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_twcrps))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(twcrps_mean_diff = mean(value[is.finite(value)]), twcrps_sd_diff = sd(value[is.finite(value)])) %>% 
  arrange(twcrps_mean_diff)
twcrps_comp

g1_params_comp <- ll_comp %>% left_join(twcrps_comp)
# save for use in markdown to create tables
# saveRDS(g1_params_comp, "./figures/g1_parameter_comparison.RDS")

### Phase two: pick best dataset on best G1 permutation  -----
## load in G1 model (climate dataset) scores extracted at time of model run on HPC 
score_files <- paste0("./scores_traceplots/extracted_scores/",
                      list.files("./scores_traceplots/extracted_scores/", 
                                 "g1_sigma-ri"))
g1_model_names <- basename(score_files) %>% 
  str_remove(., "_\\d{2}\\w{3}2023_\\d{4}_scores.RDS") %>%
  str_remove(., "g1_sigma-ri_xi-ri_")

for(i in seq_along(g1_model_names)) {
  read_scores(score_files[i], g1_model_names[i])
}

## aggregate scores across the four datasets
nfits <- length(g1_model_names)
holdout_loglik <- vector("list", nfits)
holdout_twcrps <- vector("list", nfits)
for(i in seq_along(g1_model_names)) {
  holdout_loglik[[i]] <- get(g1_model_names[i])[["holdout_loglik"]] %>%
    as_draws_df() %>%
    select(-c(".iteration", ".chain")) %>%
    pivot_longer(cols = !".draw") %>%
    rename(draw = ".draw") %>%
    group_by(draw) %>%
    summarize(loglik = sum(value)) %>%
    mutate(model = g1_model_names[i])
  holdout_twcrps[[i]] <- get(g1_model_names[i])[["holdout_twcrps"]] %>%
    as_draws_df() %>%
    select(-c(".iteration", ".chain")) %>%
    pivot_longer(cols = !".draw") %>%
    rename(draw = ".draw") %>%
    group_by(draw) %>%
    summarize(twcrps = mean(value)) %>%
    mutate(model = g1_model_names[i])
}

holdout_loglik <- bind_rows(holdout_loglik)
holdout_twcrps <- bind_rows(holdout_twcrps)

## determine best model based on log-likelihood scores
ll_ranked <- holdout_loglik %>%
  group_by(model) %>% 
  summarize(med_train_ll = median(loglik)) %>% 
  arrange(-med_train_ll)
top_mod_ll <- as.character(ll_ranked$model[1])
ll_comp <- holdout_loglik %>%
  select(c(draw, model, loglik)) %>% 
  pivot_wider(names_from = model, values_from = loglik, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_ll))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(ll_med_diff = median(value[is.finite(value)]), ll_sd_diff = sd(value[is.finite(value)])) %>% 
  arrange(-ll_med_diff)
ll_comp

## determine best model based on twCRPS
twcrps_ranked <- holdout_twcrps %>% 
  group_by(model) %>% 
  summarize(mean_twcrps = mean(twcrps, na.rm = TRUE)) %>% 
  arrange(mean_twcrps)
top_mod_twcrps <- as.character(twcrps_ranked$model[1])
twcrps_comp <- holdout_twcrps %>% mutate(twcrps = twcrps * 10^2) %>%
  select(c(draw, model, twcrps)) %>% 
  pivot_wider(names_from = model, values_from = twcrps, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_twcrps))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(twcrps_mean_diff = mean(value[is.finite(value)]), twcrps_sd_diff = sd(value[is.finite(value)])) %>% 
  arrange(twcrps_mean_diff)
twcrps_comp

g1_dataset_comp <- ll_comp %>% left_join(twcrps_comp)
# save for use in markdown to create tables
# saveRDS(g1_dataset_comp, "./figures/g1_dataset_comparison.RDS")

### Phase three: pick best burned sizes model run on best dataset -----
# NB: while 'ERC-FWI' is the best dataset, the G3 and G4 carrier families did not 
# mix well on this dataset, so scores here are shown for 'ERC' dataset
# G3 and G4 were run on a different machine since they take longer than 
# seven days to run, so those scores were not extracted like G1, G2, and lognorm
# models. These will need to be read directly from the csvs
score_files <- paste0("./scores_traceplots/extracted_scores/",
                      list.files("./scores_traceplots/extracted_scores/", 
                                 "erc_fwi"))
score_files <- score_files[!grepl("zinb", score_files)]
sizes_model_names <- basename(score_files) %>% 
  str_remove(., "_erc_fwi_\\d{2}\\w{3}2023_\\d{4}_scores.RDS")

for(i in seq_along(sizes_model_names)) {
  read_scores(score_files[i], sizes_model_names[i])
}

g3_files <- paste0("./models/sizes/g3/csv_fits/",
                   list.files("./models/sizes/g3/csv_fits/", 
                              "erc_fwi"))
g3_files <- g3_files[!grepl("nu", g3_files)]
g3_name <- basename(g3_files[1]) %>% str_remove(., "_erc_fwi_\\d{2}\\w{3}2023_\\d{4}_1.csv")
extraction(g3_files, g3_name)

g4_files <- paste0("./models/sizes/g4/csv_fits/",
                   list.files("./models/sizes/g4/csv_fits/", 
                              "erc_fwi"))
g4_name <- basename(g4_files[1]) %>% str_remove(., "_erc_fwi_\\d{2}\\w{3}2023_\\d{4}_1.csv")
extraction(g4_files, g4_name)

size_model_names <- c(sizes_model_names, g3_name, g4_name)

## aggregate scores across the four datasets
nfits <- length(size_model_names)
holdout_loglik <- vector("list", nfits)
holdout_twcrps <- vector("list", nfits)
for(i in seq_along(size_model_names)) {
  holdout_loglik[[i]] <- get(size_model_names[i])[["holdout_loglik"]] %>%
    as_draws_df() %>%
    select(-c(".iteration", ".chain")) %>%
    pivot_longer(cols = !".draw") %>%
    rename(draw = ".draw") %>%
    group_by(draw) %>%
    summarize(loglik = sum(value)) %>%
    mutate(model = size_model_names[i])
  holdout_twcrps[[i]] <- get(size_model_names[i])[["holdout_twcrps"]] %>%
    as_draws_df() %>%
    select(-c(".iteration", ".chain")) %>%
    pivot_longer(cols = !".draw") %>%
    rename(draw = ".draw") %>%
    group_by(draw) %>%
    summarize(twcrps = mean(value)) %>%
    mutate(model = size_model_names[i])
}

holdout_loglik <- bind_rows(holdout_loglik)
holdout_twcrps <- bind_rows(holdout_twcrps)

## determine best model based on log-likelihood scores
ll_ranked <- holdout_loglik %>%
  group_by(model) %>% 
  summarize(med_train_ll = median(loglik)) %>% 
  arrange(-med_train_ll)
top_mod_ll <- as.character(ll_ranked$model[1])
ll_comp <- holdout_loglik %>%
  select(c(draw, model, loglik)) %>% 
  pivot_wider(names_from = model, values_from = loglik, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_ll))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(ll_med_diff = median(value[is.finite(value)]), ll_sd_diff = sd(value[is.finite(value)])) %>% 
  arrange(-ll_med_diff)
ll_comp

## determine best model based on twCRPS
twcrps_ranked <- holdout_twcrps %>% 
  group_by(model) %>% 
  summarize(mean_twcrps = mean(twcrps, na.rm = TRUE)) %>% 
  arrange(mean_twcrps)
top_mod_twcrps <- as.character(twcrps_ranked$model[1])
twcrps_comp <- holdout_twcrps %>% mutate(twcrps = twcrps * 10^2) %>%
  select(c(draw, model, twcrps)) %>% 
  pivot_wider(names_from = model, values_from = twcrps, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_twcrps))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(twcrps_mean_diff = mean(value[is.finite(value)]), twcrps_sd_diff = sd(value[is.finite(value)])) %>% 
  arrange(twcrps_mean_diff)
twcrps_comp

sizes_models_comp <- ll_comp %>% left_join(twcrps_comp)
# save for use in markdown to create tables
# saveRDS(sizes_models_comp, "./figures/sizes_model_comparison.RDS")

# Occurrences submodel comparisons and pipeline  ----------------------
### Phase one: pick best count model based on climate covariates dataset ### -----
score_files <- paste0("./scores_traceplots/extracted_scores/",
                      list.files("./scores_traceplots/extracted_scores/", 
                                 "zi"))
count_model_names <- basename(score_files) %>% 
  str_remove(., "_climate_\\d{2}\\w{3}2023_\\d{4}_scores.RDS")

for(i in seq_along(count_model_names)) {
  read_scores(score_files[i], count_model_names[i])
}

## aggregate scores across the six models
nfits <- length(count_model_names)
holdout_loglik <- vector("list", nfits)
for(i in seq_along(count_model_names)) {
  holdout_loglik[[i]] <- get(count_model_names[i])[["holdout_loglik"]] %>%
    as_draws_df() %>%
    select(-c(".iteration", ".chain")) %>%
    pivot_longer(cols = !".draw") %>%
    rename(draw = ".draw") %>%
    group_by(draw) %>%
    summarize(loglik = sum(value)) %>%
    mutate(model = count_model_names[i])
}
holdout_loglik <- bind_rows(holdout_loglik)

## determine best model based on log-likelihood scores 
ll_ranked <- holdout_loglik %>%
  group_by(model) %>% 
  summarize(med_train_ll = median(loglik)) %>% 
  arrange(-med_train_ll)
top_mod_ll <- as.character(ll_ranked$model[1])
ll_comp <- holdout_loglik %>%
  select(c(draw, model, loglik)) %>% 
  pivot_wider(names_from = model, values_from = loglik, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_ll))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(ll_med_diff = median(value[is.finite(value)]), ll_sd_diff = sd(value[is.finite(value)])) %>% 
  arrange(-ll_med_diff)
ll_comp

# saveRDS(ll_comp, "./figures/counts_model_comparison.RDS")

### Phase two: pick best dataset on best countmodel ----------------------------
score_files <- paste0("./scores_traceplots/extracted_scores/",
                      list.files("./scores_traceplots/extracted_scores/", 
                                 "zinb_er_pi-ri"))
count_model_names <- basename(score_files) %>% 
  str_remove(., "_\\d{2}\\w{3}2023_\\d{4}_scores.RDS") %>%
  str_remove(., "zinb_er_pi-ri_")

for(i in seq_along(count_model_names)) {
  read_scores(score_files[i], count_model_names[i])
}

## aggregate scores across the six models
nfits <- length(count_model_names)
holdout_loglik <- vector("list", nfits)
for(i in seq_along(count_model_names)) {
  holdout_loglik[[i]] <- get(count_model_names[i])[["holdout_loglik"]] %>%
    as_draws_df() %>%
    select(-c(".iteration", ".chain")) %>%
    pivot_longer(cols = !".draw") %>%
    rename(draw = ".draw") %>%
    group_by(draw) %>%
    summarize(loglik = sum(value)) %>%
    mutate(model = count_model_names[i])
}
holdout_loglik <- bind_rows(holdout_loglik)

## determine best model based on log-likelihood scores
ll_ranked <- holdout_loglik %>%
  group_by(model) %>% 
  summarize(med_train_ll = median(loglik)) %>% 
  arrange(-med_train_ll)
top_mod_ll <- as.character(ll_ranked$model[1])
ll_comp <- holdout_loglik %>%
  select(c(draw, model, loglik)) %>% 
  pivot_wider(names_from = model, values_from = loglik, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_ll))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(ll_med_diff = median(value[is.finite(value)]), ll_sd_diff = sd(value[is.finite(value)])) %>% 
  arrange(-ll_med_diff)
ll_comp

# saveRDS(ll_comp, "./figures/zinb_er_dataset_comparison.RDS")