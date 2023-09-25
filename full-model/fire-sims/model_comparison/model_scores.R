library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(tidyverse)
library(stringr)
library(posterior)

# load in scores extracted at time of model run on Alpine
score_files <- paste0("full-model/fire-sims/model_comparison/extracted_values/",
                      list.files("full-model/fire-sims/model_comparison/extracted_values/", 
                                 "GQ.csv"))
# process scores only from "normal" qos - run for 24 hrs or less
# score_files <- score_files[grepl("normal", score_files)]
burn_files <- c(score_files[grepl("g1", score_files)], 
                score_files[grepl("lognorm", score_files)], 
                score_files[grepl("g2", score_files)])
# burn_files_climate <- burn_files[grepl("climate", burn_files)]
count_files <- c(score_files[grepl("zinb", score_files)], score_files[grepl("zip", score_files)])
count_files_climate <- count_files[grepl("climate", count_files)]

## remove burn models that either didn't mix or had >10% divergences ------
# g2: xi-ri, kappa-ri_xi-ri_nu, all-reg
# lognorm: all-reg
# g1: all-reg
burn_files_climate <- burn_files_climate[!grepl("all-reg", burn_files_climate)]
burn_files_climate <- burn_files_climate[!grepl("g2_xi-ri", burn_files_climate)]
burn_files_climate <- burn_files_climate[!grepl("g2_kappa-ri_xi-ri", burn_files_climate)]
burn_model_names <- str_remove(str_remove(str_remove(basename(burn_files_climate), "_\\d{2}\\w{3}2023_\\d{4}"), "_climate"), "_scores.RDS")
count_model_names <- str_remove(str_remove(basename(count_files_climate), "_\\d{2}\\w{3}2023_\\d{4}"), "_scores.RDS")

# function read each score file and appropriately rename
read_scores <- function(file, model_name) {
  temp <- readRDS(file)
  assign(model_name, temp, parent.frame())
  rm(temp)
  gc()
}

# pull scores from each BURN model and aggregate -----
for(i in seq_along(burn_model_names)) {
  read_scores(burn_files_climate[i], burn_model_names[i])
}

nfits <- length(burn_model_names)
holdout_loglik <- vector("list", nfits)
train_loglik <- vector("list", nfits)
train_twcrps <- vector("list", nfits)
holdout_twcrps <- vector("list", nfits)
for(i in seq_along(burn_model_names)) {
  holdout_loglik[[i]] <- get(burn_model_names[i])[["holdout_loglik"]] %>%
    as_draws_df() %>%
    select(-c(".iteration", ".chain")) %>% 
    pivot_longer(cols = !".draw") %>%
    rename(draw = ".draw") %>%
    group_by(draw) %>% 
    summarize(loglik = sum(value)) %>%
    mutate(model = burn_model_names[i],
           train = FALSE)
  holdout_twcrps[[i]] <- get(burn_model_names[i])[["holdout_twcrps"]] %>%
    as_draws_df() %>%
    select(-c(".iteration", ".chain")) %>% 
    pivot_longer(cols = !".draw") %>%
    rename(draw = ".draw") %>%
    group_by(draw) %>% 
    summarize(twcrps = mean(value)) %>%
    mutate(model = burn_model_names[i],
           train = FALSE)
  train_loglik[[i]] <- get(burn_model_names[i])[["train_loglik"]] %>%
    as_draws_df() %>%
    select(-c(".iteration", ".chain")) %>% 
    pivot_longer(cols = !".draw") %>%
    rename(draw = ".draw") %>%
    group_by(draw) %>% 
    summarize(loglik = sum(value)) %>%
    mutate(model = burn_model_names[i],
           train = TRUE)
  train_twcrps[[i]] <- get(burn_model_names[i])[["train_twcrps"]] %>%
    as_draws_df() %>%
    select(-c(".iteration", ".chain")) %>% 
    pivot_longer(cols = !".draw") %>%
    rename(draw = ".draw") %>%
    group_by(draw) %>% 
    summarize(twcrps = mean(value)) %>%
    mutate(model = burn_model_names[i],
           train = TRUE)
}

holdout_loglik_burns <- bind_rows(holdout_loglik)
train_loglik_burns <- bind_rows(train_loglik)
holdout_twcrps <- bind_rows(holdout_twcrps)
train_twcrps <- bind_rows(train_twcrps)

## log-likelihood aggregating and analysis (burns) -------
ll_full_burns <- holdout_loglik_burns %>%
  full_join(train_loglik_burns) %>% filter(grepl("long", model))

train_ll_burns_ranked <- ll_full_burns %>% filter(train == TRUE) %>% 
  group_by(model) %>% 
  summarize(med_train_ll = median(loglik[is.finite(loglik)])) %>% 
  arrange(-med_train_ll)
top_mod_burns_train <- as.character(train_ll_burns_ranked$model[1])

train_ll_burns_comp <- ll_full_burns %>% filter(train == TRUE) %>% 
  select(c(draw, model, loglik)) %>% 
  pivot_wider(names_from = model, values_from = loglik, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_burns_train))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(med_diff = median(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% 
  arrange(-med_diff)

test_ll_burns_ranked <- ll_full_burns %>% filter(train == FALSE) %>% 
  group_by(model) %>% 
  summarize(med_train_ll = median(loglik)) %>% 
  arrange(-med_train_ll)
top_mod_burns_test <- as.character(test_ll_burns_ranked$model[1])
test_ll_burns_comp <- ll_full_burns %>% filter(train == FALSE) %>% 
  select(c(draw, model, loglik)) %>% 
  pivot_wider(names_from = model, values_from = loglik, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_burns_test))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(med_diff = median(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% 
  arrange(-med_diff)
test_ll_burns_comp

## twCRPS aggregation and analysis ---------
twcrps_full <- holdout_twcrps %>%
  full_join(train_twcrps) %>% filter(grepl("long", model))

train_twcrps_ranked <- twcrps_full %>% filter(train == TRUE) %>% 
  group_by(model) %>% 
  summarize(mean_twcrps = mean(twcrps, na.rm = TRUE)) %>% 
  arrange(mean_twcrps)
top_mod_train_twcrps <- as.character(train_twcrps_ranked$model[1])

train_twcrps_comp <- twcrps_full %>% filter(train == TRUE) %>% 
  select(c(draw, model, twcrps)) %>% 
  pivot_wider(names_from = model, values_from = twcrps, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_train_twcrps))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(mean_diff = mean(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% 
  arrange(mean_diff)

test_twcrps_ranked <- twcrps_full %>% filter(train == FALSE) %>% 
  group_by(model) %>% 
  summarize(mean_twcrps = mean(twcrps, na.rm = TRUE)) %>% 
  arrange(mean_twcrps)
top_mod_test_twcrps <- as.character(test_twcrps_ranked$model[1])
test_twcrps_comp <- twcrps_full %>% filter(train == FALSE) %>% 
  select(c(draw, model, twcrps)) %>% 
  pivot_wider(names_from = model, values_from = twcrps, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_test_twcrps))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(mean_diff = mean(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% 
  arrange(mean_diff)
test_twcrps_comp

# pull scores from each COUNT model and aggregate -----
for(i in seq_along(count_model_names)) {
  read_scores(count_files_climate[i], count_model_names[i])
}

nfits <- length(count_model_names)
holdout_loglik <- vector("list", nfits)
train_loglik <- vector("list", nfits)
for(i in seq_along(count_model_names)) {
  holdout_loglik[[i]] <- get(count_model_names[i])[["holdout_loglik"]] %>%
    as_draws_df() %>%
    select(-c(".iteration", ".chain")) %>% 
    pivot_longer(cols = !".draw") %>%
    rename(draw = ".draw") %>%
    group_by(draw) %>% 
    summarize(loglik = sum(value)) %>%
    mutate(model = count_model_names[i],
           train = FALSE)
  train_loglik[[i]] <- get(count_model_names[i])[["train_loglik"]] %>%
    as_draws_df() %>%
    select(-c(".iteration", ".chain")) %>% 
    pivot_longer(cols = !".draw") %>%
    rename(draw = ".draw") %>%
    group_by(draw) %>% 
    summarize(loglik = sum(value)) %>%
    mutate(model = count_model_names[i],
           train = TRUE)
}

holdout_loglik_counts <- bind_rows(holdout_loglik)
train_loglik_counts <- bind_rows(train_loglik)

## log-likelihood aggregating and analysis (counts) -------
ll_full_counts <- holdout_loglik_counts %>%
  full_join(train_loglik_counts) %>% filter(grepl("long", model))

train_ll_counts_ranked <- ll_full_counts %>% filter(train == TRUE) %>% 
  group_by(model) %>% 
  summarize(med_train_ll = median(loglik)) %>% 
  arrange(-med_train_ll)
top_mod_counts_train <- as.character(train_ll_counts_ranked$model[1])

train_ll_counts_comp <- ll_full_counts %>% filter(train == TRUE) %>% 
  select(c(draw, model, loglik)) %>% 
  pivot_wider(names_from = model, values_from = loglik, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_counts_train))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(med_diff = median(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% 
  arrange(-med_diff)

test_ll_counts_ranked <- ll_full_counts %>% filter(train == FALSE) %>% 
  group_by(model) %>% 
  summarize(med_test_ll = median(loglik)) %>% 
  arrange(-med_test_ll)
top_mod_counts_test <- as.character(test_ll_counts_ranked$model[1])
test_ll_counts_comp <- ll_full_counts %>% filter(train == FALSE) %>% 
  select(c(draw, model, loglik)) %>% 
  pivot_wider(names_from = model, values_from = loglik, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_counts_test))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(med_diff = median(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% 
  arrange(-med_diff)
test_ll_counts_comp

## picking best dataset within the best model ------
best_burn_files <- burn_files[grepl("g1_sigma-ri_xi", burn_files)]
best_burn_files <- best_burn_files[grepl("long", best_burn_files)]
best_burn_names <- (str_remove(basename(best_burn_files), "_\\d{2}\\w{3}2023_\\d{4}_long_scores.RDS"))

best_count_files <- count_files[grepl("zinb_er_pi-ri", count_files)]
best_count_files <- best_count_files[grepl("long", best_count_files)]
best_count_names <- (str_remove(basename(best_count_files), "_\\d{2}\\w{3}2023_\\d{4}_long_scores.RDS"))

for(i in 2:length(best_burn_names)) {
  read_scores(best_burn_files[i], best_burn_names[i])
}

nfits <- length(best_burn_names)
holdout_loglik <- vector("list", nfits)
train_loglik <- vector("list", nfits)
train_twcrps <- vector("list", nfits)
holdout_twcrps <- vector("list", nfits)
for(i in seq_along(best_burn_names)) {
  holdout_loglik[[i]] <- get(best_burn_names[i])[["holdout_loglik"]] %>%
    as_draws_df() %>%
    select(-c(".iteration", ".chain")) %>% 
    pivot_longer(cols = !".draw") %>%
    rename(draw = ".draw") %>%
    group_by(draw) %>% 
    summarize(loglik = sum(value)) %>%
    mutate(model = best_burn_names[i],
           train = FALSE)
  holdout_twcrps[[i]] <- get(best_burn_names[i])[["holdout_twcrps"]] %>%
    as_draws_df() %>%
    select(-c(".iteration", ".chain")) %>% 
    pivot_longer(cols = !".draw") %>%
    rename(draw = ".draw") %>%
    group_by(draw) %>% 
    summarize(twcrps = mean(value)) %>%
    mutate(model = best_burn_names[i],
           train = FALSE)
  train_loglik[[i]] <- get(best_burn_names[i])[["train_loglik"]] %>%
    as_draws_df() %>%
    select(-c(".iteration", ".chain")) %>% 
    pivot_longer(cols = !".draw") %>%
    rename(draw = ".draw") %>%
    group_by(draw) %>% 
    summarize(loglik = sum(value)) %>%
    mutate(model = best_burn_names[i],
           train = TRUE)
  train_twcrps[[i]] <- get(best_burn_names[i])[["train_twcrps"]] %>%
    as_draws_df() %>%
    select(-c(".iteration", ".chain")) %>% 
    pivot_longer(cols = !".draw") %>%
    rename(draw = ".draw") %>%
    group_by(draw) %>% 
    summarize(twcrps = mean(value)) %>%
    mutate(model = best_burn_names[i],
           train = TRUE)
}

holdout_loglik_bestburnsdata <- bind_rows(holdout_loglik)
train_loglik_bestburnsdata <- bind_rows(train_loglik)
holdout_twcrps_bestdata <- bind_rows(holdout_twcrps)
train_twcrps_bestdata <- bind_rows(train_twcrps)

## log-likelihood aggregating and analysis (burns) -------
ll_full_burns_bestdata <- holdout_loglik_bestburnsdata %>%
  full_join(train_loglik_bestburnsdata)

train_ll_burns_ranked_data <- ll_full_burns_bestdata %>% filter(train == TRUE) %>% 
  group_by(model) %>% 
  summarize(med_train_ll = median(loglik[is.finite(loglik)])) %>% 
  arrange(-med_train_ll)
top_mod_burns_train_data <- as.character(train_ll_burns_ranked_data$model[1])

train_ll_burns_comp_data <- ll_full_burns_bestdata %>% filter(train == TRUE) %>% 
  select(c(draw, model, loglik)) %>% 
  pivot_wider(names_from = model, values_from = loglik, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_burns_train_data))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(med_diff = median(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% 
  arrange(-med_diff)

test_ll_burns_ranked_data <- ll_full_burns_bestdata %>% filter(train == FALSE) %>% 
  group_by(model) %>% 
  summarize(med_test_ll = median(loglik)) %>% 
  arrange(-med_test_ll)
top_mod_burns_test_data <- as.character(test_ll_burns_ranked_data$model[1])
test_ll_burns_comp_data <- ll_full_burns_bestdata %>% filter(train == FALSE) %>% 
  select(c(draw, model, loglik)) %>% 
  pivot_wider(names_from = model, values_from = loglik, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_burns_test_data))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(med_diff = median(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% 
  arrange(-med_diff)
test_ll_burns_comp_data

## twCRPS aggregation and analysis ---------
twcrps_full_data <- holdout_twcrps_bestdata %>%
  full_join(train_twcrps_bestdata)

train_twcrps_ranked_data <- twcrps_full_data %>% filter(train == TRUE) %>% 
  group_by(model) %>% 
  summarize(mean_twcrps = mean(twcrps, na.rm = TRUE)) %>% 
  arrange(mean_twcrps)
top_mod_train_twcrps_data <- as.character(train_twcrps_ranked_data$model[1])

train_twcrps_comp_data <- twcrps_full_data %>% filter(train == TRUE) %>% 
  select(c(draw, model, twcrps)) %>% 
  pivot_wider(names_from = model, values_from = twcrps, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_train_twcrps_data))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(mean_diff = mean(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% 
  arrange(mean_diff)

test_twcrps_ranked_data <- twcrps_full_data %>% filter(train == FALSE) %>% 
  group_by(model) %>% 
  summarize(mean_twcrps = mean(twcrps, na.rm = TRUE)) %>% 
  arrange(mean_twcrps)
top_mod_test_twcrps_data <- as.character(test_twcrps_ranked_data$model[1])
test_twcrps_comp_data <- twcrps_full_data %>% filter(train == FALSE) %>% 
  select(c(draw, model, twcrps)) %>% 
  pivot_wider(names_from = model, values_from = twcrps, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_test_twcrps_data))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(mean_diff = mean(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% 
  arrange(mean_diff)
test_twcrps_comp_data

for(i in 2:length(best_count_names)) {
  read_scores(best_count_files[i], best_count_names[i])
}

nfits <- length(best_count_names)
holdout_loglik <- vector("list", nfits)
train_loglik <- vector("list", nfits)
for(i in seq_along(best_count_names)) {
  holdout_loglik[[i]] <- get(best_count_names[i])[["holdout_loglik"]] %>%
    as_draws_df() %>%
    select(-c(".iteration", ".chain")) %>% 
    pivot_longer(cols = !".draw") %>%
    rename(draw = ".draw") %>%
    group_by(draw) %>% 
    summarize(loglik = sum(value)) %>%
    mutate(model = best_count_names[i],
           train = FALSE)
  train_loglik[[i]] <- get(best_count_names[i])[["train_loglik"]] %>%
    as_draws_df() %>%
    select(-c(".iteration", ".chain")) %>% 
    pivot_longer(cols = !".draw") %>%
    rename(draw = ".draw") %>%
    group_by(draw) %>% 
    summarize(loglik = sum(value)) %>%
    mutate(model = best_count_names[i],
           train = TRUE)
}

holdout_loglik_counts_data <- bind_rows(holdout_loglik)
train_loglik_counts_data <- bind_rows(train_loglik)

ll_full_counts_data <- holdout_loglik_counts_data %>%
  full_join(train_loglik_counts_data)

train_ll_counts_ranked_data <- ll_full_counts_data %>% filter(train == TRUE) %>% 
  group_by(model) %>% 
  summarize(med_train_ll = median(loglik)) %>% 
  arrange(-med_train_ll)
top_mod_counts_train_data <- as.character(train_ll_counts_ranked_data$model[1])

train_ll_counts_comp_data <- ll_full_counts_data %>% filter(train == TRUE) %>% 
  select(c(draw, model, loglik)) %>% 
  pivot_wider(names_from = model, values_from = loglik, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_counts_train_data))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(med_diff = median(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% 
  arrange(-med_diff)

test_ll_counts_ranked_data <- ll_full_counts_data %>% filter(train == FALSE) %>% 
  group_by(model) %>% 
  summarize(med_test_ll = median(loglik)) %>% 
  arrange(-med_test_ll)
top_mod_counts_test_data <- as.character(test_ll_counts_ranked_data$model[1])
test_ll_counts_comp_data <- ll_full_counts_data %>% filter(train == FALSE) %>% 
  select(c(draw, model, loglik)) %>% 
  pivot_wider(names_from = model, values_from = loglik, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_counts_test_data))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(med_diff = median(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% 
  arrange(-med_diff)
test_ll_counts_comp_data
