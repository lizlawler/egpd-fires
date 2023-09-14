library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(tidyverse)
library(stringr)
library(posterior)

# load in scores extracted at time of model run on Alpine
score_files <- paste0("full-model/fire-sims/model_comparison/extracted_values/",
                      list.files("full-model/fire-sims/model_comparison/extracted_values/", 
                                 "scores.RDS"))
# process scores only from "normal" qos - run for 24 hrs or less
score_files <- score_files[grepl("normal", score_files)]
burn_files <- c(score_files[grepl("g1", score_files)], score_files[grepl("lognorm", score_files)])
count_files <- c(score_files[grepl("zinb", score_files)], score_files[grepl("zip", score_files)])

burn_model_names <- str_remove(basename(burn_files), "_\\d{2}\\w{3}2023_\\d{4}_normal_scores.RDS")
count_model_names <- str_remove(basename(count_files), "_\\d{2}\\w{3}2023_\\d{4}_normal_scores.RDS")

# function read each score file and appropriately rename
read_scores <- function(file, model_name) {
  temp <- readRDS(file)
  assign(model_name, temp, parent.frame())
  rm(temp)
  gc()
}

for(i in seq_along(burn_model_names)) {
  read_scores(burn_files[i], burn_model_names[i])
}

# pull scores from each model and aggreagate
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

ll_full_burns <- holdout_loglik_burns %>%
  full_join(train_loglik_burns)

train_ll_burns_ranked <- ll_full_burns %>% filter(train == TRUE) %>% 
  group_by(model) %>% 
  summarize(med_train_ll = median(loglik)) %>% arrange(-med_train_ll)
top_mod_burns_train <- as.character(train_ll_burns_ranked$full_name[1])

train_ll_burns_comp <- ll_full_burns %>% filter(train == TRUE) %>% 
  select(c(draw, full_name, loglik)) %>% 
  pivot_wider(names_from = full_name, values_from = loglik, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_burns_train))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(med_diff = median(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% arrange(-med_diff)

test_ll_burns_ranked <- ll_full_burns %>% filter(train == FALSE) %>% 
  group_by(model, rand_effect, link, full_name) %>% 
  summarize(med_train_ll = median(loglik)) %>% arrange(-med_train_ll)
top_mod_burns_test <- as.character(test_ll_burns_ranked$full_name[1])
test_ll_burns_comp <- ll_full_burns %>% filter(train == FALSE) %>% 
  select(c(draw, full_name, loglik)) %>% 
  pivot_wider(names_from = full_name, values_from = loglik, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_burns_test))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(med_diff = median(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% arrange(-med_diff)
test_ll_burns_comp


## twCRPS calculations ---------
twcrps_full <- holdout_twcrps %>%
  full_join(train_twcrps) %>% mutate(full_name = paste(model, rand_effect, link, sep = "_"))

train_twcrps_ranked <- twcrps_full %>% filter(train == TRUE) %>% 
  group_by(model, rand_effect, link, full_name) %>% 
  summarize(mean_twcrps = mean(twcrps, na.rm = TRUE)) %>% arrange(mean_twcrps)
top_mod_train_twcrps <- as.character(train_twcrps_ranked$full_name[1])

train_twcrps_comp <- twcrps_full %>% filter(train == TRUE) %>% 
  select(c(draw, full_name, twcrps)) %>% 
  pivot_wider(names_from = full_name, values_from = twcrps, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_train_twcrps))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(mean_diff = mean(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% arrange(mean_diff)

test_twcrps_ranked <- twcrps_full %>% filter(train == FALSE) %>% 
  group_by(model, rand_effect, link, full_name) %>% 
  summarize(mean_twcrps = mean(twcrps, na.rm = TRUE)) %>% arrange(mean_twcrps)
top_mod_test_twcrps <- as.character(test_twcrps_ranked$full_name[1])
test_twcrps_comp <- twcrps_full %>% filter(train == FALSE) %>% 
  select(c(draw, full_name, twcrps)) %>% 
  pivot_wider(names_from = full_name, values_from = twcrps, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_test_twcrps))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(mean_diff = mean(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% arrange(mean_diff)
test_twcrps_comp
