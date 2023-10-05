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
burn_files_climate <- burn_files[grepl("climate", burn_files)]
count_files <- paste0("full-model/fire-sims/model_comparison/extracted_values/",
                      list.files("full-model/fire-sims/model_comparison/extracted_values/", 
                                 "scores"))
count_files <- c(count_files[grepl("zinb", count_files)], count_files[grepl("zip", count_files)])
count_files_climate <- count_files[grepl("climate", count_files)]
  

## remove burn models that either didn't mix or had >10% divergences ------
# g2: xi-ri, kappa-ri_xi-ri_nu, all-reg
# lognorm: all-reg
# g1: all-reg
burn_files_climate <- burn_files_climate[!grepl("all-reg", burn_files_climate)]
burn_files_climate <- burn_files_climate[!grepl("g2_xi-ri", burn_files_climate)]
burn_files_climate <- burn_files_climate[!grepl("g2_kappa-ri_xi-ri", burn_files_climate)]
count_model_names <- str_remove(str_remove(basename(count_files_climate), "_\\d{2}\\w{3}2023_\\d{4}"), "_scores.RDS")

# function to read each score file and appropriately rename
read_scores <- function(file, model_name) {
  temp <- readRDS(file)
  assign(model_name, temp, parent.frame())
  rm(temp)
  gc()
}

# function to read each burn GQ csv and appropriately rename
extraction <- function(file_group, model) {
  train_loglik <- read_cmdstan_csv(file_group, variables = "train_loglik")
  holdout_loglik <- read_cmdstan_csv(file_group, variables = "holdout_loglik")
  train_twcrps <- read_cmdstan_csv(file_group, variables = "train_twcrps")
  holdout_twcrps <- read_cmdstan_csv(file_group, variables = "holdout_twcrps")
  temp <- list(train_loglik = train_loglik, holdout_loglik = holdout_loglik, 
               train_twcrps = train_twcrps, holdout_twcrps = holdout_twcrps)
  assign(model, temp, parent.frame())
  gc()
}

nfits <- length(burn_files_climate)/3
fit_groups <- vector(mode = "list", nfits)
for(i in 1:nfits) {
  fit_groups[[i]] <- burn_files_climate[(3*i-2):(3*i)]
}

burn_model_names <- lapply(fit_groups, function(x) {
  basename(x[[1]]) %>% 
    str_remove(., "_climate_\\d{2}\\w{3}2023_\\d{4}") %>%
    str_remove(., "_\\d{1}_GQ.csv")
  }) %>% 
  unlist()


# pull scores from each BURN model and aggregate -----
for(i in seq_along(count_model_names)) {
  read_scores(count_files_climate[i], count_model_names[i])
}

for(i in seq_along(burn_model_names)) {
  extraction(fit_groups[[i]], burn_model_names[i])
}

# read in G3 model that finished and is analogous to G1 "best" model
g3_files <- paste0("full-model/fire-sims/burns/g3/csv-fits/",
                      list.files("full-model/fire-sims/burns/g3/csv-fits/", 
                                 "erc_xi-ri_nu"))
g3_files <- g3_files[c(1,3)]
g3_name <- str_remove(basename(g3_files[1]), "_\\d{2}Sep2023_\\d{4}_\\d{1}.csv")
extraction(list(g3_files)[[1]], g3_name)


nfits <- length(burn_model_names) + 1
holdout_loglik <- vector("list", nfits)
train_loglik <- vector("list", nfits)
train_twcrps <- vector("list", nfits)
holdout_twcrps <- vector("list", nfits)
for(i in seq_along(burn_model_names)) {
  holdout_loglik[[i]] <- get(burn_model_names[i])[["holdout_loglik"]]$generated_quantities %>%
    as_draws_df() %>%
    select(-c(".iteration", ".chain")) %>% 
    pivot_longer(cols = !".draw") %>%
    rename(draw = ".draw") %>%
    group_by(draw) %>% 
    summarize(loglik = sum(value)) %>%
    mutate(model = burn_model_names[i],
           train = FALSE)
  holdout_twcrps[[i]] <- get(burn_model_names[i])[["holdout_twcrps"]]$generated_quantities %>%
    as_draws_df() %>%
    select(-c(".iteration", ".chain")) %>% 
    pivot_longer(cols = !".draw") %>%
    rename(draw = ".draw") %>%
    group_by(draw) %>% 
    summarize(twcrps = mean(value)) %>%
    mutate(model = burn_model_names[i],
           train = FALSE)
  train_loglik[[i]] <- get(burn_model_names[i])[["train_loglik"]]$generated_quantities %>%
    as_draws_df() %>%
    select(-c(".iteration", ".chain")) %>% 
    pivot_longer(cols = !".draw") %>%
    rename(draw = ".draw") %>%
    group_by(draw) %>% 
    summarize(loglik = sum(value)) %>%
    mutate(model = burn_model_names[i],
           train = TRUE)
  train_twcrps[[i]] <- get(burn_model_names[i])[["train_twcrps"]]$generated_quantities %>%
    as_draws_df() %>%
    select(-c(".iteration", ".chain")) %>% 
    pivot_longer(cols = !".draw") %>%
    rename(draw = ".draw") %>%
    group_by(draw) %>% 
    summarize(twcrps = mean(value)) %>%
    mutate(model = burn_model_names[i],
           train = TRUE)
}

# add in g3 scores
holdout_loglik[[13]] <- get(g3_name)[["holdout_loglik"]]$post_warmup_draws %>%
  as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>% 
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  group_by(draw) %>% 
  summarize(loglik = sum(value)) %>%
  mutate(model = g3_name,
         train = FALSE)
holdout_twcrps[[13]] <- get(g3_name)[["holdout_twcrps"]]$post_warmup_draws %>%
  as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>% 
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  group_by(draw) %>% 
  summarize(twcrps = mean(value)) %>%
  mutate(model = g3_name,
         train = FALSE)
train_loglik[[13]] <- get(g3_name)[["train_loglik"]]$post_warmup_draws %>%
  as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>% 
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  group_by(draw) %>% 
  summarize(loglik = sum(value)) %>%
  mutate(model = g3_name,
         train = TRUE)
train_twcrps[[13]] <- get(g3_name)[["train_twcrps"]]$post_warmup_draws %>%
  as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>% 
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  group_by(draw) %>% 
  summarize(twcrps = mean(value)) %>%
  mutate(model = g3_name,
         train = TRUE)

holdout_loglik_burns <- bind_rows(holdout_loglik)
train_loglik_burns <- bind_rows(train_loglik)
holdout_twcrps <- bind_rows(holdout_twcrps)
train_twcrps <- bind_rows(train_twcrps)

## log-likelihood aggregating and analysis (burns) -------
ll_full_burns <- holdout_loglik_burns %>%
  full_join(train_loglik_burns)

train_ll_burns_ranked <- ll_full_burns %>% filter(grepl("long", model) | grepl("g3", model)) %>%
  filter(train == TRUE) %>% 
  group_by(model) %>% 
  summarize(med_train_ll = median(loglik[is.finite(loglik)])) %>% 
  arrange(-med_train_ll)
top_mod_burns_train <- as.character(train_ll_burns_ranked$model[1])

train_ll_burns_comp <- ll_full_burns %>% filter(grepl("long", model) | grepl("g3", model)) %>%
  filter(train == TRUE) %>% 
  select(c(draw, model, loglik)) %>% 
  pivot_wider(names_from = model, values_from = loglik, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_burns_train))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(med_diff = median(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% 
  arrange(-med_diff)

test_ll_burns_ranked <- ll_full_burns %>% filter(grepl("long", model) | grepl("g3", model)) %>%
  filter(train == FALSE) %>% 
  group_by(model) %>% 
  summarize(med_train_ll = median(loglik)) %>% 
  arrange(-med_train_ll)
top_mod_burns_test <- as.character(test_ll_burns_ranked$model[1])
test_ll_burns_comp <- ll_full_burns %>% filter(grepl("long", model) | grepl("g3", model)) %>%
  filter(train == FALSE) %>% 
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
  full_join(train_twcrps) %>% filter(grepl("long", model) | grepl("g3", model))

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

# pull scores from G1 "best" model run on different datasets
g1_files <- paste0("full-model/fire-sims/model_comparison/extracted_values/",
                   list.files("full-model/fire-sims/model_comparison/extracted_values/", 
                              "g1_sigma-ri"))
g1_files <- g1_files[!grepl("RDS", g1_files)]

nfits <- length(g1_files)/3
fit_groups <- vector(mode = "list", nfits)
for(i in 1:nfits) {
  fit_groups[[i]] <- g1_files[(3*i-2):(3*i)]
}

g1_model_names <- lapply(fit_groups, function(x) {
  basename(x[[1]]) %>% 
    str_remove(., "_\\d{2}\\w{3}2023_\\d{4}") %>%
    str_remove(., "_\\d{1}_GQ.csv") %>%
    str_remove(., "sigma-ri_xi-ri_")
}) %>% 
  unlist()

extraction <- function(file_group, model) {
  temp <- read_cmdstan_csv(file_group)$generated_quantities
  assign(model, temp, parent.frame())
  rm(temp)
  gc()
}

for(i in seq_along(g1_model_names)) {
  tryCatch(extraction(fit_groups[[i]], g1_model_names[i]), error=function(e) NULL)
}

g1_g3_mod_names <- c(g1_model_names[1:4], g3_name)

nfits <- length(g1_g3_mod_names)
holdout_loglik <- vector("list", nfits)
train_loglik <- vector("list", nfits)
train_twcrps <- vector("list", nfits)
holdout_twcrps <- vector("list", nfits)
for(i in seq_along(g1_model_names[1:4])) {
  draws_df <- get(g1_model_names[i]) %>% 
    as_draws_df() %>%
    select(-c(".iteration", ".chain")) %>% 
    pivot_longer(cols = !".draw") %>%
    rename(draw = ".draw")
  holdout_loglik[[i]] <- draws_df %>% filter(grepl("holdout_loglik", name)) %>%
    group_by(draw) %>% 
    summarize(loglik = sum(value)) %>%
    mutate(model = g1_model_names[i],
           train = FALSE)
  holdout_twcrps[[i]] <- draws_df %>% filter(grepl("holdout_twcrps", name)) %>%
    group_by(draw) %>% 
    summarize(twcrps = mean(value)) %>%
    mutate(model = g1_model_names[i],
           train = FALSE)
  train_loglik[[i]] <- draws_df %>% filter(grepl("train_loglik", name)) %>%
    group_by(draw) %>% 
    summarize(loglik = sum(value)) %>%
    mutate(model = g1_model_names[i],
           train = TRUE)
  train_twcrps[[i]] <- draws_df %>% filter(grepl("train_twcrps", name)) %>%
    group_by(draw) %>% 
    summarize(twcrps = mean(value)) %>%
    mutate(model = g1_model_names[i],
           train = TRUE)
}

# add in g3 scores
holdout_loglik[[5]] <- get(g3_name)[["holdout_loglik"]]$post_warmup_draws %>%
  as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>% 
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  group_by(draw) %>% 
  summarize(loglik = sum(value)) %>%
  mutate(model = g3_name,
         train = FALSE)
holdout_twcrps[[5]] <- get(g3_name)[["holdout_twcrps"]]$post_warmup_draws %>%
  as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>% 
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  group_by(draw) %>% 
  summarize(twcrps = mean(value)) %>%
  mutate(model = g3_name,
         train = FALSE)
train_loglik[[5]] <- get(g3_name)[["train_loglik"]]$post_warmup_draws %>%
  as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>% 
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  group_by(draw) %>% 
  summarize(loglik = sum(value)) %>%
  mutate(model = g3_name,
         train = TRUE)
train_twcrps[[5]] <- get(g3_name)[["train_twcrps"]]$post_warmup_draws %>%
  as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>% 
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  group_by(draw) %>% 
  summarize(twcrps = mean(value)) %>%
  mutate(model = g3_name,
         train = TRUE)

holdout_loglik_burns_data <- bind_rows(holdout_loglik)
train_loglik_burns_data <- bind_rows(train_loglik)
holdout_twcrps_data <- bind_rows(holdout_twcrps)
train_twcrps_data <- bind_rows(train_twcrps)

## log-likelihood aggregating and analysis (burns) -------
ll_full_burns_data <- holdout_loglik_burns_data %>%
  full_join(train_loglik_burns_data)

train_ll_burns_data_ranked <- ll_full_burns_data %>%
  filter(train == TRUE) %>% 
  group_by(model) %>% 
  summarize(med_train_ll = median(loglik[is.finite(loglik)])) %>% 
  arrange(-med_train_ll)
top_mod_burns_data_train <- as.character(train_ll_burns_data_ranked$model[1])

train_ll_burns_data_comp <- ll_full_burns_data %>%
  filter(train == TRUE) %>% 
  select(c(draw, model, loglik)) %>% 
  pivot_wider(names_from = model, values_from = loglik, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_burns_data_train))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(med_diff = median(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% 
  arrange(-med_diff)

test_ll_burns_data_ranked <- ll_full_burns_data %>%
  filter(train == FALSE) %>% 
  group_by(model) %>% 
  summarize(med_train_ll = median(loglik)) %>% 
  arrange(-med_train_ll)
top_mod_burns_data_test <- as.character(test_ll_burns_data_ranked$model[1])
test_ll_burns_data_comp <- ll_full_burns_data %>%
  filter(train == FALSE) %>% 
  select(c(draw, model, loglik)) %>% 
  pivot_wider(names_from = model, values_from = loglik, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_burns_data_test))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(med_diff = median(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% 
  arrange(-med_diff)
test_ll_burns_data_comp

## twCRPS aggregation and analysis ---------
twcrps_data_full <- holdout_twcrps_data %>%
  full_join(train_twcrps_data)

train_twcrps_data_ranked <- twcrps_data_full %>% filter(train == TRUE) %>% 
  group_by(model) %>% 
  summarize(mean_twcrps = mean(twcrps, na.rm = TRUE)) %>% 
  arrange(mean_twcrps)
top_mod_train_twcrps_data <- as.character(train_twcrps_data_ranked$model[1])

train_twcrps_data_comp <- twcrps_data_full %>% filter(train == TRUE) %>% 
  select(c(draw, model, twcrps)) %>% 
  pivot_wider(names_from = model, values_from = twcrps, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_train_twcrps_data))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(mean_diff = mean(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% 
  arrange(mean_diff)

test_twcrps_data_ranked <- twcrps_data_full %>% filter(train == FALSE) %>% 
  group_by(model) %>% 
  summarize(mean_twcrps = mean(twcrps, na.rm = TRUE)) %>% 
  arrange(mean_twcrps)
top_mod_test_twcrps_data <- as.character(test_twcrps_data_ranked$model[1])
test_twcrps_data_comp <- twcrps_data_full %>% filter(train == FALSE) %>% 
  select(c(draw, model, twcrps)) %>% 
  pivot_wider(names_from = model, values_from = twcrps, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_test_twcrps_data))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(mean_diff = mean(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% 
  arrange(mean_diff)
test_twcrps_data_comp


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
best_count_files <- count_files[grepl("zinb_er_pi-ri", count_files)]
best_count_files <- best_count_files[grepl("long", best_count_files)]
best_count_names <- (str_remove(basename(best_count_files), "_\\d{2}\\w{3}2023_\\d{4}_long_scores.RDS"))

for(i in seq_along(best_count_names)) {
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

## log-likelihood aggregating and analysis (burns) -------
ll_full_counts_data <- holdout_loglik_counts_data %>%
  full_join(train_loglik_counts_data)

train_ll_counts_ranked_data <- ll_full_counts_data %>% filter(train == TRUE) %>% 
  group_by(model) %>% 
  summarize(med_train_ll = median(loglik[is.finite(loglik)])) %>% 
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
