library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(tidyverse)
library(stringr)
library(posterior)

# load in scores extracted at time of model run on Alpine
score_files <- paste0("full-model/fire-sims/model_comparison/extracted_values/",
                      list.files("full-model/fire-sims/model_comparison/extracted_values/", 
                                 "g1"))
# process scores from g1 runs only
g1_files <- score_files[grepl("scores.RDS", score_files)]
g1_climate <- g1_files[grepl("climate", g1_files)]

# process scores only from "normal" qos - run for 24 hrs or less
# score_files <- score_files[grepl("normal", score_files)]
# erc_burn_files <- score_files[!grepl("fwi", score_files)]
# erc_burn_files <- erc_burn_files[!grepl("g3", erc_burn_files)]
# erc_burn_files <- erc_burn_files[!grepl("g4", erc_burn_files)]
# erc_burn_files <- erc_burn_files[grepl("long", erc_burn_files)]
# erc_burn_files <- erc_burn_files[grepl("GQ", erc_burn_files)]
# g1_burn_files <- erc_burn_files[grepl("g1", erc_burn_files)]
# g2_files <- erc_burn_files[grepl("nu_erc", erc_burn_files)]
# lognorm_files <- erc_burn_files[grepl("lognorm", erc_burn_files)]
# burn_files <- c(g2_files, lognorm_files)
# 
# erc_burn_files <- erc_burn_files[grepl("14Sep", erc_burn_files)]
# g2_gq_files <- score_files[grepl("g1", score_files)]
# g2_gq_files <- g1_gq_files[grepl("_GQ", g1_gq_files)]
# g2_gq_files <- g1_gq_files[grepl("long", g1_gq_files)]
# 
# burn_files <- c(score_files[grepl("g1", score_files)], 
#                 score_files[grepl("lognorm", score_files)], 
#                 score_files[grepl("g2", score_files)])
# burn_files_climate <- burn_files[grepl("climate", burn_files)]
# count_files <- paste0("full-model/fire-sims/model_comparison/extracted_values/",
#                       list.files("full-model/fire-sims/model_comparison/extracted_values/", 
#                                  "scores"))
# count_files <- c(count_files[grepl("zinb", count_files)], count_files[grepl("zip", count_files)])
# count_files_climate <- count_files[grepl("climate", count_files)]
  

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

g1_climate_names <- str_remove(str_remove(str_remove(basename(g1_climate), "_climate_\\d{2}\\w{3}2023_\\d{4}"), "_scores.RDS"), "g1_")

for(i in seq_along(g1_climate_names)) {
  read_scores(g1_climate[i], g1_climate_names[i])
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

nfits <- length(burn_files)/3
fit_groups <- vector(mode = "list", nfits)
for(i in 1:nfits) {
  fit_groups[[i]] <- burn_files[(3*i-2):(3*i)]
}

model_names <- lapply(fit_groups, function(x) {
  basename(x[[1]]) %>% 
    str_remove(., "_erc_\\d{2}\\w{3}2023_\\d{4}_long_\\d{1}_GQ.csv")
  }) %>% 
  unlist()


# pull scores from each BURN model and aggregate -----
for(i in seq_along(model_names)) {
  extraction(fit_groups[[i]], model_names[i])
}



model_names <- c("g1_sigma-ri_xi-ri", model_names)
extraction(g1_burn_files, model_names[1])

# read in G3 model that finished and is analogous to G1 "best" model
g3_files <- paste0("full-model/fire-sims/burns/g3/csv-fits/",
                      list.files("full-model/fire-sims/burns/g3/csv-fits/",
                                 "erc_xi-ri_nu"))
g3_files <- g3_files[c(1,3)]
g3_name <- str_remove(basename(g3_files[1]), "_\\d{2}Sep2023_\\d{4}_\\d{1}.csv")
extraction(list(g3_files)[[1]], g3_name)

burn_model_names <- c(g3_name, "g4_erc_sigma-ri_xi-ri_nu_scores")

g3_holdout_loglik <- `g3_erc_xi-ri_nu`$holdout_loglik$post_warmup_draws %>%
  as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>% 
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  group_by(draw) %>% 
  summarize(loglik = sum(value)) %>%
  mutate(model = g3_name,
         train = FALSE)

g2_chains <- `g2_sigma-ri_xi-ri_nu`$holdout_loglik$generated_quantities %>%
  as_draws_df() %>%
  select(-c(".iteration")) %>% 
  pivot_longer(cols = !c(".draw",".chain")) %>%
  rename(draw = ".draw", chain = ".chain") %>%
  group_by(chain, draw) %>% 
  summarize(loglik = sum(value)) %>% ungroup()

unique((g2_chains %>% filter(chain == 1))$loglik) # chain 1 didn't mix

g4_holdout_loglik <- `g4_erc_sigma-ri_xi-ri_nu_scores`$holdout_loglik %>% 
  as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>% 
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  group_by(draw) %>% 
  summarize(loglik = sum(value)) %>%
  mutate(model = "g4_erc_sigma-ri_xi-ri_nu",
         train = FALSE)

g1_holdout_loglik <- g1_data_loglik %>% filter(model == "erc", train == FALSE) %>% mutate(model = 'g1_sigma-xi_xi-ri_erc')
g1_g3_g4_erc_loglik <- rbind(g1_holdout_loglik, g3_holdout_loglik, g4_holdout_loglik)

test_ll_erc_ranked <- g1_g3_g4_erc_loglik %>% 
  group_by(model) %>% 
  summarize(med_train_ll = median(loglik)) %>% 
  arrange(-med_train_ll)
top_mod_erc_test <- as.character(test_ll_erc_ranked$model[1])
test_ll_erc_comp <- g1_g3_g4_erc_loglik %>% 
  select(c(draw, model, loglik)) %>% 
  pivot_wider(names_from = model, values_from = loglik, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_erc_test))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(med_diff = mean(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% 
  arrange(-med_diff)
test_ll_erc_comp

g3_holdout_twcrps <- `g3_erc_xi-ri_nu`$holdout_twcrps$post_warmup_draws %>%
  as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>% 
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  group_by(draw) %>% 
  summarize(twcrps = mean(value)) %>%
  mutate(model = g3_name,
         train = FALSE)

g4_holdout_twcrps <- `g4_erc_sigma-ri_xi-ri_nu_scores`$holdout_twcrps %>% 
  as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>% 
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  group_by(draw) %>% 
  summarize(twcrps = mean(value)) %>%
  mutate(model = "g4_erc_sigma-ri_xi-ri_nu",
         train = FALSE)

g1_holdout_twcrps <- g1_data_twcrps %>% filter(model == "erc", train == FALSE) %>% mutate(model = 'g1_sigma-xi_xi-ri_erc')
g1_g3_g4_erc_twcrps <- rbind(g1_holdout_twcrps, g3_holdout_twcrps, g4_holdout_twcrps)

## twCRPS aggregation and analysis ---------
test_twcrps_ranked <- g1_g3_g4_erc_twcrps %>% 
  group_by(model) %>% 
  summarize(mean_twcrps = mean(twcrps, na.rm = TRUE)) %>% 
  arrange(mean_twcrps)
top_mod_test_twcrps <- as.character(test_twcrps_ranked$model[1])
test_twcrps_comp <- g1_g3_g4_erc_twcrps %>% 
  select(c(draw, model, twcrps)) %>% mutate(twcrps = twcrps * 10^2) %>%
  pivot_wider(names_from = model, values_from = twcrps, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_test_twcrps))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(mean_diff = mean(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% 
  arrange(mean_diff)
test_twcrps_comp

nfits <- length(g1_climate_names)
holdout_loglik <- vector("list", nfits)
# train_loglik <- vector("list", nfits+2)
# train_twcrps <- vector("list", nfits+2)
holdout_twcrps <- vector("list", nfits)
for(i in seq_along(g1_climate_names)) {
  holdout_loglik[[i]] <- get(g1_climate_names[i])[["holdout_loglik"]] %>%
    as_draws_df() %>%
    select(-c(".iteration", ".chain")) %>%
    pivot_longer(cols = !".draw") %>%
    rename(draw = ".draw") %>%
    group_by(draw) %>%
    summarize(loglik = sum(value)) %>%
    mutate(model = g1_climate_names[i],
           train = FALSE)
  holdout_twcrps[[i]] <- get(g1_climate_names[i])[["holdout_twcrps"]] %>%
    as_draws_df() %>%
    select(-c(".iteration", ".chain")) %>%
    pivot_longer(cols = !".draw") %>%
    rename(draw = ".draw") %>%
    group_by(draw) %>%
    summarize(twcrps = mean(value)) %>%
    mutate(model = g1_climate_names[i],
           train = FALSE)
  # train_loglik[[i]] <- get(model_names[i])[["train_loglik"]]$generated_quantities %>%
  #   as_draws_df() %>%
  #   select(-c(".iteration", ".chain")) %>%
  #   pivot_longer(cols = !".draw") %>%
  #   rename(draw = ".draw") %>%
  #   group_by(draw) %>%
  #   summarize(loglik = sum(value)) %>%
  #   mutate(model = model_names[i],
  #          train = TRUE)
  # train_twcrps[[i]] <- get(model_names[i])[["train_twcrps"]]$generated_quantities %>%
  #   as_draws_df() %>%
  #   select(-c(".iteration", ".chain")) %>%
  #   pivot_longer(cols = !".draw") %>%
  #   rename(draw = ".draw") %>%
  #   group_by(draw) %>%
  #   summarize(twcrps = mean(value)) %>%
  #   mutate(model = model_names[i],
  #          train = TRUE)
}

# add in g3 scores
g3_name <- 'g3_erc_xi-ri_nu'
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
# train_loglik[[4]] <- get(g3_name)[["train_loglik"]]$post_warmup_draws %>%
#   as_draws_df() %>%
#   select(-c(".iteration", ".chain")) %>%
#   pivot_longer(cols = !".draw") %>%
#   rename(draw = ".draw") %>%
#   group_by(draw) %>%
#   summarize(loglik = sum(value)) %>%
#   mutate(model = g3_name,
#          train = TRUE)
# train_twcrps[[4]] <- get(g3_name)[["train_twcrps"]]$post_warmup_draws %>%
#   as_draws_df() %>%
#   select(-c(".iteration", ".chain")) %>%
#   pivot_longer(cols = !".draw") %>%
#   rename(draw = ".draw") %>%
#   group_by(draw) %>%
#   summarize(twcrps = mean(value)) %>%
#   mutate(model = g3_name,
#          train = TRUE)

g2_name <- "g2_sigma-ri_xi-ri_nu"
holdout_loglik[[2]] <- get(g2_name)[["holdout_loglik"]]$generated_quantities %>%
  as_draws_df() %>%
  select(-c(".iteration")) %>%
  pivot_longer(cols = !c(".draw", ".chain")) %>%
  rename(draw = ".draw", chain = ".chain") %>%
  filter(chain != 1) %>% select(-chain) %>%
  group_by(draw) %>%
  summarize(loglik = sum(value)) %>%
  mutate(model = g2_name,
         train = FALSE)
holdout_twcrps[[2]] <- get(g2_name)[["holdout_twcrps"]]$generated_quantities %>%
  as_draws_df() %>%
  select(-c(".iteration")) %>%
  pivot_longer(cols = !c(".draw", ".chain")) %>%
  rename(draw = ".draw", chain = ".chain") %>%
  filter(chain != 1) %>% select(-chain) %>%
  group_by(draw) %>%
  summarize(twcrps = mean(value)) %>%
  mutate(model = g2_name,
         train = FALSE)
train_loglik[[1]] <- get(g2_name)[["train_loglik"]]$generated_quantities %>%
  as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>%
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  group_by(draw) %>%
  summarize(loglik = sum(value)) %>%
  mutate(model = g2_name,
         train = TRUE)
train_twcrps[[1]] <- get(g2_name)[["train_twcrps"]]$generated_quantities %>%
  as_draws_df() %>%
  select(-c(".iteration", ".chain")) %>%
  pivot_longer(cols = !".draw") %>%
  rename(draw = ".draw") %>%
  group_by(draw) %>%
  summarize(twcrps = mean(value)) %>%
  mutate(model = g2_name,
         train = TRUE)


# holdout_loglik_burns <- bind_rows(holdout_loglik) %>% rbind(g4_holdout_loglik)
holdout_loglik_burns <- bind_rows(holdout_loglik)
# train_loglik_burns <- bind_rows(train_loglik)
# holdout_twcrps_df <- bind_rows(holdout_twcrps) %>% rbind(g4_holdout_twcrps)
holdout_twcrps_df <- bind_rows(holdout_twcrps)
# train_twcrps <- bind_rows(train_twcrps)

## log-likelihood aggregating and analysis (burns) -------
ll_full_burns <- holdout_loglik_burns %>%
  full_join(train_loglik_burns) %>% filter(grepl("long", model))

train_ll_burns_ranked <- ll_full_burns %>% filter(grepl("long", model)) %>%
  filter(train == TRUE) %>% 
  group_by(model) %>% 
  summarize(med_train_ll = median(loglik[is.finite(loglik)])) %>% 
  arrange(-med_train_ll)
top_mod_burns_train <- as.character(train_ll_burns_ranked$model[1])

train_ll_burns_comp <- ll_full_burns %>% filter(grepl("long", model)) %>%
  filter(train == TRUE) %>% 
  select(c(draw, model, loglik)) %>% 
  pivot_wider(names_from = model, values_from = loglik, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_burns_train))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(med_diff = median(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% 
  arrange(-med_diff)

test_ll_burns_ranked <- holdout_loglik_burns %>%
  group_by(model) %>% 
  summarize(med_train_ll = median(loglik)) %>% 
  arrange(-med_train_ll)
top_mod_burns_test <- as.character(test_ll_burns_ranked$model[1])
test_ll_burns_comp <- holdout_loglik_burns %>%
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

test_twcrps_ranked <- holdout_twcrps_df %>% 
  group_by(model) %>% 
  summarize(mean_twcrps = mean(twcrps, na.rm = TRUE)) %>% 
  arrange(mean_twcrps)
top_mod_test_twcrps <- as.character(test_twcrps_ranked$model[1])
test_twcrps_comp <- holdout_twcrps_df %>% mutate(twcrps = twcrps * 10^2) %>%
  select(c(draw, model, twcrps)) %>% 
  pivot_wider(names_from = model, values_from = twcrps, values_fill = NA) %>% 
  mutate(across(.cols = -draw, ~ .x - get(top_mod_test_twcrps))) %>% 
  pivot_longer(cols = -draw, names_to = "model") %>%
  group_by(model) %>%
  summarize(mean_diff = mean(value[is.finite(value)]), sd_diff = sd(value[is.finite(value)])) %>% 
  arrange(mean_diff)
test_twcrps_comp

burn_scores <- list(twcrps = twcrps_full, loglik = ll_full_burns)
saveRDS(burn_scores, "full-model/fire-sims/model_comparison/burn_scores.RDS")
# pull scores from G1 "best" model run on different datasets
g1_files <- paste0("full-model/fire-sims/model_comparison/extracted_values/",
                   list.files("full-model/fire-sims/model_comparison/extracted_values/", 
                              "g1_sigma-ri"))
g1_files <- g1_files[!grepl("RDS", g1_files)]
g1_files <- g1_files[grepl("long", g1_files)]

nfits <- length(g1_files)/3
fit_groups <- vector(mode = "list", nfits)
for(i in 1:nfits) {
  fit_groups[[i]] <- g1_files[(3*i-2):(3*i)]
}

g1_model_names <- lapply(fit_groups, function(x) {
  basename(x[[1]]) %>% 
    str_remove(., "_\\d{2}\\w{3}2023_\\d{4}") %>%
    str_remove(., "_\\d{1}_GQ.csv") %>%
    str_remove(., "g1_sigma-ri_xi-ri_") %>%
    str_remove(., "_long")
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

holdout_loglik <- vector("list", nfits)
train_loglik <- vector("list", nfits)
train_twcrps <- vector("list", nfits)
holdout_twcrps <- vector("list", nfits)
for(i in seq_along(g1_model_names)) {
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

saveRDS(ll_full_counts, file = "full-model/fire-sims/model_comparison/count_scores.RDS")

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

## log-likelihood aggregating and analysis (counts) -------
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

saveRDS(ll_full_counts_data, "full-model/fire-sims/model_comparison/count_scores_dataset.RDS")

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
