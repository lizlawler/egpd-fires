
# Model comparisons for counts --------------------------------------------
library(cowplot)
library(tidyverse)
library(patchwork)
library(ggrepel)
library(rstan)

model_fits <- list.files(pattern = "*.RDS", recursive = TRUE)
count_fits <- grep(model_fits, pattern = '22-Oct-2022', value = TRUE, invert = FALSE)
count_fits <- grep(count_fits, pattern = "counts", value = TRUE, invert = FALSE)

# # data frame for plotting colors
# cols <- c('Poisson' = 'green3',
#           'Negative binomial' = 'skyblue',
#           'ZI Poisson' = 'orange',
#           'ZI Negative binomial' = 'hotpink1')

# for each model fit, produce a vector of the holdout log likelihoods
zip_scores <- rstan::extract(read_rds(count_fits[3]), pars = c("train_loglik", "holdout_loglik"))
zinb_deltaER_scores <- rstan::extract(read_rds(count_fits[2]), pars = c("train_loglik", "holdout_loglik"))
zinb_constant_delta_scores <- rstan::extract(read_rds(count_fits[1]), pars = c("train_loglik", "holdout_loglik"))

zip_train <- zip_scores$train_loglik %>% unlist() %>% bind_cols()

holdout_loglik <- list()
train_loglik <- list()
model_name <- c("zinb_constant-delta", "zinb_delta-byER", "zip")
                           
                           
for (i in seq_along(count_fits)) {
  post <- rstan::extract(read_rds(count_fits[i]))
  holdout_loglik[[i]] <- post$holdout_loglik %>%
    reshape2::melt(varnames = c('iter', 'idx')) %>%
    as_tibble %>%
    group_by(iter) %>%
    summarize(value = sum(value)) %>%
    mutate(model = count_fits[i])
  train_loglik[[i]] <- post$train_loglik %>%
    reshape2::melt(varnames = c('iter', 'idx')) %>%
    as_tibble %>%
    group_by(iter) %>%
    summarize(value = sum(value)) %>%
    mutate(model = count_fits[i])
  gc()
}

holdout_c_loglik <- bind_rows(holdout_c_loglik) %>%
  mutate(train = FALSE)

train_c_loglik <- bind_rows(train_c_loglik) %>%
  mutate(train = TRUE)

ll_d <- holdout_c_loglik %>%
  full_join(train_c_loglik) %>%
  mutate(train = ifelse(train == TRUE, 'train', 'test'),
         Distribution = case_when(
           grepl('^nb', .$model) ~ 'Negative binomial',
           grepl('pois', .$model) ~ 'Poisson',
           grepl('zip', .$model) ~ 'ZI Poisson',
           grepl('zinb', .$model) ~ 'ZI Negative binomial')) %>%
  spread(train, value)

loglik_df <- ll_d %>%
  group_by(Distribution) %>%
  summarize(mean_test = median(test),
            sd_test = sd(test)) %>%
  arrange(-mean_test) %>%
  mutate(pretty = paste0(format(mean_test, digits = 0, scientific = FALSE), 
                         ' (', trimws(format(sd_test, digits = 0, scientific = FALSE)), ')')) %>%
  rename(`Holdout log likelihood` = pretty, 
         Model = Distribution) %>%
  select(-ends_with('test')) 