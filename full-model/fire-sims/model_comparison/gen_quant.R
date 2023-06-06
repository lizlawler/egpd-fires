library(cmdstanr)
library(posterior)
library(tidyverse)
library(stringr)
gq_files <- paste0("full-model/fire-sims/burns/g1/csv-fits/",
                   list.files(path="full-model/fire-sims/burns/g1/csv-fits/", pattern = "gq"))
nfits <- length(gq_files)/3
fit_groups <- vector(mode = "list", nfits)
for(i in 1:nfits) {
  fit_groups[[i]] <- gq_files[(3*i-2):(3*i)]
}
extraction <- function(file_group) {
  train_ll <- read_cmdstan_csv(file_group, variables = "train_loglik")
  holdout_ll <- read_cmdstan_csv(file_group, variables = "holdout_loglik")
  holdout_twcrps <- read_cmdstan_csv(file_group, variables = "holdout_twcrps")
  file <- basename(file_group[1])
  model <- str_remove(str_remove(str_remove(file, "_\\d{2}\\w{3}2023_\\d{4}_\\d{1}.csv"), "gq_"), "cfcns_")
  temp <- list(train_loglik = train_ll, holdout_loglik = holdout_ll, holdout_twcrps = holdout_twcrps)
  assign(model, temp, parent.frame())
  gc()
}

extraction(fit_groups[[20]])
save(list=ls(pattern="g1"), file = "full-model/fire-sims/model_comparison/g1_sqrt_xi-ri.RData")
rm(list=ls(pattern="g1"))
gc()

for(i in 1:nfits) {
  extraction(fit_groups[[i]])
}
