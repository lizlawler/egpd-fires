library(cmdstanr)
set_cmdstan_path(path = "/projects/eslawler@colostate.edu/software/anaconda/envs/lawler/bin/cmdstan") # this is only relevant to Alpine
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(magrittr)
library(stringr)
library(posterior)

# following code is for extracting from the actual model fit -------
burn_fits <- paste0("full-model/fire-sims/burns/",
                    list.files(path = "full-model/fire-sims/burns/",
                               pattern = "Jul2023", recursive = TRUE))
burn_fits <- burn_fits[-which(grepl("g2_", burn_fits))]
nfits <- length(burn_fits)/3
fit_groups <- vector(mode = "list", nfits)
for(i in 1:nfits) {
  fit_groups[[i]] <- burn_fits[(3*i-2):(3*i)]
}

burn_names <- lapply(fit_groups, function(x) str_remove(basename(x[1]), "_\\d{2}\\w{3}2023_\\d{4}_\\d{1}.csv")) %>% unlist()

# model scoring
extraction <- function(file_group, model_name) {
  model_object <- as_cmdstan_fit(file_group)
  train_loglik <- model_object$draws(variables = "train_loglik")
  holdout_loglik <- model_object$draws(variables = "holdout_loglik")
  train_twcrps <- model_object$draws(variables = "train_twcrps")
  holdout_twcrps <- model_object$draws(variables = "holdout_twcrps")
  rm(model_object)
  temp <- list(train_loglik, holdout_loglik, train_twcrps, holdout_twcrps)
  names(temp) <- c("train_loglik", "holdout_loglik", "train_twcrps", "holdout_twcrps")
  assign(model_name, temp, parent.frame())
  print(paste0(model_name, " complete"))
  gc()
}

for(i in 1:nfits) {
  extraction(fit_groups[[i]], burn_names[i])
}

save(list=c(ls(pattern="g1"), ls(pattern="lognorm"), ls(pattern = "burn_names")), 
     file = "full-model/fire-sims/model_comparison/gq_newdata_burns.RData")
