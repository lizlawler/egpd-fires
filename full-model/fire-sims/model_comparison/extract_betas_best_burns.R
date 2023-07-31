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
burn_fits <- burn_fits[which(grepl(pattern = "g1", burn_fits))]
burn_fits <- burn_fits[-which(grepl(pattern = "all-reg", burn_fits))]
nfits <- length(burn_fits)/3
fit_groups <- vector(mode = "list", nfits)
for(i in 1:nfits) {
  fit_groups[[i]] <- burn_fits[(3*i-2):(3*i)]
}

burn_names <- lapply(fit_groups, function(x) str_remove(basename(x[1]), "_\\d{2}\\w{3}2023_\\d{4}_\\d{1}.csv")) %>% unlist()

extraction <- function(file_group, model_name) {
  object <- as_cmdstan_fit(file_group)
  betas <- object$draws(variables = "beta")
  temp <- list(betas)
  names(temp) <- c("betas")
  assign(model_name, temp, parent.frame())
  rm(object)
  print(paste0(model_name, " complete"))
  gc()
}

for(i in seq_along(burn_names)) {
  extraction(fit_groups[[i]], burn_names[i])
}
# 
save(list=c(ls(pattern="g1"), ls(pattern = "burn_names")), 
     file = "full-model/fire-sims/model_comparison/g1_er_betas.RData")
