library(cmdstanr)
set_cmdstan_path(path = "/projects/eslawler@colostate.edu/software/anaconda/envs/lawler/bin/cmdstan") # this is only relevant to Alpine
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(magrittr)
library(stringr)
library(posterior)

# following code is for extracting from the actual model fit -------
count_fits <- paste0("full-model/fire-sims/counts/",
                     list.files(path = "full-model/fire-sims/counts/",
                                pattern = "Jul2023", recursive = TRUE))
count_fits <- count_fits[which(grepl(pattern = "zinb_er", count_fits))]
nfits <- length(count_fits)/3
fit_groups <- vector(mode = "list", nfits)
for(i in 1:nfits) {
  fit_groups[[i]] <- count_fits[(3*i-2):(3*i)]
}

count_name <- lapply(fit_groups, function(x) str_remove(basename(x[1]), "_\\d{2}\\w{3}2023_\\d{4}_\\d{1}.csv")) %>% unlist()

extraction <- function(file_group, count_name) {
  object <- as_cmdstan_fit(file_group)
  betas <- object$draws(variables = "beta")
  temp <- list(betas)
  names(temp) <- c("betas")
  assign(burn_name, temp, parent.frame())
  rm(object)
  gc()
}

for(i in seq_along(count_names)) {
  extraction(fit_groups[[i]], count_names[i])
}
# 
save(list=c(ls(pattern="zinb_er"), ls(pattern = "count_names")), 
     file = "full-model/fire-sims/model_comparison/zinb_er_betas.RData")
