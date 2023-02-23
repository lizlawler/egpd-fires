suffix <- "og"
params <- "all-reg"

library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(MCMCvis)
library(tidyverse)

st_time <- format(as.POSIXlt(Sys.time(), "America/Denver"), "%d-%b-%Y_%H%M")
stan_data <- readRDS(paste0("./full-model/data/stan_data_", suffix, ".rds"))
egpd_model <- cmdstan_model(paste0('./full-model/fire-sims/burns/g1/stan/g1_', params, '.stan'), compile = TRUE, stanc_options = list("O1"))
egpd_fit <- egpd_model$sample(data = stan_data, 
                              iter_warmup = 1000,
                              iter_sampling = 2000,
                              thin = 2,
                              chains = 3,
                              parallel_chains = 3,
                              init = 0.01,
                              output_dir = "full-model/fire-sims/burns/g1/csv-fits/")

end_time <- format(as.POSIXlt(Sys.time(), "America/Denver"), "%H%M")

# save CmdStanMCMC object
file_name <- paste0("./full-model/fire-sims/burns/g1/cmd-stan-fits/g1_", 
                    params, suffix, st_time, "_", end_time, ".RDS")
egpd_fit$save_object(file = file_name)

# convert CmdStanMCMC object to mcmc list for use in MCMCtrace
egpd_mcmc <- as_mcmc.list(egpd_fit)

# save traceplots
MCMCtrace(egpd_mcmc, params = "rho",
          ind = TRUE,
          open_pdf = FALSE,
          filename = paste0('./full-model/figures/g1/trace/g1_', params, '_', 
                            suffix, '_rhos_', st_time, "_", end_time, ".pdf"))
  MCMCtrace(egpd_mcmc, params = "beta",
            ind = TRUE,
            open_pdf = FALSE,
            filename = paste0('./full-model/figures/g1/trace/g1_', params, '_', 
                              suffix,'_betas_', st_time, "_", end_time, ".pdf"))
  MCMCtrace(egpd_mcmc, params = "phi",
            ind = TRUE,
            open_pdf = FALSE,
            filename = paste0('./full-model/figures/g1/trace/g1_', params, '_', 
                              suffix,'_phis_', st_time, "_", end_time, ".pdf"))
