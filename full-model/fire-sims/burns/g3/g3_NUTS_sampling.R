# pass from command line ----
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) stop("Pass in suffix (sqrt or og), params (all-reg, xi-ri)", call.=FALSE)
if (!(args[1] %in% c('sqrt', 'og'))) stop("Pass in the response data type (sqrt or og)", call.=FALSE)
if (!(args[2] %in% c("all-reg", "xi-ri"))) 
  stop("Pass in the parameter combination (all-reg, xi-ri)", call.=FALSE)

suffix <- args[1]
params <- args[2]

library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(MCMCvis)
library(tidyverse)

st_time <- format(as.POSIXlt(Sys.time(), "America/Denver"), "%d-%b-%Y_%H%M")
stan_data <- readRDS(paste0("./full-model/data/stan_data_", suffix, ".rds"))
egpd_model <- cmdstan_model(paste0('./full-model/fire-sims/burns/g3/stan/g3_', params, '.stan'), compile = TRUE, stanc_options = list("O1"))
egpd_fit <- egpd_model$sample(data = stan_data, 
                              iter_warmup = 1000,
                              iter_sampling = 2000,
                              thin = 2,
                              chains = 3,
                              parallel_chains = 3,
                              init = 0.01,
                              output_dir = "full-model/fire-sims/burns/g3/csv-fits/")

end_time <- format(as.POSIXlt(Sys.time(), "America/Denver"), "%H%M")

# save CmdStanMCMC object
file_name <- paste0("./full-model/fire-sims/burns/g3/cmd-stan-fits/g3_", 
                    params, suffix, st_time, "_", end_time, ".RDS")
egpd_fit$save_object(file = file_name)

# convert CmdStanMCMC object to mcmc list for use in MCMCtrace
egpd_mcmc <- as_mcmc.list(egpd_fit)

# save traceplots
MCMCtrace(egpd_mcmc, params = "rho",
          ind = TRUE,
          open_pdf = FALSE,
          filename = paste0('./full-model/figures/g3/trace/g3_', params, '_', 
                            suffix, '_rhos_', st_time, "_", end_time, ".pdf"))
MCMCtrace(egpd_mcmc, params = "beta",
          ind = TRUE,
          open_pdf = FALSE,
          filename = paste0('./full-model/figures/g3/trace/g3_', params, '_', 
                            suffix,'_betas_', st_time, "_", end_time, ".pdf"))
MCMCtrace(egpd_mcmc, params = "phi",
          ind = TRUE,
          open_pdf = FALSE,
          filename = paste0('./full-model/figures/g3/trace/g3_', params, '_',
                            suffix,'_phis_', st_time, "_", end_time, ".pdf"))