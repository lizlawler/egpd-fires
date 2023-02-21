# pass from command line ----
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) stop("Pass in suffix (sqrt or og), params (all-reg, xi-ri, sigma-ri_xi-ri, nu-ri_xi-ri, kappa-ri_xi-ri)", call.=FALSE)
if (!(args[1] %in% c('sqrt', 'og'))) stop("Pass in the response data type (sqrt or og)", call.=FALSE)
if (!(args[2] %in% c("all-reg", "xi-ri", "sigma-ri_xi-ri", "nu-ri_xi-ri", "kappa-ri_xi-ri"))) 
  stop("Pass in the parameter combination (all-reg, xi-ri, sigma-ri_xi-ri, nu-ri_xi-ri, kappa-ri_xi-ri)", call.=FALSE)

suffix <- args[1]
params <- args[2]

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

egpd_stan_fit <- rstan::read_stan_csv(egpd_fit$output_files())

saveRDS(egpd_stan_fit, 
        file = paste0("./full-model/fire-sims/burns/g1/stan-fits/g1_", params, suffix, 
                      st_time, "_", end_time, ".RDS"))

# save traceplots
MCMCtrace(egpd_stan_fit, params = "rho",
          ind = TRUE,
          open_pdf = FALSE,
          filename = paste0('./full-model/figures/g1/trace/g1_', params, '_', suffix, '_rhos_',
                              st_time, "_", end_time, ".pdf"))
  MCMCtrace(egpd_stan_fit, params = "beta",
            ind = TRUE,
            open_pdf = FALSE,
            filename = paste0('./full-model/figures/g1/trace/g1_', params, '_', suffix,'_betas_',
                              st_time, "_", end_time, ".pdf"))
  MCMCtrace(egpd_stan_fit, params = "phi",
            ind = TRUE,
            open_pdf = FALSE,
            filename = paste0('./full-model/figures/g1/trace/g1_', params, '_', suffix,'_phis_',
                              st_time, "_", end_time, ".pdf"))
