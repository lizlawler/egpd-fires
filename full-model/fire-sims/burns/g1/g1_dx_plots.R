# pass from command line ----
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) stop("Pass in suffix (sqrt or og), 
                            params (all-reg, xi-ri, sigma-ri_xi-ri, nu-ri_xi-ri, kappa-ri_xi-ri)", call.=FALSE)
if (!(args[1] %in% c('sqrt', 'og'))) stop("Pass in the response data type (sqrt or og)", call.=FALSE)
if (!(args[2] %in% c("all-reg", "xi-ri", "sigma-ri_xi-ri", "nu-ri_xi-ri", "kappa-ri_xi-ri"))) 
  stop("Pass in the parameter combination (all-reg, xi-ri, sigma-ri_xi-ri, nu-ri_xi-ri, kappa-ri_xi-ri)", call.=FALSE)

suffix <- args[1]
params <- args[2]

library(cmdstanr)
set_cmdstan_path(path = "/projects/eslawler@colostate.edu/.cmdstan/cmdstan-2.31.0") # this is only relevant to Alpine
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(MCMCvis)
library(tidyverse)
library(stringr)
options(mc.cores = parallel::detectCores())

file_base <- paste0("g1_", suffix, "_", params, "*")
csv_files <- paste0("./full-model/fire-sims/burns/g1/csv-fits/", 
                    list.files("full-model/fire-sims/burns/g1/csv-fits/", 
                               pattern = file_base))
print(csv_files)
st_date <- stringr::str_extract(basename(csv_files[1]), "\\d{2}-Mar-\\d{4}")

egpd_fit <- as_cmdstan_fit(csv_files, format = "draws_list")
# egpd_fit$save_obj(file = paste0("./full-model/fire-sims/burns/g1/cmd-stan-fits/g1_", suffix, "_", params, "_", st_date, ".RDS"))
egpd_mcmc <- as_mcmc.list(egpd_fit)

# save traceplots
MCMCtrace(egpd_mcmc, params = "rho",
          ind = TRUE,
          open_pdf = FALSE,
          filename = paste0('./full-model/figures/g1/trace/g1_', params, '_', 
                            suffix, '_rhos_', st_date, ".pdf"))
MCMCtrace(egpd_mcmc, params = "beta",
          ind = TRUE,
          open_pdf = FALSE,
          filename = paste0('./full-model/figures/g1/trace/g1_', params, '_', 
                            suffix,'_betas_', st_date, ".pdf"))
MCMCtrace(egpd_mcmc, params = "phi",
          ind = TRUE,
          open_pdf = FALSE,
          filename = paste0('./full-model/figures/g1/trace/g1_', params, '_', 
                            suffix,'_phis_', st_date, ".pdf"))

