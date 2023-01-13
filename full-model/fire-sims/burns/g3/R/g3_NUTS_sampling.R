# pass from command line ----
# Rscript code/scripts/r/SpRL_sim_study.R
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) stop("Pass in suffix (sqrt or og), params (delta-ri_xi-ri, nu-reg_xi-reg, nu-reg_xi-ri, nu-ri_xi-ri)", call.=FALSE)
if (!(args[1] %in% c('sqrt', 'og'))) stop("Pass in the response data type (sqrt or og)", call.=FALSE)
if (!(args[2] %in% c("delta-ri_xi-ri", "nu-reg_xi-reg", "nu-reg_xi-ri", "nu-ri_xi-ri"))) stop("Pass in the parameter combination (delta-ri_xi-ri, nu-reg_xi-reg, nu-reg_xi-ri, nu-ri_xi-ri)", call.=FALSE)

suffix <- args[1]
params <- args[2]

library(readr)
library(rstan)
library(MCMCvis)
library(tidyverse)
library(splines)
library(spdep)
library(sf)
library(RColorBrewer)
library(patchwork)
library(classInt)
library(spatialreg)
options(mc.cores = parallel::detectCores())

st_time <- format(as.POSIXlt(Sys.time(), "America/Denver"), "%d-%b-%Y_%H%M")

stan_data <- readRDS(paste0("./full-model/data/stan_data_", suffix, ".rds"))
egpd_init <- stan_model(paste0('./full-model/fire-sims/burns/g3/stan/g3_', params, '.stan'))
egpd_fit <- sampling(egpd_init, 
                     data = stan_data, 
                     iter = 3000,
                     warmup = 1000,
                     thin = 2,
                     init_r = 0.01,
                     chains = 3)

end_time <- format(as.POSIXlt(Sys.time(), "America/Denver"), "%H%M")

# save MCMC object in case below dx plots don't save properly
# post <- rstan::extract(egpd_fit, pars = c("beta", "phi", 
#                                           "rho1", "rho2", 
#                                           "delta", "xi", "sigma",
#                                           "holdout_loglik", "train_loglik"))
saveRDS(egpd_fit, file = paste0("./full-model/fire-sims/burns/g3/stan-fits/", params,"_loglik_", suffix, "_", 
                            st_time, "_", end_time, ".RDS"))

# save traceplot
if(params == 'nu-reg_xi-reg') {
  MCMCtrace(egpd_fit, params = c("rho1_delta", "rho2_delta", "rho1_nu", "rho2_nu", "rho1_xi", "rho2_xi"),
            ind = TRUE,
            open_pdf = FALSE,
            filename = paste0('./full-model/figures/g3/trace/g3_nu-reg_xi-reg_', suffix, '_rhos_',
                              st_time, "_", end_time, ".pdf"))
  MCMCtrace(egpd_fit, params = c("beta_delta", "beta_nu", "beta_xi"),
            ind = TRUE,
            open_pdf = FALSE,
            filename = paste0('./full-model/figures/g3/trace/g3_nu-reg_xi-reg_', suffix, '_betas_',
                              st_time, "_", end_time, ".pdf"))
  MCMCtrace(egpd_fit, params = c("phi_delta", "phi_nu", "phi_xi"),
            ind = TRUE,
            open_pdf = FALSE,
            filename = paste0('./full-model/figures/g3/trace/g3_nu-reg_xi-reg_', suffix, '_phis_',
                              st_time, "_", end_time, ".pdf"))
} else if (params == 'nu-reg_xi-ri') {
  MCMCtrace(egpd_fit, params = c("rho1_delta", "rho2_delta", "rho1_nu", "rho2_nu", "rho1_xi", "rho2_xi"),
            ind = TRUE,
            open_pdf = FALSE,
            filename = paste0('./full-model/figures/g3/trace/g3_nu-reg_xi-ri_', suffix, '_rhos_',
                              st_time, "_", end_time, ".pdf"))
  MCMCtrace(egpd_fit, params = c("beta_delta", "beta_nu"),
            ind = TRUE,
            open_pdf = FALSE,
            filename = paste0('./full-model/figures/g3/trace/g3_nu-reg_xi-ri_', suffix, '_betas_',
                              st_time, "_", end_time, ".pdf"))
  MCMCtrace(egpd_fit, params = c("phi_delta", "phi_nu"),
            ind = TRUE,
            open_pdf = FALSE,
            filename = paste0('./full-model/figures/g3/trace/g3_nu-reg_xi-ri_', suffix, '_phis_',
                              st_time, "_", end_time, ".pdf"))
} else if (params == 'delta-ri_xi-ri') {
  MCMCtrace(egpd_fit, params = c("rho1_delta", "rho2_delta", "rho1_nu", "rho2_nu", "rho1_xi", "rho2_xi"),
            ind = TRUE,
            open_pdf = FALSE,
            filename = paste0('./full-model/figures/g3/trace/g3_delta-ri_xi-ri_', suffix, '_rhos_',
                              st_time, "_", end_time, ".pdf"))
  MCMCtrace(egpd_fit, params = c("beta_nu"),
            ind = TRUE,
            open_pdf = FALSE,
            filename = paste0('./full-model/figures/g3/trace/g3_delta-ri_xi-ri_', suffix, '_betas_',
                              st_time, "_", end_time, ".pdf"))
  MCMCtrace(egpd_fit, params = c("phi_nu"),
            ind = TRUE,
            open_pdf = FALSE,
            filename = paste0('./full-model/figures/g3/trace/g3_delta-ri_xi-ri_', suffix, '_phis_',
                              st_time, "_", end_time, ".pdf"))
} else {
  MCMCtrace(egpd_fit, params = c("rho1_delta", "rho2_delta", "rho1_nu", "rho2_nu", "rho1_xi", "rho2_xi"),
            ind = TRUE,
            open_pdf = FALSE,
            filename = paste0('./full-model/figures/g3/trace/g3_nu-ri_xi-ri_', suffix, '_rhos_',
                              st_time, "_", end_time, ".pdf"))
  MCMCtrace(egpd_fit, params = c("beta_delta"),
            ind = TRUE,
            open_pdf = FALSE,
            filename = paste0('./full-model/figures/g3/trace/g3_nu-ri_xi-ri_', suffix, '_betas_',
                              st_time, "_", end_time, ".pdf"))
  MCMCtrace(egpd_fit, params = c("phi_delta"),
            ind = TRUE,
            open_pdf = FALSE,
            filename = paste0('./full-model/figures/g3/trace/g3_nu-ri_xi-ri_', suffix, '_phis_',
                              st_time, "_", end_time, ".pdf"))
}
