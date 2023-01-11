# pass from command line ----
# Rscript code/scripts/r/SpRL_sim_study.R
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) stop("Pass in suffix (sqrt or og), params (nu-reg_xi-reg, nu-reg_xi-ri, nu-ri_xi-ri)", call.=FALSE)
if (!(args[1] %in% c('sqrt', 'og'))) stop("Pass in the response data type (sqrt or oh)", call.=FALSE)
if (!(args[2] %in% c("nu-reg_xi-reg", "nu-reg_xi-ri", "nu-ri_xi-ri"))) stop("Pass in the parameter combination (nu-reg_xi-reg, nu-reg_xi-ri, nu-ri_xi-ri)", call.=FALSE)

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
egpd_init <- stan_model(paste0('./full-model/fire-sims/burns/g1/stan/g1_', params, '.stan'))
egpd_fit <- sampling(egpd_init, 
                     data = stan_data, 
                     iter = 3000,
                     warmup = 1000,
                     thin = 2,
                     init_r = 0.01,
                     chains = 3)

end_time <- format(as.POSIXlt(Sys.time(), "America/Denver"), "%H%M")

# save MCMC object in case below dx plots don't save properly
post <- rstan::extract(egpd_fit, pars = c("beta", "phi", "rho1", "rho2", "holdout_loglik", "train_loglik"))
saveRDS(post, 
        file = paste0("./full-model/fire-sims/burns/g1/stan-fits/g1_", params, suffix, 
                      st_time, "_", end_time, ".RDS"))

# save traceplot
if(params == 'nu-reg_xi-reg') {
  MCMCtrace(egpd_fit, params = c("rho1_kappa", "rho2_kappa", "rho1_nu", "rho2_nu", "rho1_xi", "rho2_xi"),
            ind = TRUE,
            open_pdf = FALSE,
            filename = paste0('./full-model/figures/g1/trace/g1_nu-reg_xi-reg_', suffix, '_rhos_',
                              st_time, "_", end_time, ".pdf"))
  MCMCtrace(egpd_fit, params = c("beta_kappa", "beta_nu", "beta_xi"),
            ind = TRUE,
            open_pdf = FALSE,
            filename = paste0('./full-model/figures/g1/trace/g1_nu-reg_xi-reg_', suffix, '_betas_',
                              st_time, "_", end_time, ".pdf"))
  MCMCtrace(egpd_fit, params = c("phi_kappa", "phi_nu", "phi_xi"),
            ind = TRUE,
            open_pdf = FALSE,
            filename = paste0('./full-model/figures/g1/trace/g1_nu-reg_xi-reg_', suffix, '_phis_',
                              st_time, "_", end_time, ".pdf"))
} else if (params == 'nu-reg_xi-ri') {
  MCMCtrace(egpd_fit, params = c("rho1_kappa", "rho2_kappa", "rho1_nu", "rho2_nu", "rho1_xi", "rho2_xi"),
            ind = TRUE,
            open_pdf = FALSE,
            filename = paste0('./full-model/figures/g1/trace/g1_nu-reg_xi-ri_', suffix, '_rhos_',
                              st_time, "_", end_time, ".pdf"))
  MCMCtrace(egpd_fit, params = c("beta_kappa", "beta_nu"),
            ind = TRUE,
            open_pdf = FALSE,
            filename = paste0('./full-model/figures/g1/trace/g1_nu-reg_xi-ri_', suffix, '_betas_',
                              st_time, "_", end_time, ".pdf"))
  MCMCtrace(egpd_fit, params = c("phi_kappa", "phi_nu"),
            ind = TRUE,
            open_pdf = FALSE,
            filename = paste0('./full-model/figures/g1/trace/g1_nu-reg_xi-ri_', suffix, '_phis_',
                              st_time, "_", end_time, ".pdf"))
} else {
  MCMCtrace(egpd_fit, params = c("rho1_kappa", "rho2_kappa", "rho1_nu", "rho2_nu", "rho1_xi", "rho2_xi"),
            ind = TRUE,
            open_pdf = FALSE,
            filename = paste0('./full-model/figures/g1/trace/g1_nu-ri_xi-ri_', suffix, '_rhos_',
                              st_time, "_", end_time, ".pdf"))
  MCMCtrace(egpd_fit, params = c("beta_kappa"),
            ind = TRUE,
            open_pdf = FALSE,
            filename = paste0('./full-model/figures/g1/trace/g1_nu-ri_xi-ri_', suffix, '_betas_',
                              st_time, "_", end_time, ".pdf"))
  MCMCtrace(egpd_fit, params = c("phi_kappa"),
            ind = TRUE,
            open_pdf = FALSE,
            filename = paste0('./full-model/figures/g1/trace/g1_nu-ri_xi-ri_', suffix, '_phis_',
                              st_time, "_", end_time, ".pdf"))
}
