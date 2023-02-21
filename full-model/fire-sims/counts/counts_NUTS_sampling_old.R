# pass from command line ----
# Rscript code/scripts/r/SpRL_sim_study.R
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) stop("Pass in model (zinb, zinb_er, or zip), params (pi-reg or pi-ri)", call.=FALSE)
if (!(args[1] %in% c('zinb', 'zinb_er', 'zip'))) stop("Pass in the model specs (zinb, zinb_er, or zip)", call.=FALSE)
if (!(args[2] %in% c("pi-reg", "pi-ri"))) stop("Pass in the parameter combination (pi-reg or pi-ri)", call.=FALSE)

model <- args[1]
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

stan_data <- readRDS("./full-model/data/stan_data_og.rds")
egpd_init <- stan_model(paste0('./full-model/fire-sims/counts/', model, '/stan/', model, '_lambda-reg_', params, '.stan'))
egpd_fit <- sampling(egpd_init, 
                     data = stan_data, 
                     iter = 3000,
                     warmup = 1000,
                     thin = 2,
                     chains = 3)

end_time <- format(as.POSIXlt(Sys.time(), "America/Denver"), "%H%M")

# save MCMC object in case below dx plots don't save properly
saveRDS(egpd_fit, 
        file = paste0("./full-model/fire-sims/counts/", model, '/stan-fits/', model, '_', params, '_',
                      st_time, "_", end_time, ".RDS"))

# save traceplot
if(params == 'pi-reg') {
  MCMCtrace(egpd_fit, params = c("rho1_lambda", "rho2_lambda", "rho1_pi", "rho2_pi"),
            ind = TRUE,
            open_pdf = FALSE,
            filename = paste0('./full-model/figures/', model, '/trace/', model, '_', params, '_rhos_',
                              st_time, "_", end_time, ".pdf"))
  MCMCtrace(egpd_fit, params = c("beta_lambda", "beta_pi"),
            ind = TRUE,
            open_pdf = FALSE,
            filename = paste0('./full-model/figures/', model, '/trace/', model, '_', params, '_betas_',
                              st_time, "_", end_time, ".pdf"))
  MCMCtrace(egpd_fit, params = c("phi_lambda", "phi_pi"),
            ind = TRUE,
            open_pdf = FALSE,
            filename = paste0('./full-model/figures/', model, '/trace/', model, '_', params, '_phis_',
                              st_time, "_", end_time, ".pdf"))
} else {
  MCMCtrace(egpd_fit, params = c("rho1_lambda", "rho2_lambda", "rho1_pi", "rho2_pi"),
            ind = TRUE,
            open_pdf = FALSE,
            filename = paste0('./full-model/figures/', model, '/trace/', model, '_', params, '_rhos_',
                              st_time, "_", end_time, ".pdf"))
  MCMCtrace(egpd_fit, params = c("beta_lambda"),
            ind = TRUE,
            open_pdf = FALSE,
            filename = paste0('./full-model/figures/', model, '/trace/', model, '_', params, '_betas_',
                              st_time, "_", end_time, ".pdf"))
  MCMCtrace(egpd_fit, params = c("phi_lambda"),
            ind = TRUE,
            open_pdf = FALSE,
            filename = paste0('./full-model/figures/', model, '/trace/', model, '_', params, '_phis_',
                              st_time, "_", end_time, ".pdf"))
}
