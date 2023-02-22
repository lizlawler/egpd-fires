# pass from command line ----
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) stop("Pass in model (zinb, zinb_er, or zip), params (all-reg or pi-ri)", call.=FALSE)
if (!(args[1] %in% c('zinb', 'zinb_er', 'zip'))) stop("Pass in the model specs (zinb, zinb_er, or zip)", call.=FALSE)
if (!(args[2] %in% c("all-reg", "pi-ri"))) stop("Pass in the parameter combination (all-reg or pi-ri)", call.=FALSE)

model <- args[1]
params <- args[2]

library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(MCMCvis)
library(tidyverse)

st_time <- format(as.POSIXlt(Sys.time(), "America/Denver"), "%d-%b-%Y_%H%M")
stan_data <- readRDS("./full-model/data/stan_data_og.rds")
counts_model <- cmdstan_model(paste0('./full-model/fire-sims/counts/', model, '/stan/', model, '_', params, '.stan'), 
                             compile = TRUE, stanc_options = list("O1"))
counts_fit <- counts_model$sample(data = stan_data, 
                                  iter_warmup = 1000,
                                  iter_sampling = 2000,
                                  thin = 2,
                                  chains = 3,
                                  parallel_chains = 3,
                                  output_dir = paste0("full-model/fire-sims/counts/", model, '/csv-fits/'),
                                  output_basename = paste0(model, "_", params, "_", st_time))

end_time <- format(as.POSIXlt(Sys.time(), "America/Denver"), "%d%b%Y_%H%M")

counts_stan_fit <- rstan::read_stan_csv(counts_fit$output_files())
saveRDS(counts_stan_fit, 
        file = paste0("./full-model/fire-sims/counts/", model, '/stan-fits/', model, '_', params, '_',
                      st_time, "_", end_time, ".RDS"))

# save traceplot
MCMCtrace(counts_fit, params = "rho",
          ind = TRUE,
          open_pdf = FALSE,
          filename = paste0('./full-model/figures/', model, '/trace/', model, '_', params, '_rho_', st_time, "_", end_time, ".pdf"))
MCMCtrace(counts_fit, params = "beta", 
          ind = TRUE,
          open_pdf = FALSE,
          filename = paste0('./full-model/figures/', model, '/trace/', model, '_', params, '_beta_', st_time, "_", end_time, ".pdf"))
MCMCtrace(counts_fit, params = "phi", 
          ind = TRUE, 
          open_pdf = FALSE,
          filename = paste0('./full-model/figures/', model, '/trace/', model, '_', params, '_phi_', st_time, "_", end_time, ".pdf"))

