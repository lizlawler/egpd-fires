args <- commandArgs(trailingOnly=TRUE)
type <- args[1]
model <- args[2]
params <- args[3]
dataset <- args[4]
sttime <- args[5]

library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(MCMCvis)
library(posterior)

csvbase <- paste0("./full-model/fire-sims/", type, "/", model, "/csv-fits/")
plotbase <- paste0("./full-model/figures/", type, "/", model, "/trace/")
csvpattern <- paste0( model, "_", dataset, "_", params)
csvfiles <- paste0(csvbase, list.files(path = csvbase, pattern = csvpattern))

csvfiles <- csvfiles[c(1,3)] # keeping chains that actually completed

print("Filenames being used are:")
csvfiles

print("Creating cmdstan fit object...")
fit <- as_cmdstan_fit(csvfiles)
fit$diagnostic_summary()

print("Creating mcmc list from cmdstan object...")
fitmcmc <- as_mcmc.list(fit)

options(bitmapType='cairo') # allows MCMCtrace to work on RStudio server

print("Creating traceplot of rhos...")
MCMCtrace(fitmcmc,
          params = c('rho1', 'rho2'), 
          ind = TRUE, 
          pdf = TRUE,
          open_pdf = FALSE,
          filename = paste0(plotbase, csvpattern, "_rho.pdf"))
print("rho traceplots created")

print("Creating traceplot of betas")
MCMCtrace(fitmcmc,
          params = 'beta', 
          ind = TRUE, 
          open_pdf = FALSE, 
          filename = paste0(plotbase, csvpattern, "_beta.pdf"))
print("beta traceplots created")

print("Creating traceplot of phi...")
MCMCtrace(fitmcmc,
          params = 'phi', 
          ind = TRUE, 
          open_pdf = FALSE, 
          filename = paste0(plotbase, csvpattern, "_phi.pdf"))
print("phi traceplots created")

print("Extracting scores from model object...")
if (type == "burns") {
  train_loglik <- fit$draws(variables = "train_loglik")
  holdout_loglik <- fit$draws(variables = "holdout_loglik")
  train_twcrps <- fit$draws(variables = "train_twcrps")
  holdout_twcrps <- fit$draws(variables = "holdout_twcrps")
  scores <- list(train_loglik, holdout_loglik, train_twcrps, holdout_twcrps)
  names(scores) <- c("train_loglik", "holdout_loglik", "train_twcrps", "holdout_twcrps")
} else {
  train_loglik <- fit$draws(variables = "train_loglik")
  holdout_loglik <- fit$draws(variables = "holdout_loglik")
  scores <- list(train_loglik, holdout_loglik)
  names(scores) <- c("train_loglik", "holdout_loglik")  
}

filename <- paste0("full-model/fire-sims/model_comparison/extracted_values/", csvpattern, "_scores.RDS")
saveRDS(scores, file = filename)
print("Scores have been extracted and saved to disk")
