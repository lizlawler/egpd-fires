args <- commandArgs(trailingOnly=TRUE)
type <- args[1]
model <- args[2]
params <- args[3]
sttime <- args[4]

library(cmdstanr)
set_cmdstan_path(path = "/projects/eslawler@colostate.edu/software/anaconda/envs/lawler/bin/cmdstan") # this is only relevant to Alpine
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(MCMCvis)
library(posterior)

csvbase <- paste0("./full-model/fire-sims/", type, "/csv-fits/")
plotbase <- paste0("./full-model/figures/", type, "/trace/")
csvpattern <- paste0(type, "_", model, "_", params, "_", sttime, "_erc_fwi")
csvfiles <- paste0(csvbase, list.files(path = csvbase, pattern = csvpattern))

print("Filenames being used are:")
csvfiles

print("Creating cmdstan fit object...")
fit <- as_cmdstan_fit(csvfiles)
fit$diagnostic_summary()

print("Creating mcmc list from cmdstan object...")
fitmcmc <- as_mcmc.list(fit)

print("Creating traceplot of rhos...")
MCMCtrace(fitmcmc,
          params = c('rho1', 'rho2'), 
          ind = TRUE, 
          open_pdf = FALSE, 
          filename = paste0(plotbase, csvpattern, "_rho.pdf"))
print("rho traceplots created")

print("Creating traceplot of betas for lambda")
MCMCtrace(fitmcmc,
          params = 'beta_count', 
          ind = TRUE, 
          open_pdf = FALSE, 
          filename = paste0(plotbase, csvpattern, "_lambda-beta.pdf"))
print("Creating traceplot of betas for kappa")
MCMCtrace(fitmcmc,
          params = 'beta_burn', 
          ind = TRUE, 
          open_pdf = FALSE, 
          filename = paste0(plotbase, csvpattern, "_kappa-beta.pdf"))
print("beta traceplots created")

print("Creating traceplot of phi...")
MCMCtrace(fitmcmc,
          params = 'phi', 
          ind = TRUE, 
          open_pdf = FALSE, 
          filename = paste0(plotbase, csvpattern, "_phi.pdf"))
print("phi traceplots created")

print("Extracting scores from model object...")
train_loglik_count <- fit$draws(variables = "train_loglik_count")
holdout_loglik_count <- fit$draws(variables = "holdout_loglik_count")
train_loglik_burn <- fit$draws(variables = "train_loglik_burn")
holdout_loglik_burn <- fit$draws(variables = "holdout_loglik_burn")
train_twcrps <- fit$draws(variables = "train_twcrps")
holdout_twcrps <- fit$draws(variables = "holdout_twcrps")
scores <- list(train_loglik_count, holdout_loglik_count, train_loglik_burn, 
               holdout_loglik_burn, train_twcrps, holdout_twcrps)
names(scores) <- c("train_loglik_count", "holdout_loglik_count", "train_loglik_burn", 
                   "holdout_loglik_burn","train_twcrps", "holdout_twcrps")
filename <- paste0("full-model/fire-sims/model_comparison/extracted_values/", csvpattern, "_scores.RDS")
saveRDS(scores, file = filename)
print("Scores have been extracted and saved to disk")

