args <- commandArgs(trailingOnly=TRUE)
type <- args[1]
model <- args[2]
params <- args[3]
dataset <- args[4]
sttime <- args[5]

library(cmdstanr)
set_cmdstan_path(path = "/projects/eslawler@colostate.edu/software/anaconda/envs/lawler/bin/cmdstan") # this is only relevant to Alpine
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(MCMCvis)
library(posterior)

csvbase <- paste0("./full-model/fire-sims/", type, "/csv-fits/")
csvpattern <- paste0(type, "_", model, "_", params, "_", sttime, "_", dataset, "_\\d{1}")
files <- paste0(csvbase, list.files(path = csvbase, pattern = csvpattern))
print("Filenames being used are:")
files

print("Creating cmdstan fit object...")
fit <- as_cmdstan_fit(files)
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
# 
# print("Creating traceplot of betas for lambda")
# MCMCtrace(fitmcmc,
#           params = 'beta_count',
#           ind = TRUE,
#           open_pdf = FALSE,
#           filename = paste0(plotbase, csvpattern, "_lambda-beta.pdf"))
# print("Creating traceplot of betas for kappa")
# MCMCtrace(fitmcmc,
#           params = 'beta_burn',
#           ind = TRUE,
#           open_pdf = FALSE,
#           filename = paste0(plotbase, csvpattern, "_kappa-beta.pdf"))
# print("beta traceplots created")
# 
# print("Creating traceplot of phi...")
# MCMCtrace(fitmcmc,
#           params = 'phi', 
#           ind = TRUE, 
#           open_pdf = FALSE, 
#           filename = paste0(plotbase, csvpattern, "_phi.pdf"))
# print("phi traceplots created")
# 
# print("Creating traceplot of ri_init...")
# MCMCtrace(fitmcmc,
#           params = 'ri_init',
#           ind = TRUE,
#           open_pdf = FALSE,
#           filename = paste0(plotbase, csvpattern, "_ri-init.pdf"))
# print("sigma and xi traceplots created")

print("Extracting scores from model object...")
holdout_loglik_count <- fit$draws(variables = "holdout_loglik_count")
holdout_loglik_burn <- fit$draws(variables = "holdout_loglik_burn")
holdout_twcrps <- fit$draws(variables = "holdout_twcrps")
scores <- list(holdout_loglik_count, holdout_loglik_burn, holdout_twcrps)
names(scores) <- c("holdout_loglik_count", "holdout_loglik_burn","holdout_twcrps")
filename <- paste0("full-model/fire-sims/model_comparison/extracted_values/", csvpattern, "_scores.RDS")
saveRDS(scores, file = filename)
print("Scores have been extracted and saved to disk")

