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

csvbase <- paste0("./models/", type, "/csv_fits/")
plotbase <- paste0("./scores_traceplots/traceplots/", type, "/")
csvpattern <- paste0(type, "_", model, "_", params, "_", sttime)
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

print("Creating traceplot of betas for lambda")
MCMCtrace(fitmcmc,
          params = 'beta_count',
          ind = TRUE,
          open_pdf = FALSE,
          filename = paste0(plotbase, csvpattern, "_lambda-beta.pdf"))
print("Creating traceplot of betas for kappa")
MCMCtrace(fitmcmc,
          params = 'beta_size',
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

