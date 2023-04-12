args <- commandArgs(trailingOnly=TRUE)
type <- args[1]
model <- args[2]
suffix <- args[3]
params <- args[4]
delta <- args[5]
sttime <- args[6]

library(cmdstanr)
set_cmdstan_path(path = "/projects/eslawler@colostate.edu/software/anaconda/envs/lawler/bin/cmdstan") # this is only relevant to Alpine
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(MCMCvis)

csvbase <- paste0("./full-model/fire-sims/", type, "/", model, "/csv-fits/")
plotbase <- paste0("./full-model/figures/", model, "/trace/")
csvpattern <- paste0(model, "_", suffix, "_", params, "_", delta, "_", sttime)
csvfiles <- paste0(csvbase, list.files(path = csvbase, pattern = csvpattern))

fit <- as_cmdstan_fit(csvfiles)
fitmcmc <- as_mcmc.list(fit)

MCMCtrace(fitmcmc,
          params = 'rho', 
          ind = TRUE, 
          open_pdf = FALSE, 
          filename = paste0(plotbase, csvpattern, "_rho.pdf"))
print("rho plots created")

MCMCtrace(fitmcmc,
          params = 'beta', 
          ind = TRUE, 
          open_pdf = FALSE, 
          filename = paste0(plotbase, csvpattern, "_beta.pdf"))
print("beta plots created")

MCMCtrace(fitmcmc,
          params = 'phi', 
          ind = TRUE, 
          open_pdf = FALSE, 
          filename = paste0(plotbase, csvpattern, "_phi.pdf"))
print("phi plots created")