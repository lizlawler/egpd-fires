args <- commandArgs(trailingOnly=TRUE)
csvpattern <- args[1]

library(cmdstanr)
set_cmdstan_path(path = "/projects/eslawler@colostate.edu/software/anaconda/envs/lawler/bin/cmdstan") # this is only relevant to Alpine
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(MCMCvis)

csvbase <- "./full-model/fire-sims/burns/g1/csv-fits/"
plotbase <- "./full-model/figures/g1/trace/"
csvfiles <- paste0(csvbase, list.files(path = csvbase, pattern = csvpattern))
csvfiles

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