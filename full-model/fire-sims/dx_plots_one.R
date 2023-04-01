library(cmdstanr)
set_cmdstan_path(path = "/projects/eslawler@colostate.edu/.cmdstan/cmdstan-2.31.0") # this is only relevant to Alpine
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(brms)
library(MCMCvis)

type <- "burns"
model <- "g1"
suffix <- "og"
params <- "sigma-ri_xi-ri"
delta <- "0.81"

csvbase <- paste0("full-model/fire-sims/", type, "/", model, "/csv-fits/")
plotbase <- paste0("full-model/figures/", model, "/trace/")
csvpattern <- paste0(model, "_", suffix, "_", params, "_", delta)
csvfiles <- paste0(csvbase, list.files(path = csvbase, pattern = csvpattern))
csvfiles

model_fit <- brms:::read_csv_as_stanfit(csvfiles)
MCMCtrace(model_fit, params = "rho",
          ind = TRUE,
          open_pdf = FALSE, 
          filename = paste0(plotbase, csvpattern, "_rho.pdf"))
MCMCtrace(model_fit, params = "beta", 
          ind = TRUE,
          open_pdf = FALSE,
          filename = paste0(plotbase, csvpattern, "_beta.pdf"))
MCMCtrace(model_fit, params = "phi",
          ind = TRUE,
          open_pdf = FALSE,
          filename = paste0(plotbase, csvpattern, "_phi.pdf"))