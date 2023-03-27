args <- commandArgs(trailingOnly=TRUE)
type <- args[1]
model <- args[2]
suffix <- args[3]
params <- args[4]
delta <- args[5]

library(brms)
library(MCMCvis)

csvbase <- paste0("full-model/fire-sims/", type, "/", model, "/csv-fits/")
plotbase <- paste0("full-model/figures/", model, "/trace/")
csvpattern <- paste0(model, "_", suffix, "_", params, "_", delta)
csvfiles <- paste0(csvbase, list.files(path = csvbase, pattern = csvpattern))

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