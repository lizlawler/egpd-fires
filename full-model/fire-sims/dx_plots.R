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
csv_base <- paste0(model, "_", suffix, "_", params, "_", delta)
csvfiles <- paste0(basepath, list.files(path = basepath, pattern = csv_base))

model_fit <- brms:::read_csv_as_stanfit(csvfiles)
MCMCtrace(model_fit, params = "rho",
          ind = TRUE,
          open_pdf = FALSE, 
          filename = paste0(plotbase, csv_base, "_rho.pdf"))
MCMCtrace(model_fit, params = "beta", 
          ind = TRUE,
          open_pdf = FALSE,
          filename = paste0(plotbase, csv_base, "_beta.pdf"))
MCMCtrace(model_fit, params = "phi",
          ind = TRUE,
          open_pdf = FALSE,
          filename = paste0(plotbase, csv_base, "_phi.pdf"))