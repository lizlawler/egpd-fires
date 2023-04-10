library(cmdstanr)
library(MCMCvis)
library(stringr)

all_csvs <- list.files(pattern = "*.csv")
model_name <- unique(str_remove(string = all_csvs, pattern = "_\\d{1}.csv"))
model_groups <- vector("list", length(model_name))
names(model_groups) <- model_name
for(i in 1:length(model_name)) {
  model_groups[[i]] <- c(all_csvs[grep(model_name[i], all_csvs)])
}

plotbase <- "~/research/egpd-fires/full-model/figures/g1/trace/"
trace_plots <- function(name_string, file_string) {
  name <- name_string
  fit <- cmdstanr::as_cmdstan_fit(file_string)
  fitmcmc <- cmdstanr::as_mcmc.list(fit)
  MCMCtrace(fitmcmc,
            params = 'rho',
            ind = TRUE,
            open_pdf = FALSE,
            filename = paste0(plotbase, name_string, "_rho.pdf"))
  MCMCtrace(fitmcmc,
            params = 'beta',
            ind = TRUE,
            open_pdf = FALSE,
            filename = paste0(plotbase, name_string, "_beta.pdf"))
  MCMCtrace(fitmcmc,
            params = 'phi',
            ind = TRUE,
            open_pdf = FALSE,
            filename = paste0(plotbase, name_string, "_phi.pdf"))
}

for(i in seq_along(model_groups)) {
  tryCatch(trace_plots(names(model_groups)[i], model_groups[[i]]), error=function(e) NULL)
}
  
