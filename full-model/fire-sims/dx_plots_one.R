# library(cmdstanr)
# set_cmdstan_path(path = "/projects/eslawler@colostate.edu/software/anaconda/envs/lawler/bin/cmdstan") # this is only relevant to Alpine
# check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(rstan)
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

fit <- rstan::read_stan_csv(csvfiles)
list_of_draws <- extract(fit)
print(names(list_of_draws))
rho_fit <- as.array(fit, pars = 'rho')
beta_fit <- as.array(fit, pars = 'beta')
phi_fit <- as.array(fit, pars = 'phi')

MCMCtrace(rho_fit,
         ind = TRUE,
         open_pdf = FALSE,
         filename = paste0(plotbase, csvpattern, "_rho.pdf"))
MCMCtrace(beta_fit,
         ind = TRUE,
         open_pdf = FALSE,
         filename = paste0(plotbase, csvpattern, "_beta.pdf"))
MCMCtrace(phi_fit,
         ind = TRUE,
         open_pdf = FALSE,
         filename = paste0(plotbase, csvpattern, "_phi.pdf"))
