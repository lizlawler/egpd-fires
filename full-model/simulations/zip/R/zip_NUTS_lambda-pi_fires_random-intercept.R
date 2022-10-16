library(readr)
library(rstan)
library(MCMCvis)
library(tidyverse)
library(splines)
library(spdep)
library(sf)
library(RColorBrewer)
library(patchwork)
library(classInt)
library(spatialreg)
options(mc.cores = parallel::detectCores())

# start time, for identification purposes later
st_time <- format(as.POSIXlt(Sys.time(), "America/Denver"), "%d-%b-%Y_%H%M")
stan_data <- readRDS(file = "./full-model/simulations/zip/data/stan_data_396_inc-zero-prob.RDS")

# run sampling of model with pi as random intercept (covariance = correlation matrix)
egpd_init <- stan_model('./full-model/simulations/zip/stan/zip_lambda-pi_fires_random-intercept.stan')
egpd_fit <- sampling(egpd_init, 
                     data = stan_data, 
                     iter = 1000,
                     chains = 3,
                     init_r = 0.01,
                     refresh = 50)

end_time <- format(as.POSIXlt(Sys.time(), "America/Denver"), "%H%M")

saveRDS(egpd_fit, 
        file = paste0("./full-model/simulations/zip/stan-fits/zip_lambda-pi_fires_396_random-intercept", 
                      st_time, "_", end_time, ".RDS"))

MCMCtrace(egpd_fit, params = c("rho1_lambda", "rho2_lambda", "rho1_pi", "rho2_pi"),
          ind = TRUE,
          open_pdf = FALSE,
          filename = paste0('./full-model/figures/zip/trace/zip_trace-lambda-pi_fires_396_random-intercept', 
                            st_time, "_", end_time, ".pdf"))




