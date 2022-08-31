library(readr)
library(rstan)
library(MCMCvis)

toy_data <- read_rds('manuscript/scripts/toy_sim/data/toy_data_g2_base.rds')
rstan_options(auto_write = TRUE)

egpd_init <- stan_model('manuscript/scripts/toy_sim/toy_egpd_g2_base.stan')
egpd_fit <- sampling(egpd_init,
                     data = toy_data,
                     iter = 2000,
                     chains = 4)

MCMCtrace(egpd_fit, params = c("beta_kappa1", "beta_kappa2", "beta_nu", "beta_xi"), 
          ind = TRUE, 
          gvals = c(toy_data$truth$betas_kappa1, toy_data$truth$betas_kappa2, toy_data$truth$betas_nu, toy_data$truth$betas_xi))
# write_rds(egpd_fit_base, file = 'manuscript/scripts/toy_sim/egpd_fit_base.rds')

