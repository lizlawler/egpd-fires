library(readr)
library(rstan)
egpd_data <- read_rds('data/processed/egpd_data.rds')

egpd_init <- stan_model('R/egpd_power.stan')
egpd_fit <- sampling(egpd_init,
              data = egpd_data,
              iter = 100,
              chains = 1)
write_rds(zip_fit, path = 'zip_fit.rds')
