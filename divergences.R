library(cmdstanr)
library(bayesplot)
library(ggplot2)
# setwd("~/research/egpd-fires")

files <- paste0("full-model/fire-sims/counts/zinb/csv-fits/",
                list.files(path = "full-model/fire-sims/counts/zinb/csv-fits/", 
                           pattern = "zinb_all-reg_27-Feb"))
zinb_fit <- as_cmdstan_fit(files[3], format = "draws_list")

lp_zinb <- log_posterior(zinb_fit)
np_zinb <- nuts_params(zinb_fit)
zinb_fit$diagnostic_summary()
posterior_zinb <- as.array(zinb_fit$draws())

unique(np_zinb$Parameter)
head(posterior_zinb)

rand_draws <- sample(1:346650, 10, replace = FALSE)
rand_posterior <- posterior_zinb[,,rand_draws]

mcmc_parcoord(
  rand_posterior,
  transform = function(x) {(x - mean(x)) / sd(x)},
  np = np_zinb
)
