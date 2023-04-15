library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)


model <- cmdstan_model("full-model/fire-sims/burns/g2/stan/g2_xi-ri.stan", compile = FALSE)
model$check_syntax(pedantic = TRUE)

model <- cmdstan_model("full-model/fire-sims/counts/zip/stan/zip_pi-ri.stan", compile = FALSE)
model$check_syntax(pedantic = TRUE)


# # names correspond to the data block in the Stan program
# data_list <- list(N = 10, y = c(0,1,0,0,0,0,0,0,0,1))
# 
# fit <- mod$sample(
#   data = data_list, 
#   seed = 123, 
#   chains = 4, 
#   parallel_chains = 4,
#   refresh = 500 # print update every 500 iters
# )
# 
# fit <- readRDS(temp_rds_file)
# fit$summary()
# test <- rstan::read_stan_csv(fit$output_files())
# 
# 
# mod_pedantic <- cmdstan_model(stan_file_pedantic, pedantic = TRUE)
# stan_egpd <- cmdstan_model("full-model/fire-sims/burns/g1/stan/g1_all-reg.stan",
#                            pedantic = TRUE)
# egpd_opt <- cmdstan_model("full-model/fire-sims/burns/g1/stan/g1_all-reg.stan", 
#                           stanc_options = list("Oexperimental"))
# egpd_opt$print()
# stan_egpd$format(canonicalize = TRUE)
# 
# stanc("full-model/fire-sims/burns/g1/stan/g1_all-reg.stan")
# 
# 
# # model$format(canonicalize = TRUE, overwrite_file = TRUE)
# model$print()
# 
# stan_data <- readRDS("full-model/data/stan_data_sqrt.RDS")
# model <- cmdstan_model(stan_file="full-model/fire-sims/counts/zip/stan/zip_pi-ri.stan", compile = FALSE)
# model$check_syntax(pedantic = TRUE)
# 
# 
# model$format(canonicalize = TRUE, overwrite_file = TRUE)
# model$print()
# fit <- model$sample(data = stan_data,
#                     chains = 3,
#                     parallel_chains = 3,
#                     init = 0.01,
#                     output_dir = "full-model/fire-sims/burns/g1/")
# model$check_syntax()
# 
