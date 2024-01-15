library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)

model <- cmdstan_model("models/joint/stan/joint_g1_sigma-ri.stan", compile = FALSE)
model$check_syntax(pedantic = TRUE)

