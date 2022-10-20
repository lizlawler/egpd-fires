expit <- function(x) exp(x)/(1 + exp(x))

expit(-0.01)
library(rstan)
util <- new.env()
source('~/Desktop/research/stan_utility.R', local = util)

c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

c_light_trans <- c("#DCBCBC80")
c_dark_trans <- c("#8F272780")
c_green_trans <- c("#00FF0080")

diag_files <- sapply(1:3, function(x) paste("full-model/output/zip_fires_OG_small-step_2000_diag_", x, ".csv", sep=""))
unconstrained_fit <- read_stan_csv(diag_files)

partition <- util$partition_div(unconstrained_fit)
div_samples <- partition[[1]]
nondiv_samples <- partition[[2]]

plot(nondiv_samples$rho1_lambda, nondiv_samples$rho2_lambda,
     col=c_dark_trans, pch = 16, cex = 0.8,
     xlab = "rho1", ylab = "rho2")
points(div_samples$rho1_lambda, div_samples$rho2_lambda,
       col="green", pch=16, cex=0.8)

plot(nondiv_samples$rho1_pi, nondiv_samples$rho2_pi,
     col=c_dark_trans, pch = 16, cex = 0.8,
     xlab = "rho1", ylab = "rho2")
points(div_samples$rho1_pi, div_samples$rho2_pi,
       col="green", pch=16, cex=0.8)


plot(nondiv_samples$rho1_lambda, nondiv_samples$`beta_lambda[1,1]`,
     col=c_dark_trans, pch = 16, cex = 0.8,
     xlab = "rho1", ylab = "beta[1,1]")
points(div_samples$rho1_lambda, div_samples$`beta_lambda[1,1]`,
       col="green", pch=16, cex=0.8)


plot(nondiv_samples$rho1_pi, nondiv_samples$`beta_pi[1,1]`,
     col=c_dark_trans, pch = 16, cex = 0.8,
     xlab = "rho1", ylab = "beta[1,1]")
points(div_samples$rho1_pi, div_samples$`beta_pi[1,1]`,
       col="green", pch=16, cex=0.8)

plot(nondiv_samples$rho1_lambda,
     col=c_dark_trans, pch = 16, cex = 0.8,
     xlab = "Iteration", ylab = "rho1")
points(div_samples$rho1_lambda,
       col="green", pch=16, cex=0.8)
