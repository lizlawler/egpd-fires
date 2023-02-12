# choosing weight function
burn_df <- readRDS("full-model/data/burns.RDS")
plot(ecdf(sqrt(burn_df$BurnBndAc/1000)), xlim = c(0,10))
rate <- 1/mean(sqrt(burn_df$BurnBndAc/1000))
trunc_exp_cdf <- function(x, rate) {
  num = pexp(x, rate) - pexp(1.001, rate)
  denom = 1 - pexp(1.001, rate)
  return(num/denom)
}
curve(trunc_exp_cdf(x, rate), add = TRUE, col = 'green')


# ------ 
# creating functions for twCRPS
f1_cdf <- function(x, sigma = sigma, xi = xi, kappa = kappa) {
  (1 - (1 + xi * (x/sigma))^(-1/xi))^kappa
}

f1_pdf <- function(x, sigma = sigma, xi = xi, kappa = kappa) {
  lpdf <- log(kappa) - log(sigma) - (1/xi + 1) * log(1 + xi * (x/sigma)) + 
    (kappa-1) * log(1 - (1 + xi * (x/sigma))^(-1/xi))
  return(exp(lpdf))
}

g1_pdf <- function(x, sigma = sigma, xi = xi, kappa = kappa) {
  lower <- f1_cdf(1.001, sigma, xi, kappa)
  return(f1_pdf(x, sigma, xi, kappa)/(1-lower))
}

g1_cdf <- function(x, sigma = sigma, xi = xi, kappa = kappa) {
  lower <- f1_cdf(1.001, sigma, xi, kappa)
  return((f1_cdf(x, sigma, xi, kappa) - lower)/(1-lower))
}

g1_cdf_inv <- function(u, sigma = sigma, xi = xi, kappa = kappa) {
  lower <- f1_cdf(1.001, sigma, xi, kappa)
  u_adj <- u * (1-lower) + lower
  (sigma/xi) * ((1-u_adj^(1/kappa))^-xi - 1)
}

g1_rng <- function(n, sigma, xi, kappa) {
  u = runif(n)
  return(g1_cdf_inv(u, sigma, xi, kappa))
}

test <- g1_rng(100000, 1, 0.5, 2)
hist(test, freq = FALSE, breaks = seq(0, max(test)+1, 0.1), xlim = c(0, 10))
curve(g1_pdf(x, 1, 0.5, 2), xlim = c(0,10), add = TRUE)
plot(x=xgrid, y = g1_pdf(xgrid, 1, 0.5, 2), type = "l", add = TRUE)
xgrid <- seq(1.001,20,0.001)
pdf_real <- g1_pdf(xgrid, 1, 0.5, 2)
plot(y = pdf_real, x = xgrid, type = "l", xlim = c(0.5, 10), add = TRUE)

curve(g1_cdf(x, 4, 0.7, 10), xlim = c(0,50), ylim = c(0,1))
curve(pnorm(x, 0, 1), xlim = c(-5,5))

