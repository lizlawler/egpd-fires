## Naveau parametric families

g1_pdf <- function(x, sigma = sigma, xi = xi, kappa = kappa) {
  lpdf <- log(kappa) - log(sigma) - (1/xi + 1) * log(1 + xi * (x/sigma)) + 
    (kappa-1) * log(1 - (1 + xi * (x/sigma))^(-1/xi))
  return(exp(lpdf))
}

curve(g1_pdf(x, sigma = 1, xi = 0.5, kappa = 1), xlim = c(0, 7.5))
curve(g1_pdf(x, sigma = 1, xi = 0.5, kappa = 2), xlim = c(0, 7.5), add = TRUE, col = 'blue')
curve(g1_pdf(x, sigma = 1, xi = 0.5, kappa = 5), xlim = c(0, 7.5), add = TRUE, col = 'red')


g1_cdf <- function(x, sigma = sigma, xi = xi, kappa = kappa) {
  (1 - (1 + xi * (x/sigma))^(-1/xi))^kappa
}

g1_cdf_inv <- function(u, sigma = sigma, xi = xi, kappa = kappa) {
  (sigma/xi) * ((1-u^(1/kappa))^-xi - 1)
}

u1000 <- runif(1000)
hist(g1_cdf_inv(u1000, sigma = 1, xi = 0.5, kappa = ), freq = FALSE)

g1_pdf <- function(x, sigma = sigma, xi = xi, kappa = kappa) {
  lpdf <- log(kappa) - log(sigma) - (1/xi + 1) * log(1 + xi * (x/sigma)) + 
    (kappa-1) * log(1 - (1 + xi * (x/sigma))^(-1/xi))
  return(exp(lpdf))
}

curve(g1_pdf(x, sigma = 1, xi = 0.5, kappa = 1), xlim = c(0, 7.5), add = TRUE)
curve(g1_pdf(x, sigma = 1, xi = 0.5, kappa = 2), xlim = c(0, 7.5), add = TRUE, col = 'blue')
curve(g1_pdf(x, sigma = 1, xi = 0.5, kappa = 5), xlim = c(0, 7.5), add = TRUE, col = 'red')
