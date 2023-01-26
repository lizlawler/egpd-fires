library(scoringRules)
obs_n <- c(0, 1, 2)
sample_nm <- matrix(rnorm(3e4, mean = 2, sd = 3), nrow = 3)
crps_sample(obs_n, dat = sample_nm)


crps_edf <- function(y, dat, w = NULL) {
  if (is.null(w)) {
    c_1n <- 1 / length(dat)
    x <- sort(dat)
    a <- seq.int(0.5 * c_1n, 1 - 0.5 * c_1n, length.out = length(dat))
    f <- function(s) 2 * c_1n * sum(((s < x) - a) * (x - s))
  } else {
    if (!identical(length(dat), length(w)) || any(w < 0, na.rm = TRUE)) {
      return(rep(NaN, length(y)))
    }
    ord <- order(dat)
    x <- dat[ord]
    w <- w[ord]
    p <- cumsum(w)
    P <- p[length(p)]
    a <- (p - 0.5 * w) / P
    f <- function(s) 2 / P * sum(w * ((s < x) - a) * (x - s))
  }
  sapply(y, f)
}

crps_edf <- function(y, dat, w = NULL) {
  if (is.null(w)) {
    c_1n <- 1 / length(dat)
    x <- sort(dat)
    a <- seq.int(0.5 * c_1n, 1 - 0.5 * c_1n, length.out = length(dat))
    f <- function(s) 2 * c_1n * sum(((s < x) - a) * (x - s))
  } else {
    if (!identical(length(dat), length(w)) || any(w < 0, na.rm = TRUE)) {
      return(rep(NaN, length(y)))
    }
    ord <- order(dat)
    x <- dat[ord]
    w <- w[ord]
    p <- cumsum(w)
    P <- p[length(p)]
    a <- (p - 0.5 * w) / P
    f <- function(s) 2 / P * sum(w * ((s < x) - a) * (x - s))
  }
  sapply(y, f)
}


test_x <- sort(sample_nm)
test_c <- 1/length(sample_nm)
seq.int(0.5 * 1/length(sample_nm), 1 - 0.5 * (1/length(sample_nm)), length.out = length(sample_nm) )
order(sample_nm)
 
count58 <- stan_data_og$y_train_count[,58]
nonzerocount <- count58[count58 != 0]
hist(nonzerocount, breaks = seq(min(nonzerocount), max(nonzerocount) + 1, 0.5))


# ------ 
# creating functions for custom twCRPS

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


# for (i in 1:r) {
#   reg_kappa[, i] = X_train[i] * beta_kappa[, i] + phi_kappa[idx_train_er, i];
# }
# xi_init = cholesky_decompose(corr_xi)' * Z_xi;
#   xi_matrix = rep_matrix(xi_init', t_all);
# nu_init = cholesky_decompose(corr_nu)' * Z_nu;
#   nu_matrix = rep_matrix(nu_init', t_all);
# 
# kappa = exp(to_vector(reg_kappa))[ii_tb_all];
# nu = exp(to_vector(nu_matrix[idx_train_er,]))[ii_tb_all];
# xi = exp(to_vector(xi_matrix[idx_train_er,]))[ii_tb_all];
# sigma = nu ./ (1 + xi);
