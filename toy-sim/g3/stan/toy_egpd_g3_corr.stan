functions {
  real egpd_g3_lpdf(real y, real sigma, real xi, real delta) {
    real lpdf;
    lpdf = log(1 + delta) - log(delta * sigma) - (1/xi + 1) * log(1 + xi * (y/sigma)) + 
    log(1 - (1 + xi * (y/sigma))^(-delta/xi));
    return lpdf;
  }
}

data {
  int<lower = 1> p; // # of parameters
  int<lower = 1> n; // # of observations
  int<lower = 1> r; // # of regions
  matrix[n, p] X; // design matrix
  vector[n*r] y; // response data
  
  // indicator matrices for ecoregions
  matrix[r, r] l3;
  matrix[r, r] l2;
  matrix[r, r] l1;

  // diagonal matrix for testing rhos
  matrix[p, p] iden7;
  // // indicator matrices for AR(1) process on betas
  // matrix[p, p] equal;
  // matrix[p, p] bp_lin;
  // matrix[p, p] bp_square;
  // matrix[p, p] bp_cube;
  // matrix[p, p] bp_quart;
}

parameters {
  matrix[p, r] Z_delta;
  matrix[p, r] beta_nu;
  matrix[p, r] beta_xi;
  // real<lower = 0, upper = 1> bp_init;
  real<lower = 0, upper = 1> rho2;
  real<lower=0, upper = (1-rho2)> rho1;
}

transformed parameters {
  // real<lower=0, upper = bp_init/2> bp = bp_init/2;
  
  matrix[r, r] corr_delta = l3 + rho2 * l2 + rho1 * l1;
  // matrix[p, p] cov_ar1 = equal + bp * bp_lin + bp^2 * bp_square + bp^3 * bp_cube + bp^4 * bp_quart;
  matrix[p, r] beta_delta = cholesky_decompose(iden7)' * Z_delta * cholesky_decompose(corr_delta);
  // matrix[p, r] beta_delta = Z_delta;
  vector[n*r] delta;
  vector[n*r] nu;
  vector[n*r] xi;
  vector[n*r] sigma;
  
  delta = to_vector(exp(X * beta_delta));
  nu = to_vector(exp(X * beta_nu));
  xi = to_vector(exp(X * beta_xi));
  for (i in 1:(n*r)) {
    sigma[i] = nu[i] / (1 + xi[i]);
  }
}

model {
  // priors
  to_vector(Z_delta) ~ normal(0, 1);
  to_vector(beta_nu) ~ normal(0, 1);
  to_vector(beta_xi) ~ normal(0, 1);
  
  // bp_init ~ uniform(0, 1);
  
  rho2 ~ beta(3, 4);
  rho1 ~ beta(1.5, 4);

  // likelihood
  for (i in 1:(n*r)) {
    target += egpd_g3_lpdf(y[i] | sigma[i], xi[i], delta[i]);
  }
}


