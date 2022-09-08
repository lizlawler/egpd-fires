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
  int<lower = 1> t; // # of timepoints
  int<lower = 1> r; // # of regions
  int<lower=0> N_edges;
  int<lower=1, upper = r> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper = r> node2[N_edges];  // and node1[i] < node2[i]
  
  matrix[t, p] X[r]; // design matrix; 1-D array of size r with matrices t x p; this indexing is deprecated in version 2.32 and above
  vector[t*r] y; // response data
  
  // indicator matrices for ecoregions
  matrix[r, r] l3;
  matrix[r, r] l2;
  matrix[r, r] l1;

  // indicator matrices for AR(1) process on betas
  matrix[p, p] equal;
  matrix[p, p] bp_lin;
  matrix[p, p] bp_square;
  matrix[p, p] bp_cube;
  matrix[p, p] bp_quart;
}

parameters {
  vector[r] phi_init_delta[t];
  matrix[p, r] Z_delta;
  matrix[p, r] beta_xi;
  matrix[p, r] beta_nu;
  real<lower = 0> tau_init_delta;
  real<lower = 0, upper = 1> eta_delta;
  real<lower = 0, upper = 1> bp_init_delta;
  real<lower = 0, upper = 1> rho2_delta;
  real<lower=0, upper = (1-rho2_delta)> rho1_delta;
}

transformed parameters {
  matrix[t, r] phi_delta;
  matrix[r, t] reg_delta;
  matrix[r, t] reg_xi;
  matrix[r, t] reg_nu;

  real<lower=0, upper = bp_init_delta/2> bp_delta = bp_init_delta/2;
  real<lower=0, upper = tau_init_delta/2> tau_delta = tau_init_delta/2;

  matrix[r, r] corr_delta = l3 + rho2_delta * l2 + rho1_delta * l1;
  matrix[p, p] cov_ar1_delta = equal + bp_delta * bp_lin + bp_delta^2 * bp_square + bp_delta^3 * bp_cube + bp_delta^4 * bp_quart;
  matrix[p, r] beta_delta = cholesky_decompose(cov_ar1_delta)' * Z_delta * cholesky_decompose(corr_delta);

  vector[t*r] delta;
  vector[t*r] nu;
  vector[t*r] xi;
  vector[t*r] sigma;

  phi_delta[1,] = (1/tau_delta) * phi_init_delta[1]';
  for (j in 2:t) {
    phi_delta[j,] = eta_delta * phi_delta[j-1,] + (1/tau_delta) * phi_init_delta[j]';
  }
  
  for (i in 1:r) {
    reg_delta[i, ] = (X[i] * beta_delta[, i]/3 + phi_delta[, i]/15)';
    reg_nu[i, ] = (X[i] * beta_nu[, i]/4)';
    reg_xi[i, ] = (X[i] * beta_xi[, i]/10)';
  }

  delta = to_vector(exp(reg_delta));
  nu = to_vector(exp(reg_nu));
  xi = to_vector(exp(reg_xi));
  
  for (k in 1:(t*r)) {
    sigma[k] = nu[k] / (1 + xi[k]);
  }
}

model {
  // priors
  to_vector(Z_delta) ~ normal(0, 1);
  to_vector(beta_xi) ~ normal(0, 1);
  to_vector(beta_nu) ~ normal(0, 1);
  
  bp_init_delta ~ uniform(0, 1);

  rho1_delta ~ beta(1.5, 4);
  rho2_delta ~ beta(3, 4);
  
  // IAR prior
  eta_delta ~ beta(2,8);
  tau_init_delta ~ exponential(1);
  for (j in 1:t) {
    target += -.5 * dot_self(phi_init_delta[j][node1] - phi_init_delta[j][node2]);
    sum(phi_init_delta[j]) ~ normal(0, 0.001*r);
  }
  // 
  // likelihood
  for (k in 1:(t*r)) {
    target += egpd_g3_lpdf(y[k] | sigma[k], xi[k], delta[k]);
  }
}
