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
  vector[r] phi_init_xi[t];
  matrix[p, r] Z_xi;
  matrix[p, r] beta_nu;
  matrix[p, r] beta_delta;
  real<lower = 0> tau_init_xi;
  real<lower = 0, upper = 1> eta_xi;
  real<lower = 0, upper = 1> bp_init_xi;
  real<lower = 0, upper = 1> rho2_xi;
  real<lower=0, upper = (1-rho2_xi)> rho1_xi;
}

transformed parameters {
  matrix[t, r] phi_xi;
  matrix[r, t] reg_delta;
  matrix[r, t] reg_xi;
  matrix[r, t] reg_nu;

  real<lower=0, upper = bp_init_xi/2> bp_xi = bp_init_xi/2;
  real<lower=0, upper = tau_init_xi/2> tau_xi = tau_init_xi/2;

  matrix[r, r] corr_xi = l3 + rho2_xi * l2 + rho1_xi * l1;
  matrix[p, p] cov_ar1_xi = equal + bp_xi * bp_lin + bp_xi^2 * bp_square + bp_xi^3 * bp_cube + bp_xi^4 * bp_quart;
  matrix[p, r] beta_xi = cholesky_decompose(cov_ar1_xi)' * Z_xi * cholesky_decompose(corr_xi);

  vector[t*r] delta;
  vector[t*r] nu;
  vector[t*r] xi;
  vector[t*r] sigma;

  phi_xi[1,] = (1/tau_xi) * phi_init_xi[1]';
  for (j in 2:t) {
    phi_xi[j,] = eta_xi * phi_xi[j-1,] + (1/tau_xi) * phi_init_xi[j]';
  }
  
  for (i in 1:r) {
    reg_delta[i, ] = (X[i] * beta_delta[, i] )';
    reg_nu[i, ] = (X[i] * beta_nu[, i]/5 )';
    reg_xi[i, ] = (X[i] * beta_xi[, i]/10 + phi_xi[, i]/15)';
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
  to_vector(Z_xi) ~ normal(0, 1);
  to_vector(beta_nu) ~ normal(0, 1);
  to_vector(beta_delta) ~ normal(0, 1);
  
  bp_init_xi ~ uniform(0, 1);

  rho1_xi ~ beta(1.5, 4);
  rho2_xi ~ beta(3, 4);
  
  // IAR prior
  eta_xi ~ beta(2,8);
  tau_init_xi ~ exponential(1);
  for (j in 1:t) {
    target += -.5 * dot_self(phi_init_xi[j][node1] - phi_init_xi[j][node2]);
    sum(phi_init_xi[j]) ~ normal(0, 0.001*r);
  }
  // 
  // likelihood
  for (k in 1:(t*r)) {
    target += egpd_g3_lpdf(y[k] | sigma[k], xi[k], delta[k]);
  }
}
