functions {
  real egpd_g1_lpdf(real y, real sigma, real xi, real kappa) {
    real lpdf;
    lpdf = log(kappa) - log(sigma) - (1/xi + 1) * log(1 + xi * (y/sigma)) + 
    (kappa-1) * log(1 - (1 + xi * (y/sigma))^(-1/xi));
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
  vector[r] phi_init_nu[t];
  matrix[p, r] Z_nu;
  matrix[p, r] beta_kappa;
  matrix[p, r] beta_xi;
  real<lower = 0> tau_init_nu;
  real<lower = 0, upper = 1> eta_nu;
  real<lower = 0, upper = 1> bp_init_nu;
  real<lower = 0, upper = 1> rho2_nu;
  real<lower=0, upper = (1-rho2_nu)> rho1_nu;
}

transformed parameters {
  matrix[t, r] phi_nu;
  matrix[r, t] reg_kappa;
  matrix[r, t] reg_nu;
  matrix[r, t] reg_xi;

  real<lower=0, upper = bp_init_nu/2> bp_nu = bp_init_nu/2;
  real<lower=0, upper = tau_init_nu/2> tau_nu = tau_init_nu/2;

  matrix[r, r] corr_nu = l3 + rho2_nu * l2 + rho1_nu * l1;
  matrix[p, p] cov_ar1_nu = equal + bp_nu * bp_lin + bp_nu^2 * bp_square + bp_nu^3 * bp_cube + bp_nu^4 * bp_quart;
  matrix[p, r] beta_nu = cholesky_decompose(cov_ar1_nu)' * Z_nu * cholesky_decompose(corr_nu);

  vector[t*r] kappa;
  vector[t*r] nu;
  vector[t*r] xi;
  vector[t*r] sigma;

  phi_nu[1,] = (1/tau_nu) * phi_init_nu[1]';
  for (j in 2:t) {
    phi_nu[j,] = eta_nu * phi_nu[j-1,] + (1/tau_nu) * phi_init_nu[j]';
  }
  
  for (i in 1:r) {
    reg_kappa[i, ] = (X[i] * beta_kappa[, i])';
    reg_nu[i, ] = (X[i] * beta_nu[, i]/4 + phi_nu[, i]/10)';
    reg_xi[i, ] = (X[i] * beta_xi[, i]/8)';
  }

  kappa = to_vector(exp(reg_kappa));
  nu = to_vector(exp(reg_nu));
  xi = to_vector(exp(reg_xi));
  
  for (k in 1:(t*r)) {
    sigma[k] = nu[k] / (1 + xi[k]);
  }
}

model {
  // priors
  to_vector(Z_nu) ~ normal(0, 1);
  to_vector(beta_kappa) ~ normal(0, 1);
  to_vector(beta_xi) ~ normal(0, 1);
  
  bp_init_nu ~ uniform(0, 1);

  rho1_nu ~ beta(1.5, 4);
  rho2_nu ~ beta(3, 4);
  
  // IAR prior
  eta_nu ~ beta(2,8);
  tau_init_nu ~ exponential(1);
  for (j in 1:t) {
    target += -.5 * dot_self(phi_init_nu[j][node1] - phi_init_nu[j][node2]);
    sum(phi_init_nu[j]) ~ normal(0, 0.001*r);
  }
  // 
  // likelihood
  for (k in 1:(t*r)) {
    target += egpd_g1_lpdf(y[k] | sigma[k], xi[k], kappa[k]);
  }
}
