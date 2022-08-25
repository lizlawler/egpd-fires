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
  vector[r] phi_init[t];
  matrix[p, r] Z_kappa;
  matrix[p, r] beta_xi;
  matrix[p, r] beta_nu;
  real<lower = 0> tau;
  real<lower = 0, upper = 1> eta;
  real<lower = 0, upper = 1> bp_init;
  real<lower = 0, upper = 1> rho2;
  real<lower=0, upper = (1-rho2)> rho1;
}

transformed parameters {
  matrix[t, r] phi;
  matrix[r, t] reg_kappa;
  matrix[r, t] XB_nu;
  matrix[r, t] XB_xi;
  // matrix[r, t] XB_sigma;
  
  real<lower=0, upper = bp_init/2> bp = bp_init/2;
  
  matrix[r, r] corr_kappa = l3 + rho2 * l2 + rho1 * l1;
  matrix[p, p] cov_ar1 = equal + bp * bp_lin + bp^2 * bp_square + bp^3 * bp_cube + bp^4 * bp_quart;
  matrix[p, r] beta_kappa = cholesky_decompose(cov_ar1)' * Z_kappa * cholesky_decompose(corr_kappa);

  vector[t*r] kappa;
  vector[t*r] nu;
  vector[t*r] xi;
  vector[t*r] sigma;
  
  phi[1,] = tau * phi_init[1]';
  for (j in 2:t) {
    phi[j,] = eta * phi[j-1,] + tau * phi_init[j]';
  }
  
  for (i in 1:r) {
    reg_kappa[i, ] = (X[i] * beta_kappa[, i]/3 + phi[, i]/5)';
    XB_nu[i, ]  = (X[i] * beta_nu[, i])';
    XB_xi[i, ]  = (X[i] * beta_xi[, i]/4)';
  }

  kappa = to_vector(exp(reg_kappa));
  nu = to_vector(exp(XB_nu));
  xi = to_vector(exp(XB_xi));
  
  for (k in 1:(t*r)) {
    sigma[k] = nu[k] / (1 + xi[k]);
  }
}

model {
  // priors
  to_vector(Z_kappa) ~ normal(0, 1);
  to_vector(beta_xi) ~ normal(0, 1);
  to_vector(beta_nu) ~ normal(0, 1);
  
  bp_init ~ uniform(0, 1);
  
  rho1 ~ beta(1.5, 4);
  rho2 ~ beta(3, 4);
  
  // IAR prior
  eta ~ uniform(0,1);
  tau ~ gamma(0.001, 0.001);
  for (j in 1:t) {
    target += -.5 * dot_self(phi_init[j][node1] - phi_init[j][node2]);
    sum(phi_init[j]) ~ normal(0, 0.001*r);
  }
  
  // likelihood
  for (k in 1:(t*r)) {
    target += egpd_g1_lpdf(y[k] | sigma[k], xi[k], kappa[k]);
  }
}

