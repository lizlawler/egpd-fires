data {
  int<lower = 1> p; // # of parameters
  int<lower = 1> t; // # of timepoints
  int<lower = 1> r; // # of regions
  int<lower=0> N_edges;
  int<lower=1, upper = r> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper = r> node2[N_edges];  // and node1[i] < node2[i]
  
  matrix[t, p] X[r]; // design matrix; 1-D array of size r with matrices t x p; this indexing is deprecated in version 2.32 and above
  int<lower = 0> y[t*r]; // response data
  
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
  vector[r] phi_init_lambda[t];
  vector[r] phi_init_pi[t];
  matrix[p, r] Z_lambda;
  matrix[p, r] Z_pi;
  real<lower = 0> tau_init_lambda;
  real<lower = 0, upper = 1> eta_lambda;
  real<lower = 0, upper = 1> bp_init_lambda;
  real<lower = 0, upper = 1> rho2_lambda;
  real<lower=0, upper = (1-rho2_lambda)> rho1_lambda;
  real<lower = 0> tau_init_pi;
  real<lower = 0, upper = 1> eta_pi;
  real<lower = 0, upper = 1> bp_init_pi;
  real<lower = 0, upper = 1> rho2_pi;
  real<lower=0, upper = (1-rho2_pi)> rho1_pi;
  // real<lower = 0, upper = 1> pi_zip_init;
}

transformed parameters {
  matrix[t, r] phi_lambda;
  matrix[t, r] phi_pi;
  matrix[r, t] reg_lambda;
  matrix[r, t] reg_pi;
  vector<lower = 0> [t*r] lambda;
  vector<lower = 0, upper = 1> [t*r] pi_param;
  
  real<lower=0, upper = bp_init_lambda/2> bp_lambda = bp_init_lambda/2;
  real<lower=0, upper = tau_init_lambda/2> tau_lambda = tau_init_lambda/2;
  real<lower=0, upper = bp_init_pi/2> bp_pi = bp_init_pi/2;
  real<lower=0, upper = tau_init_pi/2> tau_pi = tau_init_pi/2;
  
  matrix[r, r] corr_lambda = l3 + rho2_lambda * l2 + rho1_lambda * l1;
  matrix[p, p] cov_ar1_lambda = equal + bp_lambda * bp_lin + bp_lambda^2 * bp_square + bp_lambda^3 * bp_cube + bp_lambda^4 * bp_quart;
  matrix[p, r] beta_lambda = cholesky_decompose(cov_ar1_lambda)' * Z_lambda * cholesky_decompose(corr_lambda);
  matrix[r, r] corr_pi = l3 + rho2_pi * l2 + rho1_pi * l1;
  matrix[p, p] cov_ar1_pi = equal + bp_pi * bp_lin + bp_pi^2 * bp_square + bp_pi^3 * bp_cube + bp_pi^4 * bp_quart;
  matrix[p, r] beta_pi = cholesky_decompose(cov_ar1_pi)' * Z_pi * cholesky_decompose(corr_pi);
  

  phi_lambda[1,] = (1/tau_lambda) * phi_init_lambda[1]';
  phi_pi[1,] = (1/tau_pi) * phi_init_pi[1]';
  for (j in 2:t) {
    phi_lambda[j,] = eta_lambda * phi_lambda[j-1,] + (1/tau_lambda) * phi_init_lambda[j]';
    phi_pi[j,] = eta_pi * phi_pi[j-1,] + (1/tau_pi) * phi_init_pi[j]';
  }
  
  for (i in 1:r) {
    reg_lambda[i, ] = (X[i] * beta_lambda[, i]/5 + phi_lambda[, i]/10)';
    reg_pi[i, ] = (X[i] * beta_pi[, i]/5 + phi_pi[, i]/10)';
  }

  lambda = to_vector(exp(reg_lambda));
  pi_param = to_vector(inv_logit(reg_pi));
}

model {
  // priors
  to_vector(Z_lambda) ~ normal(0, 1);
  to_vector(Z_pi) ~ normal(0, 1);

  bp_init_lambda ~ uniform(0, 1);
  bp_init_pi ~ uniform(0, 1);
  
  rho1_lambda ~ beta(1.5, 4);
  rho2_lambda ~ beta(3, 4);
  rho1_pi ~ beta(1.5, 4);
  rho2_pi ~ beta(3, 4);
  
  // IAR prior
  tau_init_lambda ~ exponential(1);
  tau_init_pi ~ exponential(1);
  eta_lambda ~ beta(2,8);
  eta_pi ~ beta(2,8);
  for (j in 1:t) {
    target += -.5 * dot_self(phi_init_lambda[j][node1] - phi_init_lambda[j][node2]);
    sum(phi_init_lambda[j]) ~ normal(0, 0.001*r);
    target += -.5 * dot_self(phi_init_pi[j][node1] - phi_init_pi[j][node2]);
    sum(phi_init_pi[j]) ~ normal(0, 0.001*r);
  }

  // likelihood
  for (i in 1:(t*r)) {
    if (y[i] == 0) {
      target += log_sum_exp(bernoulli_lpmf(1 | pi_param[i]),
                            bernoulli_lpmf(0 | pi_param[i])
                            + poisson_lpmf(y[i] | lambda[i]));
    } else {
      target += bernoulli_lpmf(0 | pi_param[i])
                + poisson_lpmf(y[i] | lambda[i]);
    }
  }
}
