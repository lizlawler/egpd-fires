data {
  int<lower = 1> p; // # of parameters
  int<lower = 1> t; // # of timepoints
  int<lower = 1> r; // # of regions
  int<lower=0> N_edges;
  int<lower=1, upper = r> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper = r> node2[N_edges];  // and node1[i] < node2[i]
  
  matrix[t, p] X[r]; // design matrix; 1-D array of size r with matrices t x p; this indexing is deprecated in version 2.32 and above
  int<lower = 0> y[t*r]; // response data
  real<lower = 0, upper = 1> pi_val;
  
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
  matrix[p, r] Z_lambda;
  real<lower = 0> tau_init;
  real<lower = 0, upper = 1> eta;
  real<lower = 0, upper = 1> bp_init;
  real<lower = 0, upper = 1> rho2;
  real<lower=0, upper = (1-rho2)> rho1;
  // real<lower = 0, upper = 1> pi_zip_init;
}

transformed parameters {
  matrix[t, r] phi;
  matrix[r, t] reg_lambda;
  vector<lower = 0> [t*r] lambda;
  
  real<lower=0, upper = bp_init/2> bp = bp_init/2;
  real<lower=0, upper = tau_init/2> tau = tau_init/2;
  // real<lower=0, upper = pi_zip_init/2> pi_zip = pi_zip_init/2;
  
  matrix[r, r] corr_lambda = l3 + rho2 * l2 + rho1 * l1;
  matrix[p, p] cov_ar1 = equal + bp * bp_lin + bp^2 * bp_square + bp^3 * bp_cube + bp^4 * bp_quart;
  matrix[p, r] beta_lambda = cholesky_decompose(cov_ar1)' * Z_lambda * cholesky_decompose(corr_lambda);

  phi[1,] = (1/tau) * phi_init[1]';
  for (j in 2:t) {
    phi[j,] = eta * phi[j-1,] + (1/tau) * phi_init[j]';
  }
  
  for (i in 1:r) {
    reg_lambda[i, ] = (X[i] * beta_lambda[, i]/5 + phi[, i]/10)';
    // reg_lambda[i, ] = (X[i] * beta_lambda[, i])';
  }

  lambda = to_vector(exp(reg_lambda));
}

model {
  // priors
  to_vector(Z_lambda) ~ normal(0, 1);
  // pi_zip_init ~ beta(2, 8);
  bp_init ~ uniform(0, 1);
  
  rho1 ~ beta(1.5, 4);
  rho2 ~ beta(3, 4);
  
  // IAR prior
  tau_init ~ exponential(1);
  eta ~ beta(2,8);
  for (j in 1:t) {
    target += -.5 * dot_self(phi_init[j][node1] - phi_init[j][node2]);
    sum(phi_init[j]) ~ normal(0, 0.001*r);
  }

  // likelihood
  for (i in 1:(t*r)) {
    if (y[i] == 0) {
      target += log_sum_exp(bernoulli_lpmf(1 | pi_val),
                            bernoulli_lpmf(0 | pi_val)
                            + poisson_lpmf(y[i] | lambda[i]));
    } else {
      target += bernoulli_lpmf(0 | pi_val)
                + poisson_lpmf(y[i] | lambda[i]);
    }
  }
}
