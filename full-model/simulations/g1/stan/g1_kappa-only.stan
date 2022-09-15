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
  int<lower = 1> t; // # of total timepoints
  int<lower = 1> t_train; // # of training timepoints
  int<lower = 1> t_hold; // # of holdout timepoints
  int<lower = 1> r; // # of regions
  int<lower=0> N_edges;
  int<lower=1, upper = r> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper = r> node2[N_edges];  // and node1[i] < node2[i]
  
  // full data
  matrix[t, p] X[r]; // design matrix; 1-D array of size r with matrices t x p; this indexing is deprecated in version 2.32 and above
  // vector[t*r] y; // response data
  
  // training data
  int<lower = 1, upper = t> train_idx[t_train];
  matrix[t_train, p] X_train[r]; // design matrix; 1-D array of size r with matrices t_train x p; this indexing is deprecated in version 2.32 and above
  vector[t_train*r] y_train; // response data
  
  // holdout data
  int<lower = 1, upper = t>  hold_idx[t_hold];
  matrix[t_hold, p] X_hold[r]; // design matrix; 1-D array of size r with matrices t_hold x p; this indexing is deprecated in version 2.32 and above
  vector[t_hold*r] y_hold; // response data
  
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
  vector[r] phi_init_kappa[t];
  matrix[p, r] Z_kappa;
  matrix[p, r] beta_xi;
  matrix[p, r] beta_nu;
  real<lower = 0> tau_init_kappa;
  real<lower = 0, upper = 1> eta_kappa;
  real<lower = 0, upper = 1> bp_init_kappa;
  real<lower = 0, upper = 1> rho2_kappa;
  real<lower=0, upper = (1-rho2_kappa)> rho1_kappa;
}

transformed parameters {
  matrix[t, r] phi_kappa;
  matrix[r, t_train] reg_kappa;
  matrix[r, t_train] reg_xi;
  matrix[r, t_train] reg_nu;

  real<lower=0, upper = bp_init_kappa/2> bp_kappa = bp_init_kappa/2;
  real<lower=0, upper = tau_init_kappa/2> tau_kappa = tau_init_kappa/2;

  matrix[r, r] corr_kappa = l3 + rho2_kappa * l2 + rho1_kappa * l1;
  matrix[p, p] cov_ar1_kappa = equal + bp_kappa * bp_lin + bp_kappa^2 * bp_square + bp_kappa^3 * bp_cube + bp_kappa^4 * bp_quart;
  matrix[p, r] beta_kappa = cholesky_decompose(cov_ar1_kappa)' * Z_kappa * cholesky_decompose(corr_kappa);

  vector[t_train*r] kappa;
  vector[t_train*r] nu;
  vector[t_train*r] xi;
  vector[t_train*r] sigma;

  phi_kappa[1,] = (1/tau_kappa) * phi_init_kappa[1]';
  for (j in 2:t) {
    phi_kappa[j,] = eta_kappa * phi_kappa[j-1,] + (1/tau_kappa) * phi_init_kappa[j]';
  }
  
  for (i in 1:r) {
    reg_kappa[i, ] = (X_train[i] * beta_kappa[, i]/4 + phi_kappa[, i][train_idx]/15)';
    reg_nu[i, ] = (X_train[i] * beta_nu[, i]/4)';
    reg_xi[i, ] = (X_train[i] * beta_xi[, i]/8)';
  }

  kappa = to_vector(exp(reg_kappa));
  nu = to_vector(exp(reg_nu));
  xi = to_vector(exp(reg_xi));
  
  for (k in 1:(t_train*r)) {
    sigma[k] = nu[k] / (1 + xi[k]);
  }
}

model {
  // priors
  to_vector(Z_kappa) ~ normal(0, 1);
  to_vector(beta_xi) ~ normal(0, 1);
  to_vector(beta_nu) ~ normal(0, 1);
  
  bp_init_kappa ~ uniform(0, 1);

  rho1_kappa ~ beta(1.5, 4);
  rho2_kappa ~ beta(3, 4);
  
  // IAR prior
  eta_kappa ~ beta(2,8);
  tau_init_kappa ~ exponential(1);
  for (j in 1:t) {
    target += -.5 * dot_self(phi_init_kappa[j][node1] - phi_init_kappa[j][node2]);
    sum(phi_init_kappa[j]) ~ normal(0, 0.001*r);
  }
  // 
  // likelihood
  for (k in 1:(t_train*r)) {
    target += egpd_g1_lpdf(y_train[k] | sigma[k], xi[k], kappa[k]);
  }
}

generated quantities {
  vector[t_hold*r] holdout_loglik;
  vector[t_train*r] train_loglik;
  
  matrix[r, t] reg_kappa_full;  
  matrix[r, t] reg_nu_full;  
  matrix[r, t] reg_xi_full;  
  
  vector[t_hold*r] kappa_hold;
  vector[t_hold*r] nu_hold;
  vector[t_hold*r] xi_hold;
  vector[t_hold*r] sigma_hold;
  
  // expected values of parameters based on all timepoints, then cut to only be holdout parameters
  for (i in 1:r) {
    reg_kappa_full[i, ] = (X[i] * beta_kappa[, i]/4 + phi_kappa[, i]/15)';
    reg_nu_full[i, ] = (X[i] * beta_nu[, i]/4)';
    reg_xi_full[i, ] = (X[i] * beta_xi[, i]/8)';
  }

  kappa_hold = to_vector(exp(reg_kappa_full[, hold_idx]));
  nu_hold = to_vector(exp(reg_nu_full[, hold_idx]));
  xi_hold = to_vector(exp(reg_xi_full[, hold_idx]));
  
  for (k in 1:(t_hold*r)) {
    sigma_hold[k] = nu_hold[k] / (1 + xi_hold[k]);
  }

  // training log likelihoods
  for (k in 1:(t_train*r)) {
    train_loglik[k] = egpd_g1_lpdf(y_train[k] | sigma[k], xi[k], kappa[k]);
  }
  
  // holdout log likelihoods
  for (k in 1:(t_hold*r)) {
    holdout_loglik[k] = egpd_g1_lpdf(y_hold[k] | sigma_hold[k], xi_hold[k], kappa_hold[k]);
  }
}
