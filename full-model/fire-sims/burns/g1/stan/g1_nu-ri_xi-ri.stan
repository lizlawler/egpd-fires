functions {
  real egpd_g1_lpdf(real y, real sigma, real xi, real kappa) {
    real lpdf;
    real cdf;
    lpdf = log(kappa) - log(sigma) - (1/xi + 1) * log(1 + xi * (y/sigma)) + 
    (kappa-1) * log(1 - (1 + xi * (y/sigma))^(-1/xi));
    cdf = (1 - (1 + xi * (1.001/sigma))^(-1/xi))^kappa;
    return lpdf/log(1 - cdf);
  }
  
  real matnormal_lpdf(matrix y, matrix cov, matrix corr) {
    real lpdf;
    real r;
    real p;
    r = rows(corr);
    p = rows(cov);
    lpdf = -(r*p/2) * log(2 * pi()) - (p/2)*log_determinant(corr) - (r/2)*log_determinant(cov) -
          0.5 * trace(mdivide_right_spd(mdivide_left_spd(corr, y'), cov) * y);
    return lpdf;
  }
}

data {
  int<lower = 1> r; // # of regions
  int<lower = 1> p; // # of parameters
  int<lower = 1> t_all; // # of timepoints in full dataset
  int<lower = 1> t_train;
  int<lower = 1> t_hold;
  
  // covariate data
  matrix[t_all, p] X_full[r]; // design matrix; 1-D array of size r with matrices t x p; this indexing is deprecated in version 2.32 and above
  matrix[t_train, p] X_train[r]; // design matrix; 1-D array of size r with matrices t x p; this indexing is deprecated in version 2.32 and above
  
  // training data
  int<lower = 1> N_tb_obs;
  int<lower = 1> N_tb_mis;
  int<lower = 1> N_tb_all;
  real<lower = 1> y_train_obs[N_tb_obs]; // original burn area
  int<lower = 1> ii_tb_obs[N_tb_obs];
  int<lower = 1, upper = N_tb_all> ii_tb_mis[N_tb_mis];
  int<lower = 1, upper = N_tb_all> ii_tb_all[N_tb_all]; // for broadcasting
  int<lower = 1> idx_train_er[t_train];

  // holdout data
  int<lower = 1> N_hold_obs;
  int<lower = 1> N_hold_all; // includes 'missing' and observed
  int<lower = 1> ii_hold_obs[N_hold_obs]; // vector of indices for holdout data timepoints
  int<lower = 1> ii_hold_all[N_hold_all]; // vector of indices for broadcasting to entire holdout dataset
  real<lower = 0> y_hold_obs[N_hold_obs]; // 
  int<lower = 1> idx_hold_er[t_hold];

  // neighbor information
  int<lower=0> n_edges;
  int<lower=1, upper = r> node1[n_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper = r> node2[n_edges];  // and node1[i] < node2[i]

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
  real<lower = 0> y_train_mis[N_tb_mis];
  vector[r] Z_xi;
  vector[r] Z_nu;
  vector[r] phi_init_kappa[t_all];
  matrix[p, r] beta_kappa;
  real<lower = 0> tau_init_kappa;
  real<lower = 0, upper = 1> eta_kappa;
  real<lower = 0, upper = 1> bp_init_kappa;
  real<lower = 0, upper = 1> rho2_kappa;
  real<lower=0, upper = (1-rho2_kappa)> rho1_kappa;
  real<lower = 0, upper = 1> rho2_nu;
  real<lower=0, upper = (1-rho2_nu)> rho1_nu;
  real<lower = 0, upper = 1> rho2_xi;
  real<lower=0, upper = (1-rho2_xi)> rho1_xi;
}

transformed parameters {
  real<lower = 0> y_train[N_tb_all];
  matrix[t_all, r] phi_kappa;
  matrix[t_train, r] reg_kappa;
  vector[r] nu_init;
  matrix[t_all, r] nu_matrix;
  vector[r] xi_init;
  matrix[t_all, r] xi_matrix;
  
  real<lower=0, upper = bp_init_kappa/2> bp_kappa = bp_init_kappa/2;
  real<lower=0, upper = tau_init_kappa/2> tau_kappa = tau_init_kappa/2;

  matrix[r, r] corr_kappa = l3 + rho2_kappa * l2 + rho1_kappa * l1;
  matrix[p, p] cov_ar1_kappa = equal + bp_kappa * bp_lin + bp_kappa^2 * bp_square + bp_kappa^3 * bp_cube + bp_kappa^4 * bp_quart;

  matrix[r, r] corr_nu = l3 + rho2_nu * l2 + rho1_nu * l1;
  matrix[r, r] corr_xi = l3 + rho2_xi * l2 + rho1_xi * l1;
  
  vector<lower = 0>[N_tb_all] kappa;
  vector<lower = 0>[N_tb_all] nu;
  vector<lower = 0>[N_tb_all] xi;
  vector<lower = 0>[N_tb_all] sigma;

  y_train[ii_tb_obs] = y_train_obs;
  y_train[ii_tb_mis] = y_train_mis;

  phi_kappa[1,] = (1/tau_kappa) * phi_init_kappa[1]';
  for (j in 2:t_all) {
    phi_kappa[j,] = eta_kappa * phi_kappa[j-1,] + (1/tau_kappa) * phi_init_kappa[j]';
  }
  
  for (i in 1:r) {
    reg_kappa[, i] = X_train[i] * beta_kappa[, i] + phi_kappa[idx_train_er, i];
  }
  xi_init = cholesky_decompose(corr_xi)' * Z_xi;
  xi_matrix = rep_matrix(xi_init', t_all);
  nu_init = cholesky_decompose(corr_nu)' * Z_nu;
  nu_matrix = rep_matrix(nu_init', t_all);

  kappa = exp(to_vector(reg_kappa))[ii_tb_all];
  nu = exp(to_vector(nu_matrix[idx_train_er,]))[ii_tb_all];
  xi = exp(to_vector(xi_matrix[idx_train_er,]))[ii_tb_all];
  sigma = nu ./ (1 + xi);
}

model {
  // priors
  bp_init_kappa ~ uniform(0, 1);
  to_vector(Z_nu) ~ normal(0, 1);
  to_vector(Z_xi) ~ normal(0, 1);
  
  rho1_kappa ~ beta(3, 4);
  rho2_kappa ~ beta(1.5, 4);
  rho1_nu ~ beta(3, 4);
  rho2_nu ~ beta(1.5, 4);
  rho1_xi ~ beta(3, 4);
  rho2_xi ~ beta(1.5, 4);
  
  target += matnormal_lpdf(beta_kappa | cov_ar1_kappa, corr_kappa);
  
  // IAR prior
  eta_kappa ~ beta(2,8);
  tau_init_kappa ~ exponential(1);
  for (j in 1:t_all) {
    target += -.5 * dot_self(phi_init_kappa[j][node1] - phi_init_kappa[j][node2]);
    sum(phi_init_kappa[j]) ~ normal(0, 0.001*r);
  }
  // 
  // likelihood
  for (k in 1:N_tb_all) {
    target += egpd_g1_lpdf(y_train[k] | sigma[k], xi[k], kappa[k]);
  }
}

generated quantities {
  matrix[t_all, r] reg_kappa_full;

  vector<lower = 0>[N_tb_obs] kappa_train;
  vector<lower = 0>[N_tb_obs] nu_train;
  vector<lower = 0>[N_tb_obs] xi_train;
  vector<lower = 0>[N_tb_obs] sigma_train;

  vector<lower = 0>[N_hold_obs] kappa_hold;
  vector<lower = 0>[N_hold_obs] nu_hold;
  vector<lower = 0>[N_hold_obs] xi_hold;
  vector<lower = 0>[N_hold_obs] sigma_hold;

  real holdout_loglik[N_hold_obs];
  real train_loglik[N_tb_obs];

  // expected values of EGPD components based on all timepoints, then cut to only be holdout timepoints
  for (i in 1:r) {
    reg_kappa_full[, i] = X_full[i] * beta_kappa[, i] + phi_kappa[, i];
  }

  kappa_train = exp(to_vector(reg_kappa_full))[ii_tb_all][ii_tb_obs];
  nu_train = exp(to_vector(nu_matrix))[ii_tb_all][ii_tb_obs];
  xi_train = exp(to_vector(xi_matrix))[ii_tb_all][ii_tb_obs];
  sigma_train = nu_train ./ (1 + xi_train);

  kappa_hold = exp(to_vector(reg_kappa_full))[ii_hold_all][ii_hold_obs];
  nu_hold = exp(to_vector(nu_matrix))[ii_hold_all][ii_hold_obs];
  xi_hold = exp(to_vector(xi_matrix))[ii_hold_all][ii_hold_obs];
  sigma_hold = nu_hold ./ (1 + xi_hold);

  if (max(y_train_obs) < 1200) { // condition determines if the data read in are the sqrt or original burn areas
    // training log-likelihood
    for (k in 1:N_tb_obs) {
      train_loglik[k] = egpd_g1_lpdf(y_train_obs[k] | sigma_train[k], xi_train[k], kappa_train[k]) + log(0.5) - log(y_train_obs[k]);
    }
    // holdout log-likelihood
    for (k in 1:N_hold_obs) {
      holdout_loglik[k] = egpd_g1_lpdf(y_hold_obs[k] | sigma_hold[k], xi_hold[k], kappa_hold[k]) + log(0.5) - log(y_hold_obs[k]);
    }
  }
  else {
    // training log-likelihood
    for (k in 1:N_tb_obs) {
      train_loglik[k] = egpd_g1_lpdf(y_train_obs[k] | sigma_train[k], xi_train[k], kappa_train[k]);
    }
    // holdout log-likelihood
    for (k in 1:N_hold_obs) {
      holdout_loglik[k] = egpd_g1_lpdf(y_hold_obs[k] | sigma_hold[k], xi_hold[k], kappa_hold[k]);
    }
  }
}
