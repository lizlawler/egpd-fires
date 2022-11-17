functions {
  real egpd_g3_lpdf(real y, real sigma, real xi, real delta) {
    real lpdf;
    real cdf;
    lpdf = log(1 + delta) - log(delta * sigma) - (1/xi + 1) * log(1 + xi * (y/sigma)) + 
    log(1 - (1 + xi * (y/sigma))^(-delta/xi));
    cdf = 1 - beta_cdf((1 + xi * (1.001/sigma))^(-delta/xi), 1/delta, 2);
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
  vector[r] phi_init_delta[t_all];
  vector[r] phi_init_nu[t_all];
  matrix[p, r] beta_delta;
  matrix[p, r] beta_nu;
  real<lower = 0> tau_init_delta;
  real<lower = 0> tau_init_nu;
  real<lower = 0, upper = 1> eta_delta;
  real<lower = 0, upper = 1> eta_nu;
  real<lower = 0, upper = 1> bp_init_delta;
  real<lower = 0, upper = 1> bp_init_nu;
  real<lower = 0, upper = 1> rho2_delta;
  real<lower=0, upper = (1-rho2_delta)> rho1_delta;
  real<lower = 0, upper = 1> rho2_nu;
  real<lower=0, upper = (1-rho2_nu)> rho1_nu;
  real<lower = 0, upper = 1> rho2_xi;
  real<lower=0, upper = (1-rho2_xi)> rho1_xi;
}

transformed parameters {
  real<lower = 0> y_train[N_tb_all];
  matrix[t_all, r] phi_delta;
  matrix[t_all, r] phi_nu;
  matrix[t_train, r] reg_delta;
  matrix[t_train, r] reg_nu;
  vector[r] xi_init;
  matrix[t_all, r] xi_matrix;

  real<lower=0, upper = bp_init_delta/2> bp_delta = bp_init_delta/2;
  real<lower=0, upper = bp_init_nu/2> bp_nu = bp_init_nu/2;
  real<lower=0, upper = tau_init_delta/2> tau_delta = tau_init_delta/2;
  real<lower=0, upper = tau_init_nu/2> tau_nu = tau_init_nu/2;

  matrix[r, r] corr_delta = l3 + rho2_delta * l2 + rho1_delta * l1;
  matrix[p, p] cov_ar1_delta = equal + bp_delta * bp_lin + bp_delta^2 * bp_square + bp_delta^3 * bp_cube + bp_delta^4 * bp_quart;

  matrix[r, r] corr_nu = l3 + rho2_nu * l2 + rho1_nu * l1;
  matrix[p, p] cov_ar1_nu = equal + bp_nu * bp_lin + bp_nu^2 * bp_square + bp_nu^3 * bp_cube + bp_nu^4 * bp_quart;

  matrix[r, r] corr_xi = l3 + rho2_xi * l2 + rho1_xi * l1;
  
  vector<lower = 0>[N_tb_all] delta;
  vector<lower = 0>[N_tb_all] nu;
  vector<lower = 0>[N_tb_all] xi;
  vector<lower = 0>[N_tb_all] sigma;

  y_train[ii_tb_obs] = y_train_obs;
  y_train[ii_tb_mis] = y_train_mis;

  phi_delta[1,] = (1/tau_delta) * phi_init_delta[1]';
  phi_nu[1,] = (1/tau_nu) * phi_init_nu[1]';
  for (j in 2:t_all) {
    phi_delta[j,] = eta_delta * phi_delta[j-1,] + (1/tau_delta) * phi_init_delta[j]';
    phi_nu[j,] = eta_nu * phi_nu[j-1,] + (1/tau_nu) * phi_init_nu[j]';
  }
  
  for (i in 1:r) {
    reg_delta[, i] = X_train[i] * beta_delta[, i] + phi_delta[idx_train_er, i];
    reg_nu[, i] = X_train[i] * beta_nu[, i] + phi_nu[idx_train_er, i];
  }

  xi_init = cholesky_decompose(corr_xi)' * Z_xi;
  xi_matrix = rep_matrix(xi_init', t_all);

  delta = exp(to_vector(reg_delta))[ii_tb_all];
  nu = exp(to_vector(reg_nu))[ii_tb_all];
  xi = exp(to_vector(xi_matrix[idx_train_er,]))[ii_tb_all];
  sigma = nu ./ (1 + xi);
}

model {
  // priors
  bp_init_delta ~ uniform(0, 1);
  bp_init_nu ~ uniform(0, 1);
  Z_xi ~ normal(0, 1);
  
  rho1_delta ~ beta(3, 4);
  rho2_delta ~ beta(1.5, 4);
  rho1_nu ~ beta(3, 4);
  rho2_nu ~ beta(1.5, 4);
  rho1_xi ~ beta(3, 4);
  rho2_xi ~ beta(1.5, 4);
  
  target += matnormal_lpdf(beta_delta | cov_ar1_delta, corr_delta);
  target += matnormal_lpdf(beta_nu | cov_ar1_nu, corr_nu);
  
  // IAR prior
  eta_delta ~ beta(2,8);
  eta_nu ~ beta(2,8);
  tau_init_delta ~ exponential(1);
  tau_init_nu ~ exponential(1);
  for (j in 1:t_all) {
    target += -.5 * dot_self(phi_init_delta[j][node1] - phi_init_delta[j][node2]);
    sum(phi_init_delta[j]) ~ normal(0, 0.001*r);
    target += -.5 * dot_self(phi_init_nu[j][node1] - phi_init_nu[j][node2]);
    sum(phi_init_nu[j]) ~ normal(0, 0.001*r);
  }
  // 
  // likelihood
  for (k in 1:N_tb_all) {
    target += egpd_g3_lpdf(y_train[k] | sigma[k], xi[k], delta[k]);
  }
}

generated quantities {
  matrix[t_all, r] reg_delta_full;
  matrix[t_all, r] reg_nu_full;

  vector<lower = 0>[N_tb_obs] delta_train;
  vector<lower = 0>[N_tb_obs] nu_train;
  vector<lower = 0>[N_tb_obs] xi_train;
  vector<lower = 0>[N_tb_obs] sigma_train;

  vector<lower = 0>[N_hold_obs] delta_hold;
  vector<lower = 0>[N_hold_obs] nu_hold;
  vector<lower = 0>[N_hold_obs] xi_hold;
  vector<lower = 0>[N_hold_obs] sigma_hold;

  real holdout_loglik[N_hold_obs];
  real train_loglik[N_tb_obs];

  // expected values of EGPD components based on all timepoints, then cut to only be holdout timepoints
  for (i in 1:r) {
    reg_delta_full[, i] = X_full[i] * beta_delta[, i] + phi_delta[, i];
    reg_nu_full[, i] = X_full[i] * beta_nu[, i] + phi_nu[, i];
  }

  delta_train = exp(to_vector(reg_delta_full))[ii_tb_all][ii_tb_obs];
  nu_train = exp(to_vector(reg_nu_full))[ii_tb_all][ii_tb_obs];
  xi_train = exp(to_vector(xi_matrix))[ii_tb_all][ii_tb_obs];
  sigma_train = nu_train ./ (1 + xi_train);

  delta_hold = exp(to_vector(reg_delta_full))[ii_hold_all][ii_hold_obs];
  nu_hold = exp(to_vector(reg_nu_full))[ii_hold_all][ii_hold_obs];
  xi_hold = exp(to_vector(xi_matrix))[ii_hold_all][ii_hold_obs];
  sigma_hold = nu_hold ./ (1 + xi_hold);

  if (max(y_train_obs) < 1200) { // condition determines if the data read in are the sqrt or original burn areas
    // training log-likelihood
    for (k in 1:N_tb_obs) {
      train_loglik[k] = egpd_g3_lpdf(y_train_obs[k] | sigma_train[k], xi_train[k], delta_train[k]) + log(0.5) - log(y_train_obs[k]);
    }
    // holdout log-likelihood
    for (k in 1:N_hold_obs) {
      holdout_loglik[k] = egpd_g3_lpdf(y_hold_obs[k] | sigma_hold[k], xi_hold[k], delta_hold[k]) + log(0.5) - log(y_hold_obs[k]);
    }
  }
  else {
    // training log-likelihood
    for (k in 1:N_tb_obs) {
      train_loglik[k] = egpd_g3_lpdf(y_train_obs[k] | sigma_train[k], xi_train[k], delta_train[k]);
    }
    // holdout log-likelihood
    for (k in 1:N_hold_obs) {
      holdout_loglik[k] = egpd_g3_lpdf(y_hold_obs[k] | sigma_hold[k], xi_hold[k], delta_hold[k]);
    }
  }
}
