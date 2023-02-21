functions {
  real egpd_g1_lpdf(real y, real sigma, real xi, real kappa) {
    real lpdf;
    real cst;
    lpdf = log(kappa) - log(sigma) - (1/xi + 1) * log(1 + xi * (y/sigma)) + 
    (kappa-1) * log(1 - (1 + xi * (y/sigma))^(-1/xi));
    cst = (1 - (1 + xi * (1.001/sigma))^(-1/xi))^kappa;
    return lpdf - log1m(cst);
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
  real<lower = 1> y_hold_obs[N_hold_obs]; // 
  int<lower = 1> idx_hold_er[t_hold];

  // neighbor information
  int<lower=0> n_edges;
  int<lower=1, upper = r> node1[n_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper = r> node2[n_edges];  // and node1[i] < node2[i]

  // indicator matrices for ecoregions
  matrix[r, r] l3;
  matrix[r, r] l2;
  matrix[r, r] l1;

  // indicator matrices for AR(1) penalization on spline coefficients of betas
  matrix[p, p] equal;
  matrix[p, p] bp_lin;
  matrix[p, p] bp_square;
  matrix[p, p] bp_cube;
  matrix[p, p] bp_quart;
}

transformed data {
  int S = 1; // # of parameters with regression (ranges from 1 to 3)
  int C = 3; // # of parameters with correlation (always 3)
}

parameters {
  real<lower = 1> y_train_mis[N_tb_mis];
  vector[r] Z[2];
  row_vector[r] phi_init[t_all, S];
  matrix[p, r] beta[S];
  real<lower = 0> tau_init[S];
  real<lower = 0, upper = 1> eta[S];
  real<lower = 0, upper = 1> bp_init[S];
  real<lower = 0, upper = 1> rho[2, C];
}

transformed parameters {
  real<lower = 1> y_train[N_tb_all];
  matrix[t_all, r] phi[S];
  matrix[t_train, r] reg[S];
  matrix[p, p] cov_ar1[S];
  real<lower = 0, upper = 1> bp[S];
  real<lower = 0> tau[S];
  matrix[r, r] corr[C]; // 1 = kappa, 2 = nu, 3 = xi
  vector[r] ri_init[2]; // random intercept vector; 1= nu, 2=xi
  matrix[t_all, r] ri_matrix[2]; // broadcast ri_init to full matrix
  
  vector<lower = 0>[N_tb_all] kappa;
  vector<lower = 0>[N_tb_all] sigma;
  vector<lower = 0>[N_tb_all] xi;

  y_train[ii_tb_obs] = y_train_obs;
  y_train[ii_tb_mis] = y_train_mis;

  for (i in 1:C) {
    corr[i] = l3 + rho[2, i] * l2 + rho[1, i] * l1;
  }

  for (i in 1:2) {
    ri_init[i] = cholesky_decompose(corr[i+1])' * Z[i];
    ri_matrix[i] = rep_matrix(ri_init[i]', t_all);
  }

  for (i in 1:S) {
    bp[i] = bp_init[i]/2;
    tau[i] = tau_init[i]/2;
    cov_ar1[i] = equal + bp[i] * bp_lin + bp[i]^2 * bp_square + bp[i]^3 * bp_cube + bp[i]^4 * bp_quart;
    
    // ICAR variables
    phi[i][1, ] = 1/tau[i] * phi_init[1, i];
    for (j in 2:t_all) {
      phi[i][j, ] = eta[i] * phi[i][j-1, ] + 1/tau[i] * phi_init[j, i];
    }
    
    // regression for kappa
    for (k in 1:r) {
      reg[i][, k] = X_train[k] * beta[i][, k] + phi[i][idx_train_er, k];
    }
  }

  kappa = exp(to_vector(reg[1]))[ii_tb_all];
  sigma = exp(to_vector(ri_matrix[1][idx_train_er,]))[ii_tb_all];
  xi = exp(to_vector(ri_matrix[2][idx_train_er,]))[ii_tb_all];
}

model {
  Z[1] ~ normal(0, 1);
  Z[2] ~ normal(0, 1);
  // priors on rhos and AR(1) penalization of splines
  to_vector(bp_init) ~ uniform(0, 1);
  to_vector(rho[1,]) ~ beta(3, 4); // prior on rho1 for kappa, nu, and xi
  // to_vector(rho[2,]) ~ beta(1.5, 4); // prior on rho2 for kappa, nu, and xi
  
  // priors scaling constants in ICAR
  to_vector(eta) ~ beta(2,8);
  to_vector(tau_init) ~ exponential(1);
  
  for (i in 1:C) {
    // soft constraint for sum of rhos within an individual param to be <= 1 (ie rho1kappa + rho2kappa <= 1)
    sum(rho[,i]) ~ uniform(0,1); 
  }
  
  for (i in 1:S) {
    // MVN prior on betas
    target += matnormal_lpdf(beta[i] | cov_ar1[i], corr[i]);
    // ICAR prior
    for (j in 1:t_all) {
      target += -.5 * dot_self(phi_init[j, i][node1] - phi_init[j, i][node2]);
      sum(phi_init[j, i]) ~ normal(0, 0.001*r);
    }
  }

  // likelihood
  for (k in 1:N_tb_all) {
    target += egpd_g1_lpdf(y_train[k] | sigma[k], xi[k], kappa[k]);
  }
}

generated quantities {
  matrix[t_all, r] reg_full[S];

  vector<lower = 0>[N_tb_obs] kappa_train;
  vector<lower = 0>[N_tb_obs] sigma_train;
  vector<lower = 0>[N_tb_obs] xi_train;

  vector<lower = 0>[N_hold_obs] kappa_hold;
  vector<lower = 0>[N_hold_obs] sigma_hold;
  vector<lower = 0>[N_hold_obs] xi_hold;

  real holdout_loglik[N_hold_obs];
  real train_loglik[N_tb_obs];

  for (i in 1:S) {
    for (k in 1:r) {
      reg_full[i][,k] = X_full[k] * beta[i][, k] + phi[i][, k];
    }
  }

  kappa_train = exp(to_vector(reg_full[1]))[ii_tb_all][ii_tb_obs];
  sigma_train = exp(to_vector(ri_matrix[1]))[ii_tb_all][ii_tb_obs];
  xi_train = exp(to_vector(ri_matrix[2]))[ii_tb_all][ii_tb_obs];

  kappa_hold = exp(to_vector(reg_full[1]))[ii_hold_all][ii_hold_obs];
  sigma_hold = exp(to_vector(ri_matrix[1]))[ii_hold_all][ii_hold_obs];
  xi_hold = exp(to_vector(ri_matrix[2]))[ii_hold_all][ii_hold_obs];

  if (max(y_train_obs) < 50) { // condition determines if the data read in are the sqrt or original burn areas
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
