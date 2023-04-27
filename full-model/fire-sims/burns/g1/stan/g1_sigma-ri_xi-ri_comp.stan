#include g1_fcns.stan
data {
  int<lower=1> R; // # of regions
  int<lower=1> p; // # of parameters
  int<lower=1> T_all; // # of timepoints in full dataset
  int<lower=1> T_train;
  int<lower=1> T_hold;
  
  // covariate data
  array[R] matrix[T_all, p] X_full; // design matrix; 1-D array of size r with matrices t x p
  array[R] matrix[T_train, p] X_train; // design matrix; 1-D array of size r with matrices t x p
  
  // lower bound of burns
  real y_min;

  // training data
  int<lower=1> N_tb_obs;
  int<lower=1> N_tb_mis;
  int<lower=1> N_tb_all;
  array[N_tb_obs] real<lower=y_min> y_train_obs; // burn area for observed training timepoints
  array[N_tb_obs] int<lower=1> ii_tb_obs;
  array[N_tb_mis] int<lower=1, upper=N_tb_all> ii_tb_mis;
  array[N_tb_all] int<lower=1, upper=N_tb_all> ii_tb_all; // for broadcasting
  array[T_train] int<lower=1> idx_train_er;
  
  // holdout data
  int<lower=1> N_hold_obs;
  int<lower=1> N_hold_all; // includes 'missing' and observed
  array[N_hold_obs] int<lower=1> ii_hold_obs; // vector of indices for holdout data timepoints
  array[N_hold_all] int<lower=1> ii_hold_all; // vector of indices for broadcasting to entire holdout dataset
  array[N_hold_obs] real<lower=1> y_hold_obs; // burn area for observed holdout timepoints
  array[T_hold] int<lower=1> idx_hold_er;
  
  // neighbor information
  int<lower=0> n_edges;
  array[n_edges] int<lower=1, upper=R> node1; // node1[i] adjacent to node2[i]
  array[n_edges] int<lower=1, upper=R> node2; // and node1[i] < node2[i]
  
  // indicator matrices for ecoregions
  matrix[R, R] l3;
  matrix[R, R] l2;
  matrix[R, R] l1;
  
  // indicator matrices for AR(1) penalization on spline coefficients of betas
  matrix[p, p] equal;
  matrix[p, p] bp_lin;
  matrix[p, p] bp_square;
  matrix[p, p] bp_cube;
  matrix[p, p] bp_quart;
}
transformed data {
  int S = 1; // # of parameters with regression (ranges from 1 to 3)
  int C = 3; // # of parameters with correlation (either regression or random intercept)
}
parameters {
  array[N_tb_mis] real<lower=y_min> y_train_mis;
  matrix[R, 2] Z; // 1 = sigma, 2 = xi
  array[T_all, S] row_vector[R] phi_init;
  array[S] matrix[p, R] beta;
  vector<lower=0>[S] tau_init;
  vector<lower=0, upper = 1>[S] eta;
  vector<lower=0, upper = 1>[S] bp_init;
  vector<lower=0, upper = 1>[C] rho1; // 1 = kappa, 2 = sigma, 3 = xi
  vector<lower=rho1, upper = 1>[C] rho_sum;
}
transformed parameters {
  array[N_tb_all] real<lower=y_min> y_train;
  array[S] matrix[T_all, R] phi;
  array[S] matrix[T_train, R] reg;
  vector<lower=0>[S] bp = bp_init / 2;
  vector<lower=0>[S] tau = tau_init / 2;
  vector[C] rho2 = rho_sum - rho1;
  array[S] cov_matrix[p] cov_ar1;
  array[C] corr_matrix[R] corr;
  
  array[2] vector[R] ri_init; // random intercept vector
  array[2] matrix[T_all, R] ri_matrix; // broadcast ri_init to full matrix
  
  y_train[ii_tb_obs] = y_train_obs;
  y_train[ii_tb_mis] = y_train_mis;
  
  for (c in 1:C) {
    corr[c] = l3 + rho2[c] * l2 + rho1[c] * l1;
  }
  
  for (i in 1:2) {
    ri_init[i] = cholesky_decompose(corr[i+1])' * Z[,i];
    ri_matrix[i] = rep_matrix(ri_init[i]', T_all);
  }
  
  for (s in 1:S) {
    cov_ar1[s] = equal + bp[s] * bp_lin + bp[s] ^ 2 * bp_square
                 + bp[s] ^ 3 * bp_cube + bp[s] ^ 4 * bp_quart;
    
    // ICAR variables
    phi[s][1, ] = 1 / tau[s] * phi_init[1, s];
    for (t in 2:T_all) {
      phi[s][t, ] = eta[s] * phi[s][t - 1, ]
                       + 1 / tau[s] * phi_init[t, s];
    }
    
    // regression for kappa, nu, and xi
    for (r in 1:R) {
      reg[s][, r] = X_train[r] * beta[s][, r] + phi[s][idx_train_er, r];
    }
  }
}

model {
  vector[N_tb_all] kappa = exp(to_vector(reg[1]))[ii_tb_all];
  vector[N_tb_all] sigma = exp(to_vector(ri_matrix[1][idx_train_er,]))[ii_tb_all];
  vector[N_tb_all] xi = exp(to_vector(ri_matrix[2][idx_train_er,]))[ii_tb_all];
  
  to_vector(Z) ~ std_normal();
  
  // prior on AR(1) penalization of splines
  to_vector(bp_init) ~ uniform(0, 1);
  
  // priors scaling constants in ICAR
  to_vector(eta) ~ beta(3, 4);
  to_vector(tau_init) ~ exponential(1);

  // prior on rhos
  to_vector(rho1) ~ beta(3, 4);
  to_vector(rho_sum) ~ beta(8, 2);
  
  for (s in 1:S) {
    // MVN prior on betas
    target += matnormal_lpdf(beta[s] | cov_ar1[s], corr[s]);
    // ICAR prior
    for (t in 1:T_all) {
      target += -.5 * dot_self(phi_init[t, s][node1] - phi_init[t, s][node2]);
      sum(phi_init[t, s]) ~ normal(0, 0.001 * R);
    }
  }
  
  // likelihood
  for (n in 1:N_tb_all) {
    target += egpd_cens_lpdf(y_train[n] | y_min, sigma[n], xi[n], kappa[n]);
  }
}

// generated quantities {
//   array[S] matrix[T_all, R] reg_full;
//   
//   vector<lower=0>[N_tb_obs] kappa_train;
//   vector<lower=0>[N_tb_obs] sigma_train;
//   vector<lower=0>[N_tb_obs] xi_train;
//   
//   vector<lower=0>[N_hold_obs] kappa_hold;
//   vector<lower=0>[N_hold_obs] sigma_hold;
//   vector<lower=0>[N_hold_obs] xi_hold;
//   
//   array[N_tb_obs] real train_loglik;
//   array[N_hold_obs] real holdout_loglik;
//   array[N_hold_obs] real holdout_twcrps;
//   
//   // variables needed for estimation of twCRPS integral via summation
//   real interval = max(y_hold_obs) - min(y_hold_obs);
//   int n_pred = 10000;
//   real delta = interval / n_pred;
//   
//   for (s in 1:S) {
//     for (r in 1:R) {
//       reg_full[s][, r] = X_full[r] * beta[s][, r] + phi[s][, r];
//     }
//   }
//   
//   kappa_train = exp(to_vector(reg_full[1]))[ii_tb_all][ii_tb_obs];
//   sigma_train = exp(to_vector(ri_matrix[1]))[ii_tb_all][ii_tb_obs];
//   xi_train = exp(to_vector(ri_matrix[2]))[ii_tb_all][ii_tb_obs];
//   
//   kappa_hold = exp(to_vector(reg_full[1]))[ii_hold_all][ii_hold_obs];
//   sigma_hold = exp(to_vector(ri_matrix[1]))[ii_hold_all][ii_hold_obs];
//   xi_hold = exp(to_vector(ri_matrix[2]))[ii_hold_all][ii_hold_obs];
//   
//   if (max(y_train_obs) < 50) {
//     // condition determines if the data read in are the sqrt or original burn areas
//     // training log-likelihood
//     for (n in 1:N_tb_obs) {
//       train_loglik[n] = egpd_lpdf(y_train_obs[n] | sigma_train[n], xi_train[n], kappa_train[n])
//                         + log(0.5) - log(y_train_obs[n]);
//     }
//     // holdout scores
//     for (n in 1:N_hold_obs) {
//       // log-likelihood
//       holdout_loglik[n] = egpd_lpdf(y_hold_obs[n] | sigma_hold[n], xi_hold[n], kappa_hold[n])
//                           + log(0.5) - log(y_hold_obs[n]);
//       // twCRPS
//       holdout_twcrps[n] = twCRPS(y_hold_obs[n],
//                                  forecast_rng(n_pred, sigma_hold[n],
//                                               xi_hold[n], kappa_hold[n]),
//                                  delta, sqrt(21), 3);
//     }
//   } else {
//     // training log-likelihood
//     for (n in 1:N_tb_obs) {
//       train_loglik[n] = egpd_lpdf(y_train_obs[n] | sigma_train[n], xi_train[n], kappa_train[n]);
//     }
//     // holdout scores
//     for (n in 1:N_hold_obs) {
//       // log-likelihood
//       holdout_loglik[n] = egpd_lpdf(y_hold_obs[n] | sigma_hold[n], xi_hold[n], kappa_hold[n]);
//       // twCRPS 
//       holdout_twcrps[n] = twCRPS(y_hold_obs[n],
//                                  forecast_rng(n_pred, sigma_hold[n],
//                                               xi_hold[n], kappa_hold[n]),
//                                  delta, 21, 9);
//     }
//   }
// }
// 
// 
