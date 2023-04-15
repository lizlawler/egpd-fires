functions {
  // forecast_rng and egpd_lpdf vary by model
  vector forecast_rng(int n_pred, real ymin, real kappa, real sigma, real xi, real gamma) {
    vector[n_pred] forecast;
    vector[n_pred] a = rep_vector(0, n_pred);
    vector[n_pred] b = rep_vector(1, n_pred);
    array[n_pred] real u = uniform_rng(a, b);
    real alpha = 1/gamma;
    real beta = 2;
    real cst_term = (1 + xi * (ymin/sigma))^(-gamma/xi);
    real cst = 1 - beta_cdf(cst_term | alpha, beta);
    for (n in 1:n_pred) {
      real u_adj = u[n] * (1 - cst) + cst;
      forecast[n] = (sigma / xi) * ((inv_inc_beta(alpha, beta, 1-u_adj^(2/kappa))) ^ (-xi/gamma) - 1);
    }
    return forecast;
  }
  
  real egpd_lpdf(real y, real ymin, real sigma, real xi, real gamma, real kappa) {
    real alpha = 1/gamma;
    real beta = 2;
    real w = 1 + xi * (y/sigma);
    real cst_term = (1 + xi * (ymin/sigma))^(-gamma/xi);
    real cst = (1 - beta_cdf(cst_term | alpha, beta))^(kappa/2);
    real lpdf = log(kappa) - log(2) + (kappa/2 - 1) * log1m(beta_cdf(w^(-gamma/xi) | alpha, beta)) + 
      log(gamma) - log(sigma) + beta_lpdf(w^(-gamma/xi) | alpha, beta) - (gamma/xi + 1) * log(w);
    return lpdf - log1m(cst);
  }
  
  // twCRPS and matnormal_lpdf remain unchanged across models
  real twCRPS(real y, vector forecast, real delta, real w_mean, real w_sd) {
    real score;
    real summand;
    int N = rows(forecast);
    summand = 0;
    for (n in 1:N) {
      summand += (forecast[n] - step(forecast[n] - y)) ^ 2
                 * normal_cdf(forecast[n] | w_mean, w_sd);
    }
    score = summand * delta;
    return score;
  }
  
  real matnormal_lpdf(matrix y, matrix cov, matrix corr) {
    real lpdf;
    real r = rows(corr);
    real p = rows(cov);
    lpdf = -(r * p / 2) * log(2 * pi()) - (p / 2) * log_determinant(corr)
           - (r / 2) * log_determinant(cov)
           - 0.5 * trace(mdivide_right_spd(mdivide_left_spd(corr, y'), cov) * y);
    return lpdf;
  }
}
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
  int C = 4; // # of parameters with correlation (either regression or random intercept)
}
parameters {
  array[N_tb_mis] real<lower=y_min> y_train_mis;
  matrix[R, 3] Z; //1=nu, 2=xi, 3=gamma
  array[T_all, S] row_vector[R] phi_init;
  array[S] matrix[p, R] beta;
  vector<lower=0>[S] tau_init;
  vector<lower=0, upper = 1>[S] eta;
  vector<lower=0, upper = 1>[S] bp_init;
  vector<lower=0, upper = 1>[C] rho1;
  vector<lower=rho1, upper = 1>[C] rho_sum; // 1 = kappa, 2 = nu, 3 = xi, 4 = gamma
}
transformed parameters {
  array[N_tb_all] real<lower=y_min> y_train;
  array[S] matrix[T_all, R] phi;
  array[S] matrix[T_train, R] reg;
  vector<lower=0>[S] bp = bp_init / 2;
  vector<lower=0>[S] tau = tau_init / 2;
  vector[C] rho2 = rho_sum - rho1;
  array[S] cov_matrix[p] cov_ar1;
  array[C] corr_matrix[R] corr; // 1 = kappa, 2= nu, 3 = xi, 4 = gamma
  
  array[3] vector[R] ri_init; 
  array[3] matrix[T_all, R] ri_matrix; 
  
  y_train[ii_tb_obs] = y_train_obs;
  y_train[ii_tb_mis] = y_train_mis;
  
  for (c in 1:C) {
    corr[c] = l3 + rho2[c] * l2 + rho1[c] * l1;
  }
  
  for(i in 1:3) {
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
    
    // regression for gamma, nu, and xi
    for (r in 1:R) {
      reg[s][, r] = X_train[r] * beta[s][, r] + phi[s][idx_train_er, r];
    }
  }
}
model {
  vector[N_tb_all] kappa = exp(to_vector(reg[1]))[ii_tb_all];
  vector[N_tb_all] nu = exp(to_vector(ri_matrix[1][idx_train_er,]))[ii_tb_all];
  vector[N_tb_all] xi = exp(to_vector(ri_matrix[2][idx_train_er,]))[ii_tb_all];
  vector[N_tb_all] gamma = exp(to_vector(ri_matrix[3][idx_train_er,]))[ii_tb_all];
  vector[N_tb_all] sigma = nu ./ (1 + xi);
  
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
    target += egpd_lpdf(y_train[n] | y_min, sigma[n], xi[n], gamma[n], kappa[n]);
  }
}

// generated quantities {
//   array[S] matrix[T_all, R] reg_full;
//  
//   vector<lower=0>[N_tb_obs] kappa_train;
//   vector<lower=0>[N_tb_obs] nu_train;
//   vector<lower=0>[N_tb_obs] xi_train;
//   vector<lower=0>[N_tb_obs] gamma_train;
//   vector<lower=0>[N_tb_obs] sigma_train;
//  
//   vector<lower=0>[N_hold_obs] kappa_hold;
//   vector<lower=0>[N_hold_obs] nu_hold;
//   vector<lower=0>[N_hold_obs] xi_hold;
//   vector<lower=0>[N_hold_obs] gamma_hold;
//   vector<lower=0>[N_hold_obs] sigma_hold;
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
//   nu_train = exp(to_vector(ri_matrix[1]))[ii_tb_all][ii_tb_obs];
//   xi_train = exp(to_vector(ri_matrix[2]))[ii_tb_all][ii_tb_obs];
//   gamma_train = exp(to_vector(ri_matrix[3]))[ii_tb_all][ii_tb_obs];
//   sigma_train = nu_train ./ (1 + xi_train);
//   
//   kappa_hold = exp(to_vector(reg_full[1]))[ii_hold_all][ii_hold_obs];
//   nu_hold = exp(to_vector(ri_matrix[1]))[ii_hold_all][ii_hold_obs];
//   xi_hold = exp(to_vector(ri_matrix[2]))[ii_hold_all][ii_hold_obs];
//   gamma_hold = exp(to_vector(ri_matrix[3]))[ii_hold_all][ii_hold_obs];
//   sigma_hold = nu_hold ./ (1 + xi_hold);
//   
//   if (max(y_train_obs) < 50) {
//     // condition determines if the data read in are the sqrt or original burn areas
//     // training log-likelihood
//     for (n in 1:N_tb_obs) {
//       train_loglik[n] = egpd_lpdf(y_train_obs[n] | sigma_train[n], xi_train[n], gamma_train[n], kappa_train[n])
//                         + log(0.5) - log(y_train_obs[n]);
//     }
//     // holdout scores
//     for (n in 1:N_hold_obs) {
//       // log-likelihood
//       holdout_loglik[n] = egpd_lpdf(y_hold_obs[n] | sigma_hold[n], xi_hold[n], gamma_hold[n], kappa_hold[n])
//                           + log(0.5) - log(y_hold_obs[n]);
//       // twCRPS
//       holdout_twcrps[n] = twCRPS(y_hold_obs[n],
//                                  forecast_rng(n_pred, sigma_hold[n],
//                                               xi_hold[n], gamma_hold[n], kappa_hold[n]),
//                                  delta, sqrt(21), 3);
//     }
//   } else {
//     // training log-likelihood
//     for (n in 1:N_tb_obs) {
//       train_loglik[n] = egpd_lpdf(y_train_obs[n] | sigma_train[n], xi_train[n], gamma_train[n], kappa_train[n]);
//     }
//     // holdout scores
//     for (n in 1:N_hold_obs) {
//       // log-likelihood
//       holdout_loglik[n] = egpd_lpdf(y_hold_obs[n] | sigma_hold[n], xi_hold[n], gamma_hold[n], kappa_hold[n]);
//       // twCRPS 
//       holdout_twcrps[n] = twCRPS(y_hold_obs[n],
//                                  forecast_rng(n_pred, sigma_hold[n],
//                                               xi_hold[n], gamma_hold[n], kappa_hold[n]),
//                                  delta, 21, 9);
//     }
//   }
// }
// 
// 
