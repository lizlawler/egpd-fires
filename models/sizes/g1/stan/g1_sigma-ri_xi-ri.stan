functions {
  #include /../../gpd_fcns.stanfunctions
  #include g1_fcns.stanfunctions
  #include /../../../twcrps_matnorm_fcns.stanfunctions
}
#include /../../sizes_data.stan
transformed data {
  int S = 1; // # of parameters with regression (ranges from 1 to 3)
  int C = 3; // # of parameters with correlation (either regression or random intercept)
}
parameters {
  array[N_ts_mis] real<lower=y_min> y_train_mis;
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
  array[N_ts_all] real<lower=y_min> y_train;
  array[S] matrix[T_all, R] phi;
  array[S] matrix[T_train, R] reg;
  vector<lower=0>[S] bp = bp_init / 2;
  vector<lower=0>[S] tau = tau_init / 2;
  vector[C] rho2 = rho_sum - rho1;
  array[S] cov_matrix[p] cov_ar1;
  array[C] corr_matrix[R] corr;
  
  array[2] vector[R] ri_init; // random intercept vector
  array[2] matrix[T_all, R] ri_matrix; // broadcast ri_init to full matrix
  
  y_train[ii_ts_obs] = y_train_obs;
  y_train[ii_ts_mis] = y_train_mis;
  
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
  vector[N_ts_all] kappa = exp(to_vector(reg[1]))[ii_ts_all];
  vector[N_ts_all] sigma = exp(to_vector(ri_matrix[1][idx_train_er,]))[ii_ts_all];
  vector[N_ts_all] xi = exp(to_vector(ri_matrix[2][idx_train_er,]))[ii_ts_all];
  
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
  for (n in 1:N_ts_all) {
    target += egpd_trunc_lpdf(y_train[n] | y_min, sigma[n], xi[n], kappa[n]);
  }
}

generated quantities {
  array[N_hold_obs] real train_loglik;
  array[N_hold_obs] real holdout_loglik;
  array[N_hold_obs] real train_twcrps;
  array[N_hold_obs] real holdout_twcrps;

  array[S] matrix[T_all, R] reg_full;
  for (s in 1:S) {
    for (r in 1:R) {
      reg_full[s][, r] = X_full[r] * beta[s][, r] + phi[s][, r];
    }
  }
  
  // training scores
  for (n in 1:N_ts_obs) {
    real kappa_train = exp(to_vector(reg_full[1][idx_train_er,]))[ii_ts_all][ii_ts_obs][n];
    real sigma_train = exp(to_vector(ri_matrix[1][idx_train_er,]))[ii_ts_all][ii_ts_obs][n];
    real xi_train = exp(to_vector(ri_matrix[2][idx_train_er,]))[ii_ts_all][ii_ts_obs][n];

    train_loglik[n] = egpd_trunc_lpdf(y_train_obs[n] | y_min, sigma_train, xi_train, kappa_train);
    // forecasting then twCRPS, on training dataset
    vector[n_int] pred_probs_train = prob_forecast(n_int, int_pts, y_min,
                                            sigma_train, xi_train, kappa_train);
    train_twcrps[n] = twCRPS(y_train_obs[n], n_int, int_range, int_pts, pred_probs_train);
  }
  // holdout scores
  for (n in 1:N_hold_obs) {
    real kappa_hold = exp(to_vector(reg_full[1]))[ii_hold_all][ii_hold_obs][n];
    real sigma_hold = exp(to_vector(ri_matrix[1]))[ii_hold_all][ii_hold_obs][n];
    real xi_hold = exp(to_vector(ri_matrix[2]))[ii_hold_all][ii_hold_obs][n];
    
    // log-likelihood
    holdout_loglik[n] = egpd_trunc_lpdf(y_hold_obs[n] | y_min, sigma_hold, xi_hold, kappa_hold);
    // forecasting then twCRPS, on holdout dataset
    vector[n_int] pred_probs_hold = prob_forecast(n_int, int_pts, y_min, 
                                          sigma_hold, xi_hold, kappa_hold);
    holdout_twcrps[n] = twCRPS(y_hold_obs[n], n_int, int_range, int_pts, pred_probs_hold);
  }
}
