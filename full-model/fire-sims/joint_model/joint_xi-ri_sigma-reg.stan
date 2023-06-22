functions {
  #include /../../gpd_fcns.stanfunctions
  #include g1_fcns.stanfunctions
  #include /../../twcrps_matnorm_fcns.stanfunctions
}
#include /../../burns_data.stan
transformed data {
  int S = 3; // # of parameters with regression (ranges from 2 to 3)
  //ordering of S: 1 = lambda, 2 = kappa, 3 = sigma
  int C = 6; // # of parameters with correlation (either regression or random intercept)
}
parameters {
  array[N_tb_mis] real<lower=y_min> y_train_mis;
  matrix[R, 3] Z; // ordering: 1 = pi, 2 = delta, 3 = xi
  array[T_all, S] row_vector[R] phi_init;
  array[S] matrix[p, R] beta;
  vector<lower=0>[S] tau_init;
  vector<lower=0, upper = 1>[S] eta;
  vector<lower=0, upper = 1>[S] bp_init;
  // ordering of C: 1 = lambda, 2 = kappa, 3 = sigma, 4 = pi, 5 = delta, 6 = xi
  vector<lower=0, upper = 1>[C] rho1;
  vector<lower=rho1, upper = 1>[C] rho_sum;
  real<lower = 0> theta; // constant shared random effect
  real gamma; // scaling constant for linking theta; only used in burn links
}
transformed parameters {
  array[N_tb_all] real<lower=y_min> y_train;
  vector<lower = 0>[R] delta;
  matrix[T_train, R] lambda;
  vector[R] pi_prob;
  array[S] matrix[T_all, R] phi;
  array[S] matrix[T_train, R] reg;
  vector<lower=0>[S] bp = bp_init / 2;
  vector<lower=0>[S] tau = tau_init / 2;
  vector[C] rho2 = rho_sum - rho1;
  array[S] cov_matrix[p] cov_ar1;
  array[C] corr_matrix[R] corr;
  
  vector[R] ri_init; // random intercept vector
  matrix[T_all, R] ri_matrix; // broadcast ri_init to full matrix
  
  y_train[ii_tb_obs] = y_train_obs;
  y_train[ii_tb_mis] = y_train_mis;
  
  for (c in 1:C) {
    corr[c] = l3 + rho2[c] * l2 + rho1[c] * l1;
  }

  pi_prob = cholesky_decompose(corr[4])' * Z[,1];
  delta = exp(cholesky_decompose(corr[5])' * Z[,2]);  
  ri_init = cholesky_decompose(corr[6])' * Z[,3];
  ri_matrix = rep_matrix(ri_init', T_all);

  for (s in 1:S) {
    cov_ar1[s] = equal + bp[s] * bp_lin + bp[s] ^ 2 * bp_square
                 + bp[s] ^ 3 * bp_cube + bp[s] ^ 4 * bp_quart;
    
    // ICAR variables
    phi[s][1, ] = 1 / tau[s] * phi_init[1, s];
    for (t in 2:T_all) {
      phi[s][t, ] = eta[s] * phi[s][t - 1, ]
                       + 1 / tau[s] * phi_init[t, s];
    }
    
    // regression for kappa, sigma, and xi
    for (r in 1:R) {
      reg[s][, r] = X_train[r] * beta[s][, r] + phi[s][idx_train_er, r];
      lambda[, r] = reg[1][, r] + area_offset[r] + theta;
    }
  }
}
model {
  theta ~ normal(0, 5);
  gamma ~ normal(0, 5);
  
  vector[N_tb_all] kappa = exp(to_vector(reg[2] + gamma * theta))[ii_tb_all];
  vector[N_tb_all] sigma = exp(to_vector(reg[3] + gamma * theta))[ii_tb_all];
  vector[N_tb_all] xi = exp(to_vector(ri_matrix[idx_train_er,]))[ii_tb_all];
  
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
    target += egpd_trunc_lpdf(y_train[n] | y_min, sigma[n], xi[n], kappa[n]);
  }
}

generated quantities {
  array[N_tb_obs] real train_loglik;
  array[N_hold_obs] real holdout_loglik;
  array[N_tb_obs] real train_twcrps;
  array[N_hold_obs] real holdout_twcrps;  

  // condition determines if the data read in are the sqrt or original burn areas
  if (max(y_train_obs) < 30) {
    // declaring the below as "local" variables so they don't contribute to csv size
    array[S] matrix[T_all, R] reg_full;
    for (s in 1:S) {
      for (r in 1:R) {
        reg_full[s][, r] = X_full[r] * beta[s][, r] + phi[s][, r];
      }
    }
    vector[N_tb_obs] kappa_train = exp(to_vector(reg_full[1]))[ii_tb_all][ii_tb_obs];
    vector[N_tb_obs] sigma_train = exp(to_vector(reg_full[2]))[ii_tb_all][ii_tb_obs];
    vector[N_tb_obs] xi_train = exp(to_vector(ri_matrix))[ii_tb_all][ii_tb_obs];
  
    vector[N_hold_obs] kappa_hold = exp(to_vector(reg_full[1]))[ii_hold_all][ii_hold_obs];
    vector[N_hold_obs] sigma_hold = exp(to_vector(reg_full[2]))[ii_hold_all][ii_hold_obs];
    vector[N_hold_obs] xi_hold = exp(to_vector(ri_matrix))[ii_hold_all][ii_hold_obs];
    
    // training scores
    for (n in 1:N_tb_obs) {
      train_loglik[n] = egpd_trunc_lpdf(y_train_obs[n] | y_min, sigma_train[n], xi_train[n], kappa_train[n])
                        + log(0.5) - log(y_train_obs[n]);
      // forecasting then twCRPS, on training dataset
      vector[n_int] pred_probs_train = prob_forecast(n_int, sqrt(int_pts_train), y_min, 
                                                sigma_train[n], xi_train[n], kappa_train[n]);
      train_twcrps[n] = twCRPS((y_train_obs[n])^2, n_int, int_train, int_pts_train, pred_probs_train);
    }
    // holdout scores
    for (n in 1:N_hold_obs) {
      // log-likelihood
      holdout_loglik[n] = egpd_trunc_lpdf(y_hold_obs[n] | y_min, sigma_hold[n], xi_hold[n], kappa_hold[n])
                          + log(0.5) - log(y_hold_obs[n]);
      // forecasting then twCRPS, on holdout dataset
      vector[n_int] pred_probs_hold = prob_forecast(n_int, sqrt(int_pts_holdout), y_min, 
                                                sigma_hold[n], xi_hold[n], kappa_hold[n]);
      holdout_twcrps[n] = twCRPS((y_hold_obs[n])^2, n_int, int_holdout, int_pts_holdout, pred_probs_hold);
    }
  } else {
    array[S] matrix[T_all, R] reg_full;
    for (s in 1:S) {
      for (r in 1:R) {
        reg_full[s][, r] = X_full[r] * beta[s][, r] + phi[s][, r];
      }
    }
    vector[N_tb_obs] kappa_train = exp(to_vector(reg_full[1]))[ii_tb_all][ii_tb_obs];
    vector[N_tb_obs] sigma_train = exp(to_vector(reg_full[2]))[ii_tb_all][ii_tb_obs];
    vector[N_tb_obs] xi_train = exp(to_vector(ri_matrix))[ii_tb_all][ii_tb_obs];
  
    vector[N_hold_obs] kappa_hold = exp(to_vector(reg_full[1]))[ii_hold_all][ii_hold_obs];
    vector[N_hold_obs] sigma_hold = exp(to_vector(reg_full[2]))[ii_hold_all][ii_hold_obs];
    vector[N_hold_obs] xi_hold = exp(to_vector(ri_matrix))[ii_hold_all][ii_hold_obs];
    
    // training scores
    for (n in 1:N_tb_obs) {
      train_loglik[n] = egpd_trunc_lpdf(y_train_obs[n] | y_min, sigma_train[n], xi_train[n], kappa_train[n]);
      // forecasting then twCRPS, on training dataset
      vector[n_int] pred_probs_train = prob_forecast(n_int, int_pts_train, y_min, 
                                              sigma_train[n], xi_train[n], kappa_train[n]);
      train_twcrps[n] = twCRPS(y_train_obs[n], n_int, int_train, int_pts_train, pred_probs_train);
    }
    // holdout scores
    for (n in 1:N_hold_obs) {
      // log-likelihood
      holdout_loglik[n] = egpd_trunc_lpdf(y_hold_obs[n] | y_min, sigma_hold[n], xi_hold[n], kappa_hold[n]);
      // forecasting then twCRPS, on holdout dataset
      vector[n_int] pred_probs_hold = prob_forecast(n_int, int_pts_train, y_min, 
                                            sigma_hold[n], xi_hold[n], kappa_hold[n]);
      holdout_twcrps[n] = twCRPS(y_hold_obs[n], n_int, int_train, int_pts_train, pred_probs_hold);
    }
  }
}