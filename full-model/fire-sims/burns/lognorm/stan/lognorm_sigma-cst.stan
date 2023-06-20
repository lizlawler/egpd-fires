functions {
  #include lognorm_fcns.stanfunctions
}
#include /../../burns_data.stan
parameters {
  array[N_tb_mis] real<lower=y_min> y_train_mis;
  array[T_all] row_vector[R] phi_init;
  matrix[p, R] beta;
  real<lower=0> tau_init;
  real<lower=0, upper = 1> eta;
  real<lower=0, upper = 1> bp_init;
  real<lower=0, upper = 1> rho1;
  real<lower=rho1, upper = 1> rho_sum; 
  real<lower=0> sigma;
}
transformed parameters {
  array[N_tb_all] real<lower=1> y_train;
  matrix[T_all, R] phi;
  matrix[T_train, R] reg;
  real<lower=0> bp = bp_init / 2;
  real<lower=0> tau = tau_init / 2;
  real rho2 = rho_sum - rho1;
  cov_matrix[p] cov_ar1;
  corr_matrix[R] corr; // 1 = mu
  
  y_train[ii_tb_obs] = y_train_obs;
  y_train[ii_tb_mis] = y_train_mis;
  
  corr = l3 + rho2 * l2 + rho1 * l1;
  cov_ar1 = equal + bp * bp_lin + bp ^ 2 * bp_square
                 + bp ^ 3 * bp_cube + bp ^ 4 * bp_quart;
    
  // ICAR variables
  phi[1, ] = 1 / tau * phi_init[1];
  for (t in 2:T_all) {
    phi[t, ] = eta * phi[t - 1, ]
                + 1 / tau * phi_init[t];
  }
  
  // regression on mu
  for (r in 1:R) {
    reg[, r] = X_train[r] * beta[, r] + phi[idx_train_er, r];
  }
}
model {
  vector[N_tb_all] mu = to_vector(reg)[ii_tb_all];
  
  // prior on AR(1) penalization of splines
  bp_init ~ uniform(0, 1);
  
  // priors scaling constants in ICAR
  eta ~ beta(3, 4);
  tau_init ~ exponential(1);

  // prior on rhos
  rho1 ~ beta(3, 4);
  rho_sum ~ beta(8, 2);

  // MVN prior on betas
  target += matnormal_lpdf(beta | cov_ar1, corr);
  // ICAR prior
  for (t in 1:T_all) {
    target += -.5 * dot_self(phi_init[t][node1] - phi_init[t][node2]);
    sum(phi_init[t]) ~ normal(0, 0.001 * R);
  }
  
  // prior on scale parameter
  sigma ~ student_t(4, 0, 1);
  
  // likelihood
  for (n in 1:N_tb_all) {
    target += lognorm_trunc_lpdf(y_train[n] | y_min, mu[n], sigma);
  }
}
generated quantities {
  array[N_tb_obs] real train_loglik;
  array[N_hold_obs] real holdout_loglik;
  array[N_tb_obs] real train_twcrps;
  array[N_hold_obs] real holdout_twcrps;
  
  // condition determines if the data read in are the sqrt or original burn areas  
  if (max(y_train_obs) < 30) {
    matrix[T_all, R] reg_full;
    for (r in 1:R) {
      reg_full[, r] = X_full[r] * beta[, r] + phi[, r];
    }
    vector[N_tb_obs] mu_train = to_vector(reg_full)[ii_tb_all][ii_tb_obs];
    vector[N_hold_obs] mu_hold = to_vector(reg_full)[ii_hold_all][ii_hold_obs];
    
    // training scores
    for (n in 1:N_tb_obs) {
      train_loglik[n] = lognorm_trunc_lpdf(y_train_obs[n] | y_min, mu_train[n], sigma)
                        + log(0.5) - log(y_train_obs[n]);
      // forecasting then twCRPS, on training dataset
      vector[n_int] pred_probs_train = prob_forecast(n_int, sqrt(int_pts_train), 
                                        y_min, mu_train[n], sigma);
      train_twcrps[n] = twCRPS((y_train_obs[n])^2, n_int, int_train, int_pts_train, pred_probs_train);
    }
    // holdout scores
    for (n in 1:N_hold_obs) {
      // log-likelihood
      holdout_loglik[n] = lognorm_trunc_lpdf(y_hold_obs[n] | y_min, mu_hold[n], sigma)
                          + log(0.5) - log(y_hold_obs[n]);
      // forecasting then twCRPS, on holdout dataset
      vector[n_int] pred_probs_hold = prob_forecast(n_int, sqrt(int_pts_holdout), 
                                        y_min, mu_hold[n], sigma);
      holdout_twcrps[n] = twCRPS((y_hold_obs[n])^2, n_int, int_holdout, int_pts_holdout, pred_probs_hold);
    }
  } else {
    matrix[T_all, R] reg_full;
    for (r in 1:R) {
      reg_full[, r] = X_full[r] * beta[, r] + phi[, r];
    }
    vector[N_tb_obs] mu_train = to_vector(reg_full)[ii_tb_all][ii_tb_obs];
    vector[N_hold_obs] mu_hold = to_vector(reg_full)[ii_hold_all][ii_hold_obs];
    
    // training scores
    for (n in 1:N_tb_obs) {
      train_loglik[n] = lognorm_trunc_lpdf(y_train_obs[n] | y_min, mu_train[n], sigma);
      // forecasting then twCRPS, on training dataset
      vector[n_int] pred_probs_train = prob_forecast(n_int, int_pts_train, 
                                        y_min, mu_train[n], sigma);
      train_twcrps[n] = twCRPS(y_train_obs[n], n_int, int_train, int_pts_train, pred_probs_train);
    }
    // holdout scores
    for (n in 1:N_hold_obs) {
      // log-likelihood
      holdout_loglik[n] = lognorm_trunc_lpdf(y_hold_obs[n] | y_min, mu_hold[n], sigma);
      // forecasting then twCRPS, on holdout dataset
      vector[n_int] pred_probs_hold = prob_forecast(n_int, int_pts_holdout, 
                                        y_min, mu_hold[n], sigma);
      holdout_twcrps[n] = twCRPS(y_hold_obs[n], n_int, int_holdout, int_pts_holdout, pred_probs_hold);
    }
  }
}

