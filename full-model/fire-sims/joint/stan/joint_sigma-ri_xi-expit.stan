functions {
  #include /../../burns/gpd_fcns.stanfunctions
  #include /../../burns/g1/stan/g1_fcns.stanfunctions
  #include /../../burns/twcrps_matnorm_fcns.stanfunctions
}
#include /../joint_data.stan
transformed data {
  int S = 2; // # of parameters with regression (ranges from 2 to 3)
  //ordering of S: 1 = lambda, 2 = kappa
  int C = 7; // # of parameters with correlation (either regression or random intercept)
  // ordering: 1=lambda, 2=kappa, 3=pi, 4=delta, 5=sigma, 6=xi, 7 = gamma
}
parameters {
  array[N_tb_mis] real<lower=y_min> y_train_burn_mis;
  matrix[R, C-S] Z; // ordering: 1 = pi, 2 = delta, 3 = sigma, 4 = xi, 5 = gamma
  array[T_all, S] row_vector[R] phi_init;
  matrix[p, R] beta_count;
  matrix[p_burn, R] beta_burn;
  vector<lower=0>[S] tau_init;
  vector<lower=0, upper = 1>[S+1] eta; // additional eta for theta AR(1) prior
  vector<lower=0, upper = 1>[S] bp_init; // AR(1) param for covariance matrix
  vector<lower=0, upper = 1>[C] rho1;
  vector<lower=rho1, upper = 1>[C] rho_sum;
  vector[T_all] theta_init; // shared random effect, varying in time
  real<lower=0> eps; // noise for AR(1) prior on theta
}
transformed parameters {
  array[N_tb_all] real<lower=y_min> y_train_burn;
  vector<lower = 0>[R] delta;
  vector[R] gamma; // scaling parameter, varying by region
  matrix[T_train, R] lambda;
  vector[R] pi_prob;
  array[S] matrix[T_all, R] phi;
  matrix[T_train, R] reg;
  vector<lower=0>[S] bp = bp_init / 2;
  vector<lower=0>[S] tau = tau_init / 2;
  vector[C] rho2 = rho_sum - rho1;
  array[S] cov_matrix[p] cov_ar1;
  array[C] corr_matrix[R] corr;
  vector[T_all] theta;
  
  array[2] vector[R] ri_init; // random intercept vector
  array[2] matrix[T_all, R] ri_matrix; // broadcast ri_init to full matrix
  
  y_train_burn[ii_tb_obs] = y_train_burn_obs;
  y_train_burn[ii_tb_mis] = y_train_burn_mis;
  
  for (c in 1:C) {
    corr[c] = l3 + rho2[c] * l2 + rho1[c] * l1;
  }

  pi_prob = cholesky_decompose(corr[3])' * Z[,1];
  delta = exp(cholesky_decompose(corr[4])' * Z[,2]);  
  for (i in 1:2) {
    ri_init[i] = cholesky_decompose(corr[i+4])' * Z[,i+2];
    ri_matrix[i] = rep_matrix(ri_init[i]', T_all);
  }
  gamma = cholesky_decompose(corr[7])' * Z[,5];

  for (s in 1:S) {
    cov_ar1[s] = equal + bp[s] * bp_lin + bp[s] ^ 2 * bp_square
                 + bp[s] ^ 3 * bp_cube + bp[s] ^ 4 * bp_quart;
    
    // temporal CAR
    phi[s][1, ] = 1 / tau[s] * phi_init[1, s];
    for (t in 2:T_all) {
      phi[s][t, ] = eta[s] * phi[s][t - 1, ]
                       + 1 / tau[s] * phi_init[t, s];
    }
  }

  // temporal AR(1)
  theta[1] = theta_init[1];
  for (t in 2:T_all) {
    theta[t] = eta[S+1] * theta[t-1] + theta_init[t];
  }
  
  // regression link for lambda (counts) and kappa (burns)
  for (r in 1:R) {
    lambda[, r] = X_train_count[r] * beta_count[, r] + phi[1][idx_train_er, r] + area_offset[r] + theta[idx_train_er];
    reg[, r] = X_train_burn[r] * beta_burn[, r] + phi[2][idx_train_er, r] + gamma[r] * theta[idx_train_er];
  }
}
model {
  vector[N_tb_all] kappa = exp(to_vector(reg))[ii_tb_all];
  vector[N_tb_all] sigma = exp(to_vector(ri_matrix[1][idx_train_er,]))[ii_tb_all];
  vector[N_tb_all] xi = (inv_logit(to_vector(ri_matrix[2][idx_train_er,]))[ii_tb_all]) * 1.5;
  
  to_vector(Z) ~ normal(0, 1.5);
  
  // prior on AR(1) penalization of splines
  to_vector(bp_init) ~ uniform(0, 1);
  
  // priors for scaling constants in temporal CAR
  to_vector(eta) ~ beta(3, 4);
  to_vector(tau_init) ~ exponential(1);

  // initialize prior on theta
  eps ~ student_t(4, 0, 1);
  to_vector(theta_init) ~ normal(0, eps);

  // prior on rhos
  to_vector(rho1) ~ beta(3, 4);
  to_vector(rho_sum) ~ beta(8, 2);
  
  // MVN prior on betas
  // first, for lambda - has 37 covariates
  target += matnormal_lpdf(beta_count | cov_ar1[1], corr[1]);
  // then, for kappa and simga - each have 13 covariates
  target += matnormal_lpdf(beta_burn | cov_ar1[2][1:p_burn, 1:p_burn], corr[2]);
  for (s in 1:S) {
    // ICAR spatial prior
    for (t in 1:T_all) {
      target += -.5 * dot_self(phi_init[t, s][node1] - phi_init[t, s][node2]);
      sum(phi_init[t, s]) ~ normal(0, 0.001 * R);
    }
  }
  
  // burn likelihood
  for (n in 1:N_tb_all) {
    target += egpd_trunc_lpdf(y_train_burn[n] | y_min, sigma[n], xi[n], kappa[n]);
  }

  // count likelihood
  for (r in 1:R) {
    for (t in 1:T_train) {
      if (y_train_count[t, r] == 0) {
        target += log_sum_exp(bernoulli_logit_lpmf(1 | pi_prob[r]),
                              bernoulli_logit_lpmf(0 | pi_prob[r])
                              + neg_binomial_2_log_lpmf(y_train_count[t, r] | lambda[t, r], delta[r]));
      } else {
        target += bernoulli_logit_lpmf(0 | pi_prob[r])
                  + neg_binomial_2_log_lpmf(y_train_count[t, r] | lambda[t, r], delta[r]);
      }
    }
  }
}

generated quantities {
  matrix[T_train, R] train_loglik_count;
  matrix[T_hold, R] holdout_loglik_count;
  array[N_tb_obs] real train_loglik_burn;
  array[N_hold_obs] real holdout_loglik_burn;
  array[N_tb_obs] real train_twcrps;
  array[N_hold_obs] real holdout_twcrps;  

  matrix[T_all, R] reg_full;
  matrix[T_all, R] lambda_full;
  matrix[T_all, R] burn_pred;

  // expected value of all parameters based on all timepoints, then cut to only be holdout parameters
  for (r in 1:R) {
    lambda_full[, r] = X_full_count[r] * beta_count[, r] + phi[1][, r] + theta;
    reg_full[, r] = X_full_burn[r] * beta_burn[, r] + phi[2][, r] + gamma[r] * theta;
  }
  
  // burn component scores
  // training scores
  for (n in 1:N_tb_obs) {
    real kappa_train = exp(to_vector(reg_full[idx_train_er,]))[ii_tb_all][ii_tb_obs][n];
    real sigma_train = exp(to_vector(ri_matrix[1][idx_train_er,]))[ii_tb_all][ii_tb_obs][n];
    real xi_train = (inv_logit(to_vector(ri_matrix[2][idx_train_er,]))[ii_tb_all][ii_tb_obs][n]) * 1.5;
    
    train_loglik_burn[n] = egpd_trunc_lpdf(y_train_burn_obs[n] | y_min, sigma_train, xi_train, kappa_train);
    // forecasting then twCRPS, on training dataset
    vector[n_int] pred_probs_train = prob_forecast(n_int, int_pts_train, y_min, 
                                            sigma_train, xi_train, kappa_train);
    train_twcrps[n] = twCRPS(y_train_burn_obs[n], n_int, int_train, int_pts_train, pred_probs_train);
  }
  // holdout scores
  for (n in 1:N_hold_obs) {
    real kappa_hold = exp(to_vector(reg_full))[ii_hold_all][ii_hold_obs][n];
    real sigma_hold = exp(to_vector(ri_matrix[1]))[ii_hold_all][ii_hold_obs][n];
    real xi_hold = (inv_logit(to_vector(ri_matrix[2]))[ii_hold_all][ii_hold_obs][n]) * 1.5;
    
    // log-likelihood
    holdout_loglik_burn[n] = egpd_trunc_lpdf(y_hold_burn_obs[n] | y_min, sigma_hold, xi_hold, kappa_hold);
      // forecasting then twCRPS, on holdout dataset
    vector[n_int] pred_probs_hold = prob_forecast(n_int, int_pts_holdout, y_min, 
                                          sigma_hold, xi_hold, kappa_hold);
    holdout_twcrps[n] = twCRPS(y_hold_burn_obs[n], n_int, int_holdout, int_pts_holdout, pred_probs_hold);
  }

  // count component scores
  // training log-likelihood
  for (r in 1:R) {
    for (t in 1:T_train) {
      if (y_train_count[t, r] == 0) {
        train_loglik_count[t, r] = log_sum_exp(bernoulli_logit_lpmf(1 | pi_prob[r]),
                                         bernoulli_logit_lpmf(0 | pi_prob[r])
                                         + neg_binomial_2_log_lpmf(y_train_count[t, r] | lambda[t, r], delta[r]));
      } else {
        train_loglik_count[t, r] = bernoulli_logit_lpmf(0 | pi_prob[r])
                             + neg_binomial_2_log_lpmf(y_train_count[t, r] | lambda[t, r], delta[r]);
      }
    }
  }
  
  // holdout log-likelihood
  for (r in 1:R) {
    vector[T_hold] lambda_hold = lambda_full[idx_hold_er, r];
    for (t in 1:T_hold) {
      if (y_hold_count[t, r] == 0) {
        holdout_loglik_count[t, r] = log_sum_exp(bernoulli_logit_lpmf(1 | pi_prob[r]),
                                           bernoulli_logit_lpmf(0 | pi_prob[r])
                                           + neg_binomial_2_log_lpmf(y_hold_count[t, r] | lambda_hold[t], delta[r]));
      } else {
        holdout_loglik_count[t, r] = bernoulli_logit_lpmf(0 | pi_prob[r])
                               + neg_binomial_2_log_lpmf(y_hold_count[t, r] | lambda_hold[t], delta[r]);
      }
    }
  }
  
  // generate burn area predictions, based on how many fires occurred (draw from ZINB)
  for (r in 1:R) {
    for (t in 1:T_all) {
      vector[500] burn_draws;
      real sigma = exp(ri_init[1][r]);
      real xi = inv_logit(ri_init[2][r]) * 1.5;
      real kappa = exp(reg_full[t, r]);
      for (i in 1:500) {
        int zero = bernoulli_logit_rng(pi_prob[r]);
        int count_draw = (1 - zero) * neg_binomial_2_log_rng(lambda_full[t, r], delta[r]);
        if (count_draw == 0) {
          burn_draws[i] = 0; 
        } else {
          burn_draws[i] = sum(egpd_rng(count_draw, y_min, sigma, xi, kappa));
        }
      }
      burn_pred[t, r] = mean(burn_draws);
    }
  }
}
