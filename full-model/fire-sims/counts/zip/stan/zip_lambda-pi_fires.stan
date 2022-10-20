functions{
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
  int<lower = 1> p; // # of parameters
  int<lower = 1> r; // # of regions
  int<lower=0> n_edges;
  int<lower=1, upper = r> node1[n_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper = r> node2[n_edges];  // and node1[i] < node2[i]
  
  // training dataset
  int<lower = 1> t_tc; // # of training timepoints
  int<lower = 1> idx_tc_er[t_tc]; // vector of indices for training data timepoints
  matrix[t_tc, p] X_tc[r]; // design matrix; 1-D array of size r with matrices t x p; this indexing is deprecated in version 2.32 and above
  int<lower = 0> y_tc[t_tc, r]; // response data
  real area_offset[r]; // known offset vector of areas
  
  // holdout dataset
  int<lower = 1> t_ho; // # of holdout timepoints
  int<lower = 1> idx_hold_er[t_ho]; // vector of indices for holdout data timepoints
  int<lower = 0> y_hold[t_ho, r]; // response data
  
  // full dataset
  int<lower = 1> t_all; // # timepoints in entire dataset
  matrix[t_all, p] X_all_tmpt[r]; // design matrix; 1-D array of size r with matrices t x p; this indexing is deprecated in version 2.32 and above
  
  // indicator matrices for ecoregions
  matrix[r, r] l3;
  matrix[r, r] l2;
  matrix[r, r] l1;

  // indicator matrices for Ar(1) process on betas
  matrix[p, p] equal;
  matrix[p, p] bp_lin;
  matrix[p, p] bp_square;
  matrix[p, p] bp_cube;
  matrix[p, p] bp_quart;
}

parameters {
  vector[r] phi_init_lambda[t_all];
  vector[r] phi_init_pi[t_all];
  matrix[p, r] beta_lambda;
  matrix[p, r] beta_pi;
  real<lower = 0> tau_init_lambda;
  real<lower = 0, upper = 1> eta_lambda;
  real<lower = 0, upper = 1> bp_init_lambda;
  real<lower = 0, upper = 1> rho2_lambda;
  real<lower=0, upper = (1-rho2_lambda)> rho1_lambda;
  real<lower = 0> tau_init_pi;
  real<lower = 0, upper = 1> eta_pi;
  real<lower = 0, upper = 1> bp_init_pi;
  real<lower = 0, upper = 1> rho2_pi;
  real<lower=0, upper = (1-rho2_pi)> rho1_pi;
}

transformed parameters {
  matrix[t_all, r] phi_lambda;
  matrix[t_all, r] phi_pi;
  matrix[t_tc, r] lambda;
  matrix[t_tc, r] pi_param;
  
  real<lower=0, upper = bp_init_lambda/2> bp_lambda = bp_init_lambda/2;
  real<lower=0, upper = tau_init_lambda/2> tau_lambda = tau_init_lambda/2;
  real<lower=0, upper = bp_init_pi/2> bp_pi = bp_init_pi/2;
  real<lower=0, upper = tau_init_pi/2> tau_pi = tau_init_pi/2;
  
  matrix[r, r] corr_lambda = l3 + rho2_lambda * l2 + rho1_lambda * l1;
  matrix[p, p] cov_ar1_lambda = equal + bp_lambda * bp_lin + bp_lambda^2 * bp_square + bp_lambda^3 * bp_cube + bp_lambda^4 * bp_quart;
  matrix[r, r] corr_pi = l3 + rho2_pi * l2 + rho1_pi * l1;
  matrix[p, p] cov_ar1_pi = equal + bp_pi * bp_lin + bp_pi^2 * bp_square + bp_pi^3 * bp_cube + bp_pi^4 * bp_quart;

  phi_lambda[1,] = (1/tau_lambda) * phi_init_lambda[1]';
  phi_pi[1,] = (1/tau_pi) * phi_init_pi[1]';
  for (j in 2:t_all) {
    phi_lambda[j,] = eta_lambda * phi_lambda[j-1,] + (1/tau_lambda) * phi_init_lambda[j]';
    phi_pi[j,] = eta_pi * phi_pi[j-1,] + (1/tau_pi) * phi_init_pi[j]';
  }
  
  for (i in 1:r) {
    lambda[, i] = X_tc[i] * beta_lambda[, i] + phi_lambda[idx_tc_er, i] + area_offset[i];
    pi_param[, i] = X_tc[i] * beta_pi[, i] + phi_pi[idx_tc_er, i];
  }
}

model {
  // priors
  bp_init_lambda ~ uniform(0, 1);
  bp_init_pi ~ uniform(0, 1);
  
  rho1_lambda ~ beta(1.5, 4);
  rho2_lambda ~ beta(3, 4);
  rho1_pi ~ beta(1.5, 4);
  rho2_pi ~ beta(3, 4);
  
  target += matnormal_lpdf(beta_lambda | cov_ar1_lambda, corr_lambda);
  target += matnormal_lpdf(beta_pi | cov_ar1_pi, corr_pi);
  
  // IAR prior
  tau_init_lambda ~ exponential(1);
  tau_init_pi ~ exponential(1);
  eta_lambda ~ beta(2,8);
  eta_pi ~ beta(2,8);
  for (j in 1:t_tc) {
    target += -.5 * dot_self(phi_init_lambda[j][node1] - phi_init_lambda[j][node2]);
    sum(phi_init_lambda[j]) ~ normal(0, 0.001*r);
    target += -.5 * dot_self(phi_init_pi[j][node1] - phi_init_pi[j][node2]);
    sum(phi_init_pi[j]) ~ normal(0, 0.001*r);
  }

  // likelihood
  for (i in 1:r) {
    for (j in 1:t_tc) {
       if (y_tc[j, i] == 0) {
        target += log_sum_exp(bernoulli_logit_lpmf(1 | pi_param[j, i]),
                            bernoulli_logit_lpmf(0 | pi_param[j, i])
                            + poisson_log_lpmf(y_tc[j, i] | lambda[j, i]));
      } else {
        target += bernoulli_logit_lpmf(0 | pi_param[j, i])
                + poisson_log_lpmf(y_tc[j, i] | lambda[j, i]);
      }
    }
  }
}

generated quantities {
  matrix[t_all, r] lambda_full;  
  matrix[t_all, r] pi_param_full; 
  matrix[t_ho, r] lambda_hold;  
  matrix[t_ho, r] pi_param_hold;
  
  matrix[t_ho, r] holdout_loglik;
  matrix[t_tc, r] train_loglik;
 
  // expected values of parameters based on all timepoints, then cut to only be holdout parameters
  for (i in 1:r) {
    lambda_full[, i] = X_all_tmpt[i] * beta_lambda[, i] + phi_lambda[, i] + area_offset[i];
    pi_param_full[, i] = X_all_tmpt[i] * beta_pi[, i] + phi_pi[, i];
  }
  lambda_hold = lambda_full[idx_hold_er, ];
  pi_param_hold = pi_param_full[idx_hold_er, ];
  
  // training log-likelihood
  for (i in 1:r) {
    for (j in 1:t_tc) {
       if (y_tc[j, i] == 0) {
        train_loglik[j, i] = log_sum_exp(bernoulli_logit_lpmf(1 | pi_param[j, i]),
                             bernoulli_logit_lpmf(0 | pi_param[j, i])
                             + poisson_log_lpmf(y_tc[j, i] | lambda[j, i]));
      } else {
        train_loglik[j, i] = bernoulli_logit_lpmf(0 | pi_param[j, i]) 
                             + poisson_log_lpmf(y_tc[j, i] | lambda[j, i]);
      }
    }
  }
  
  // holdout log-likelihood
  for (i in 1:r) {
    for (j in 1:t_ho) {
       if (y_hold[j, i] == 0) {
        holdout_loglik[j, i] = log_sum_exp(bernoulli_logit_lpmf(1 | pi_param_hold[j, i]),
                             bernoulli_logit_lpmf(0 | pi_param_hold[j, i])
                             + poisson_log_lpmf(y_hold[j, i] | lambda_hold[j, i]));
      } else {
        holdout_loglik[j, i] = bernoulli_logit_lpmf(0 | pi_param_hold[j, i]) 
                             + poisson_log_lpmf(y_hold[j, i] | lambda_hold[j, i]);
      }
    }
  }
}
