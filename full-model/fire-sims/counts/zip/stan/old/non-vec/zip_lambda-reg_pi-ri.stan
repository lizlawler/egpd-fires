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
  int<lower = 1> r; // # of regions
  int<lower = 1> p; // # of parameters
  int<lower = 1> t_all; // # of timepoints in full dataset
  int<lower = 1> t_train;
  int<lower = 1> t_hold;
  
  // covariate data
  matrix[t_all, p] X_full[r]; // design matrix; 1-D array of size r with matrices t x p; this indexing is deprecated in version 2.32 and above
  matrix[t_train, p] X_train[r]; // design matrix; 1-D array of size r with matrices t x p; this indexing is deprecated in version 2.32 and above
  
  // area offset
  real area_offset[r]; // known offset vector of areas
  
  // training dataset
  int<lower = 1> idx_train_er[t_train]; // vector of indices for training data timepoints
  int<lower = 0> y_train_count[t_train, r]; // response data
  
  // holdout dataset
  int<lower = 1> idx_hold_er[t_hold]; // vector of indices for holdout data timepoints
  int<lower = 0> y_hold_count[t_hold, r]; // response data

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
  vector[r] Z_pi;
  vector[r] phi_init_lambda[t_all];
  matrix[p, r] beta_lambda;
  real<lower = 0> tau_init_lambda;
  real<lower = 0, upper = 1> eta_lambda;
  real<lower = 0, upper = 1> bp_init_lambda;
  real<lower = 0, upper = 1> rho2_lambda;
  real<lower=0, upper = (1-rho2_lambda)> rho1_lambda;
  real<lower = 0, upper = 1> rho2_pi;
  real<lower=0, upper = (1-rho2_pi)> rho1_pi;
}

transformed parameters {
  matrix[t_all, r] phi_lambda;
  matrix[t_train, r] lambda;
  vector[r] pi_param;
  
  real<lower=0, upper = bp_init_lambda/2> bp_lambda = bp_init_lambda/2;
  real<lower=0, upper = tau_init_lambda/2> tau_lambda = tau_init_lambda/2;
  
  matrix[r, r] corr_lambda = l3 + rho2_lambda * l2 + rho1_lambda * l1;
  matrix[p, p] cov_ar1_lambda = equal + bp_lambda * bp_lin + bp_lambda^2 * bp_square + bp_lambda^3 * bp_cube + bp_lambda^4 * bp_quart;
  matrix[r, r] corr_pi = l3 + rho2_pi * l2 + rho1_pi * l1;

  phi_lambda[1,] = (1/tau_lambda) * phi_init_lambda[1]';
  for (j in 2:t_all) {
    phi_lambda[j,] = eta_lambda * phi_lambda[j-1,] + (1/tau_lambda) * phi_init_lambda[j]';
  }
  
  for (i in 1:r) {
    lambda[, i] = X_train[i] * beta_lambda[, i] + phi_lambda[idx_train_er, i] + area_offset[i];
  }
  pi_param = cholesky_decompose(corr_pi)' * Z_pi;
}

model {
  // priors
  bp_init_lambda ~ uniform(0, 1);
  Z_pi ~ normal(0, 1);
  
  rho1_lambda ~ beta(1.5, 4);
  rho2_lambda ~ beta(3, 4);
  rho1_pi ~ beta(1.5, 4);
  rho2_pi ~ beta(3, 4);
  
  target += matnormal_lpdf(beta_lambda | cov_ar1_lambda, corr_lambda);
  
  // IAR prior
  tau_init_lambda ~ exponential(1);
  eta_lambda ~ beta(2,8);
  for (j in 1:t_all) {
    target += -.5 * dot_self(phi_init_lambda[j][node1] - phi_init_lambda[j][node2]);
    sum(phi_init_lambda[j]) ~ normal(0, 0.001*r);
  }

  // likelihood
  for (i in 1:r) {
    for (j in 1:t_train) {
       if (y_train_count[j, i] == 0) {
        target += log_sum_exp(bernoulli_logit_lpmf(1 | pi_param[i]),
                            bernoulli_logit_lpmf(0 | pi_param[i])
                            + poisson_log_lpmf(y_train_count[j, i] | lambda[j, i]));
      } else {
        target += bernoulli_logit_lpmf(0 | pi_param[i])
                + poisson_log_lpmf(y_train_count[j, i] | lambda[j, i]);
      }
    }
  }
}

generated quantities {
  matrix[t_all, r] lambda_full;  
  matrix[t_hold, r] lambda_hold;  
  
  matrix[t_hold, r] holdout_loglik;
  matrix[t_train, r] train_loglik;
 
  // expected values of parameters based on all timepoints, then cut to only be holdout parameters
  for (i in 1:r) {
    lambda_full[, i] = X_full[i] * beta_lambda[, i] + phi_lambda[, i] + area_offset[i];
  }
  lambda_hold = lambda_full[idx_hold_er, ];
  
  // training log-likelihood
  for (i in 1:r) {
    for (j in 1:t_train) {
       if (y_train_count[j, i] == 0) {
        train_loglik[j, i] = log_sum_exp(bernoulli_logit_lpmf(1 | pi_param[i]),
                             bernoulli_logit_lpmf(0 | pi_param[i])
                             + poisson_log_lpmf(y_train_count[j, i] | lambda[j, i]));
      } else {
        train_loglik[j, i] = bernoulli_logit_lpmf(0 | pi_param[i]) 
                             + poisson_log_lpmf(y_train_count[j, i] | lambda[j, i]);
      }
    }
  }
  
  // holdout log-likelihood
  for (i in 1:r) {
    for (j in 1:t_hold) {
       if (y_hold_count[j, i] == 0) {
        holdout_loglik[j, i] = log_sum_exp(bernoulli_logit_lpmf(1 | pi_param[i]),
                             bernoulli_logit_lpmf(0 | pi_param[i])
                             + poisson_log_lpmf(y_hold_count[j, i] | lambda_hold[j, i]));
      } else {
        holdout_loglik[j, i] = bernoulli_logit_lpmf(0 | pi_param[i]) 
                             + poisson_log_lpmf(y_hold_count[j, i] | lambda_hold[j, i]);
      }
    }
  }
}
