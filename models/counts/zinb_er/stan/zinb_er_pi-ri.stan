functions {
  #include /../../../twcrps_matnorm_fcns.stanfunctions
}
#include /../../counts_data.stan
transformed data {
  int S = 1; // # of parameters with regression (either 1 or 2)
  int C = 3; // # of parameters with correlation (includes random intercept only)
}
parameters {
  matrix[R, 2] Z;
  array[T_all, S] row_vector[R] phi_init;
  array[S] matrix[p, R] beta;
  vector<lower=0>[S] tau_init;
  vector<lower=0, upper = 1>[S] eta;
  vector<lower=0, upper = 1>[S] bp_init;
  vector<lower=0, upper = 1>[C] rho1;
  vector<lower=rho1, upper = 1>[C] rho_sum; // ordering: 1 = lambda, 2 = pi, 3 = delta
}
transformed parameters {
  vector<lower = 0>[R] delta;
  array[S] matrix[T_all, R] phi;
  array[S] matrix[T_train, R] reg;
  vector<lower=0>[S] bp = bp_init / 2;
  vector<lower=0>[S] tau = tau_init / 2;
  vector[C] rho2 = rho_sum - rho1;
  array[S] cov_matrix[p] cov_ar1;
  array[C] corr_matrix[R] corr;
  
  matrix[T_train, R] lambda;
  vector[R] pi_prob;
  
  for (c in 1:C) {
    corr[c] = l3 + rho2[c] * l2 + rho1[c] * l1;
  }
  
  for (s in 1:S) {
    cov_ar1[s] = equal + bp[s] * bp_lin + bp[s] ^ 2 * bp_square
                 + bp[s] ^ 3 * bp_cube + bp[s] ^ 4 * bp_quart;
    
    // ICAR variables
    phi[s][1, ]= 1 / tau[s] * phi_init[1, s];
    for (t in 2:T_all) {
      phi[s][t, ]= eta[s] * phi[s][t - 1, ]
                       + 1 / tau[s] * phi_init[t, s];
    }
    
    // regression for lambda (and pi if necessary)
    for (r in 1:R) {
      reg[s][, r] = X_train[r] * beta[s][, r] + phi[s][idx_train_er, r];
      lambda[, r] = reg[1][, r] + area_offset[r];
    }
  }
  pi_prob = cholesky_decompose(corr[2])' * Z[,1];
  delta = exp(cholesky_decompose(corr[3])' * Z[,2]);
}

model {
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
  array[S] matrix[T_all, R] reg_full;
  matrix[T_hold, R] lambda_hold;
  
  matrix[T_hold, R] holdout_loglik;
  matrix[T_train, R] train_loglik;
  
  // expected values of parameters based on all timepoints, then cut to only be holdout parameters
  for (s in 1:S) {
    for (r in 1:R) {
      reg_full[s][, r] = X_full[r] * beta[s][, r] + phi[s][, r];
      lambda_hold[, r] = reg_full[1][idx_hold_er, r] + area_offset[r];
    }
  }
  
  // training log-likelihood
  for (r in 1:R) {
    for (t in 1:T_train) {
      if (y_train_count[t, r] == 0) {
        train_loglik[t, r] = log_sum_exp(bernoulli_logit_lpmf(1 | pi_prob[r]),
                                         bernoulli_logit_lpmf(0 | pi_prob[r])
                                         + neg_binomial_2_log_lpmf(y_train_count[t, r] | lambda[t, r], delta[r]));
      } else {
        train_loglik[t, r] = bernoulli_logit_lpmf(0 | pi_prob[r])
                             + neg_binomial_2_log_lpmf(y_train_count[t, r] | lambda[t, r], delta[r]);
      }
    }
  }
  
  // holdout log-likelihood
  for (r in 1:R) {
    for (t in 1:T_hold) {
      if (y_hold_count[t, r] == 0) {
        holdout_loglik[t, r] = log_sum_exp(bernoulli_logit_lpmf(1 | pi_prob[r]),
                                           bernoulli_logit_lpmf(0 | pi_prob[r])
                                           + neg_binomial_2_log_lpmf(y_hold_count[t, r] | lambda_hold[t, r], delta[r]));
      } else {
        holdout_loglik[t, r] = bernoulli_logit_lpmf(0 | pi_prob[r])
                               + neg_binomial_2_log_lpmf(y_hold_count[t, r] | lambda_hold[t, r], delta[r]);
      }
    }
  }
}
