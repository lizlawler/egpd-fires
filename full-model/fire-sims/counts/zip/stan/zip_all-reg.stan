functions {
  real matnormal_lpdf(matrix y, matrix cov, matrix corr) {
    real lpdf;
    real r;
    real p;
    r = rows(corr);
    p = rows(cov);
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
  array[R] matrix[T_all, p] X_full; // design matrix; 1-D array of size r with matrices t x p; this indexing is deprecated in version 2.32 and above
  array[R] matrix[T_train, p] X_train; // design matrix; 1-D array of size r with matrices t x p; this indexing is deprecated in version 2.32 and above
  
  // area offset
  array[R] real area_offset; // known offset vector of areas
  
  // training dataset
  array[T_train] int<lower=1> idx_train_er; // vector of indices for training data timepoints
  array[T_train, R] int<lower=0> y_train_count; // response data
  
  // holdout dataset
  array[T_hold] int<lower=1> idx_hold_er; // vector of indices for holdout data timepoints
  array[T_hold, R] int<lower=0> y_hold_count; // response data
  
  // neighbor information
  int<lower=0> n_edges;
  array[n_edges] int<lower=1, upper=R> node1; // node1[i] adjacent to node2[i]
  array[n_edges] int<lower=1, upper=R> node2; // and node1[i] < node2[i]
  
  // indicator matrices for ecoregions
  matrix[R, R] l3;
  matrix[R, R] l2;
  matrix[R, R] l1;
  
  // indicator matrices for AR(1) process on betas
  matrix[p, p] equal;
  matrix[p, p] bp_lin;
  matrix[p, p] bp_square;
  matrix[p, p] bp_cube;
  matrix[p, p] bp_quart;
}
transformed data {
  int S = 2; // # of parameters with regression (either 1 or 2)
  int C = 2; // # of parameters with correlation (includes random intercept only)
}
parameters {
  array[T_all, S] row_vector[R] phi_init;
  array[S] matrix[p, R] beta;
  array[S] real<lower=0> tau_init;
  array[S] real<lower=0, upper=1> eta;
  array[S] real<lower=0, upper=1> bp_init;
  array[2, C] real<lower=0, upper=1> rho; // ordering: 1 = lambda, 2 = pi
}
transformed parameters {
  array[S] matrix[T_all, R] phi;
  array[S] matrix[T_train, R] reg;
  array[S] matrix[p, p] cov_ar1;
  array[S] real<lower=0, upper=1> bp;
  array[S] real<lower=0> tau;
  array[C] matrix[R, R] corr;
  matrix[T_train, R] lambda;
  matrix[T_train, R] pi_prob;
  
  for (c in 1:C) {
    corr[c] = l3 + rho[2, c] * l2 + rho[1, c] * l1;
  }
  
  for (s in 1:S) {
    bp[s] = bp_init[s] / 2;
    tau[s] = tau_init[s] / 2;
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
  pi_prob = reg[2];
}
model {
  // priors on rhos and AR(1) penalization of splines
  to_vector(bp_init) ~ uniform(0, 1);
  to_vector(rho[1, ]) ~ beta(3, 4); // prior on rho1 for lambda and pi
  // to_vector(rho[2,]) ~ beta(1.5, 4); // prior on rho2 for lambda and pi
  
  // priors scaling constants in ICAR
  to_vector(eta) ~ beta(2, 8);
  to_vector(tau_init) ~ exponential(1);
  
  for (c in 1:C) {
    // soft constraint for sum of rhos within an individual param to be <= 1 (ie rho1kappa + rho2kappa <= 1)
    sum(rho[, c]) ~ uniform(0, 1);
  }
  
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
        target += log_sum_exp(bernoulli_logit_lpmf(1 | pi_prob[t, r]),
                              bernoulli_logit_lpmf(0 | pi_prob[t, r])
                              + poisson_log_lpmf(y_train_count[t, r] | lambda[t, r]));
      } else {
        target += bernoulli_logit_lpmf(0 | pi_prob[t, r])
                  + poisson_log_lpmf(y_train_count[t, r] | lambda[t, r]);
      }
    }
  }
}
generated quantities {
  array[S] matrix[T_all, R] reg_full;
  matrix[T_hold, R] lambda_hold;
  matrix[T_hold, R] pi_prob_hold;
  
  matrix[T_hold, R] holdout_loglik;
  matrix[T_train, R] train_loglik;
  
  // expected values of parameters based on all timepoints, then cut to only be holdout parameters
  for (s in 1:S) {
    for (r in 1:R) {
      reg_full[s][, r] = X_full[r] * beta[s][, r] + phi[s][, r];
      lambda_hold[, r] = reg_full[1][idx_hold_er, r] + area_offset[r];
    }
  }
  pi_prob_hold = reg_full[2][idx_hold_er, ];
  
  // training log-likelihood
  for (r in 1:R) {
    for (t in 1:T_train) {
      if (y_train_count[t, r] == 0) {
        train_loglik[t, r] = log_sum_exp(bernoulli_logit_lpmf(1 | pi_prob[t, r]),
                                         bernoulli_logit_lpmf(0 | pi_prob[t, r])
                                         + poisson_log_lpmf(y_train_count[t, r] | lambda[t, r]));
      } else {
        train_loglik[t, r] = bernoulli_logit_lpmf(0 | pi_prob[t, r])
                             + poisson_log_lpmf(y_train_count[t, r] | lambda[t, r]);
      }
    }
  }
  
  // holdout log-likelihood
  for (r in 1:R) {
    for (t in 1:T_hold) {
      if (y_hold_count[t, r] == 0) {
        holdout_loglik[t, r] = log_sum_exp(bernoulli_logit_lpmf(1 | pi_prob_hold[t, r]),
                                           bernoulli_logit_lpmf(0 | pi_prob_hold[t, r])
                                           + poisson_log_lpmf(y_hold_count[t, r] | lambda_hold[t, r]));
      } else {
        holdout_loglik[t, r] = bernoulli_logit_lpmf(0 | pi_prob_hold[t, r])
                               + poisson_log_lpmf(y_hold_count[t, r] | lambda_hold[t, r]);
      }
    }
  }
}


