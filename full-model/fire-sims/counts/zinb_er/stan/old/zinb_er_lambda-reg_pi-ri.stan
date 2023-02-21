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

transformed data {
  int S = 1; // of parameters with regression (either 1 or 2)
  int C = 3; // # of parameters with correlation (includes random intercept only)
}

parameters {
  vector[r] Z[2];
  row_vector[r] phi_init[t_all, S];
  matrix[p, r] beta[S];
  real<lower = 0> tau_init[S];
  real<lower = 0, upper = 1> eta[S];
  real<lower = 0, upper = 1> bp_init[S];
  real<lower = 0, upper = 1> rho[2, C]; // ordering: 1 = lambda, 2 = pi, 3 = delta
}

transformed parameters {
  vector<lower = 0>[r] delta;
  matrix[t_all, r] phi[S];
  matrix[t_train, r] reg[S];
  matrix[p, p] cov_ar1[S];
  real<lower = 0, upper = 1> bp[S];
  real<lower = 0> tau[S];
  matrix[r, r] corr[C];
  matrix[t_train, r] lambda;
  vector[r] pi_prob;
  
  for (i in 1:C) {
    corr[i] = l3 + rho[2, i] * l2 + rho[1, i] * l1;
  }
  
  for (i in 1:S) {
    bp[i] = bp_init[i]/2;
    tau[i] = tau_init[i]/2;
    cov_ar1[i] = equal + bp[i] * bp_lin + bp[i]^2 * bp_square + bp[i]^3 * bp_cube + bp[i]^4 * bp_quart;
    
    // ICAR variables
    phi[i][1, ] = 1/tau[i] * phi_init[1, i];
    for (j in 2:t_all) {
      phi[i][j, ] = eta[i] * phi[i][j-1, ] + 1/tau[i] * phi_init[j, i];
    }
    
    // regression for lambda (and pi if necessary)
    for (k in 1:r) {
      reg[i][, k] = X_train[k] * beta[i][, k] + phi[i][idx_train_er, k];
      lambda[,k] = reg[1][, k] + area_offset[k];
    }
  }
  pi_prob = exp(cholesky_decompose(corr[2])' * Z[1]);
  delta = exp(cholesky_decompose(corr[3])' * Z[2]);
}

model {
  Z[1] ~ normal(0, 1);
  Z[2] ~ normal(0, 1);
  // priors on rhos and AR(1) penalization of splines
  to_vector(bp_init) ~ uniform(0, 1);
  to_vector(rho[1,]) ~ beta(3, 4); // prior on rho1 for lambda, pi, and delta
  // to_vector(rho[2,]) ~ beta(1.5, 4); // prior on rho2 for lambda, pi, and delta
  
  // priors scaling constants in ICAR
  to_vector(eta) ~ beta(2,8);
  to_vector(tau_init) ~ exponential(1);
  
  for (i in 1:C) {
    // soft constraint for sum of rhos within an individual param to be <= 1 (ie rho1kappa + rho2kappa <= 1)
    sum(rho[,i]) ~ uniform(0,1); 
  }
  
  for (i in 1:S) {
    // MVN prior on betas
    target += matnormal_lpdf(beta[i] | cov_ar1[i], corr[i]);
    // ICAR prior
    for (j in 1:t_all) {
      target += -.5 * dot_self(phi_init[j, i][node1] - phi_init[j, i][node2]);
      sum(phi_init[j, i]) ~ normal(0, 0.001*r);
    }
  }

  // likelihood
  for (i in 1:r) {
    for (j in 1:t_train) {
       if (y_train_count[j, i] == 0) {
        target += log_sum_exp(bernoulli_logit_lpmf(1 | pi_prob[i]),
                            bernoulli_logit_lpmf(0 | pi_prob[i])
                            + neg_binomial_2_log_lpmf(y_train_count[j, i] | lambda[j, i], delta[i]));
      } else {
        target += bernoulli_logit_lpmf(0 | pi_prob[i])
                + neg_binomial_2_log_lpmf(y_train_count[j, i] | lambda[j, i], delta[i]);
      }
    }
  }
}

generated quantities {
  matrix[t_all, r] reg_full[S];  
  matrix[t_hold, r] lambda_hold;
  
  matrix[t_hold, r] holdout_loglik;
  matrix[t_train, r] train_loglik;
 
  // expected values of parameters based on all timepoints, then cut to only be holdout parameters
  for (i in 1:S) {
    for (k in 1:r) {
      reg_full[i][, k] = X_full[k] * beta[i][, k] + phi[i][, k];
      lambda_hold[,k] = reg_full[1][idx_hold_er, k] + area_offset[k];
    }
  }
  
  // training log-likelihood
  for (i in 1:r) {
    for (j in 1:t_train) {
       if (y_train_count[j, i] == 0) {
        train_loglik[j, i] = log_sum_exp(bernoulli_logit_lpmf(1 | pi_prob[i]),
                             bernoulli_logit_lpmf(0 | pi_prob[i])
                             + neg_binomial_2_log_lpmf(y_train_count[j, i] | lambda[j, i], delta[i]));
      } else {
        train_loglik[j, i] = bernoulli_logit_lpmf(0 | pi_prob[i]) 
                             + neg_binomial_2_log_lpmf(y_train_count[j, i] | lambda[j, i], delta[i]);
      }
    }
  }
  
  // holdout log-likelihood
  for (i in 1:r) {
    for (j in 1:t_hold) {
       if (y_hold_count[j, i] == 0) {
        holdout_loglik[j, i] = log_sum_exp(bernoulli_logit_lpmf(1 | pi_prob[i]),
                             bernoulli_logit_lpmf(0 | pi_prob[i])
                             + neg_binomial_2_log_lpmf(y_hold_count[j, i] | lambda_hold[j, i], delta[i]));
      } else {
        holdout_loglik[j, i] = bernoulli_logit_lpmf(0 | pi_prob[i]) 
                             + neg_binomial_2_log_lpmf(y_hold_count[j, i] | lambda_hold[j, i], delta[i]);
      }
    }
  }
}
