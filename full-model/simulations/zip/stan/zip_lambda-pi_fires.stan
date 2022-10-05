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
  int<lower = 1> T; // # of timepoints
  int<lower = 1> R; // # of regions
  int<lower = 1> N; // total # of observations (TxR)
  int<lower=0> n_edges;
  int<lower=1, upper = R> node1[n_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper = R> node2[n_edges];  // and node1[i] < node2[i]
  
  matrix[T, p] X[R]; // design matrix; 1-D array of size r with matrices t x p; this indexing is deprecated in version 2.32 and above
  int<lower = 0> y[N]; // response data
  
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

parameters {
  vector[R] phi_init_lambda[T];
  vector[R] phi_init_pi[T];
  matrix[p, R] beta_lambda;
  matrix[p, R] beta_pi;
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
  // real<lower = 0, upper = 1> pi_zip_init;
}

transformed parameters {
  matrix[T, R] phi_lambda;
  matrix[T, R] phi_pi;
  matrix[T, R] reg_lambda;
  matrix[T, R] reg_pi;
  vector<lower = 0> [N] lambda;
  vector<lower = 0, upper = 1> [N] pi_param;
  
  real<lower=0, upper = bp_init_lambda/2> bp_lambda = bp_init_lambda/2;
  real<lower=0, upper = tau_init_lambda/2> tau_lambda = tau_init_lambda/2;
  real<lower=0, upper = bp_init_pi/2> bp_pi = bp_init_pi/2;
  real<lower=0, upper = tau_init_pi/2> tau_pi = tau_init_pi/2;
  
  matrix[R, R] corr_lambda = l3 + rho2_lambda * l2 + rho1_lambda * l1;
  matrix[p, p] cov_ar1_lambda = equal + bp_lambda * bp_lin + bp_lambda^2 * bp_square + bp_lambda^3 * bp_cube + bp_lambda^4 * bp_quart;
  matrix[R, R] corr_pi = l3 + rho2_pi * l2 + rho1_pi * l1;
  matrix[p, p] cov_ar1_pi = equal + bp_pi * bp_lin + bp_pi^2 * bp_square + bp_pi^3 * bp_cube + bp_pi^4 * bp_quart;

  phi_lambda[1,] = (1/tau_lambda) * phi_init_lambda[1]';
  phi_pi[1,] = (1/tau_pi) * phi_init_pi[1]';
  for (j in 2:T) {
    phi_lambda[j,] = eta_lambda * phi_lambda[j-1,] + (1/tau_lambda) * phi_init_lambda[j]';
    phi_pi[j,] = eta_pi * phi_pi[j-1,] + (1/tau_pi) * phi_init_pi[j]';
  }
  
  for (i in 1:R) {
    reg_lambda[, i] = X[i] * beta_lambda[, i] + phi_lambda[, i];
    reg_pi[, i] = X[i] * beta_pi[, i] + phi_pi[, i];
  }

  lambda = to_vector(exp(reg_lambda));
  pi_param = to_vector(inv_logit(reg_pi));
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
  for (j in 1:T) {
    target += -.5 * dot_self(phi_init_lambda[j][node1] - phi_init_lambda[j][node2]);
    sum(phi_init_lambda[j]) ~ normal(0, 0.001*R);
    target += -.5 * dot_self(phi_init_pi[j][node1] - phi_init_pi[j][node2]);
    sum(phi_init_pi[j]) ~ normal(0, 0.001*R);
  }

  // likelihood
  for (i in 1:N) {
    if (y[i] == 0) {
      target += log_sum_exp(bernoulli_lpmf(1 | pi_param[i]),
                            bernoulli_lpmf(0 | pi_param[i])
                            + poisson_lpmf(y[i] | lambda[i]));
    } else {
      target += bernoulli_lpmf(0 | pi_param[i])
                + poisson_lpmf(y[i] | lambda[i]);
    }
  }
}
