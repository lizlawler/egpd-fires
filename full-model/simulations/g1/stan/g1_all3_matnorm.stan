functions {
  real egpd_g1_lpdf(real y, real sigma, real xi, real kappa) {
    real lpdf;
    lpdf = log(kappa) - log(sigma) - (1/xi + 1) * log(1 + xi * (y/sigma)) + 
    (kappa-1) * log(1 - (1 + xi * (y/sigma))^(-1/xi));
    return lpdf;
  }
  
  real matnormal_lpdf(matrix y, matrix u, matrix v) {
    real lpdf;
    real n;
    real p;
    n = rows(u);
    p = rows(v);
    lpdf = -(n*p/2) * log(2 * pi()) - (n/2)*log_determinant(v) - (p/2)*log_determinant(u) -
          0.5 * trace(mdivide_right_spd(mdivide_left_spd(v, y'), u) * y);
    return lpdf;
  }
}

data {
  int<lower = 1> p; // # of parameters
  int<lower = 1> T; // # of timepoints
  int<lower = 1> R; // # of regions
  int<lower = 1> N_obs;
  // int<lower = 1> N_miss;
  // int<lower = 1, upper = N_obs + N_miss> N;
  int<lower = 1> ii_obs[N_obs];
  // int<lower = 1, upper = N> ii_miss[N_miss];
  // int<lower = 1, upper = N> ii_full[N];
  int<lower=0> n_edges;
  int<lower=1, upper = R> node1[n_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper = R> node2[n_edges];  // and node1[i] < node2[i]
  
  matrix[T, p] X[R]; // design matrix; 1-D array of size r with matrices t x p; this indexing is deprecated in version 2.32 and above
  real<lower = 0> y_obs[N_obs]; // observed response data
  
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

// should consider transforming X and y so they're centered and scaled

parameters {
  // real y_miss[N_miss];
  vector[R] phi_init_kappa[T];
  vector[R] phi_init_nu[T];
  vector[R] phi_init_xi[T];
  matrix[R, p] beta_kappa;
  matrix[R, p] beta_nu;
  matrix[R, p] beta_xi;
  real<lower = 0> tau_init_kappa;
  real<lower = 0> tau_init_nu;
  real<lower = 0> tau_init_xi;
  real<lower = 0, upper = 1> eta_kappa;
  real<lower = 0, upper = 1> eta_nu;
  real<lower = 0, upper = 1> eta_xi;
  real<lower = 0, upper = 1> bp_init_kappa;
  real<lower = 0, upper = 1> bp_init_nu;
  real<lower = 0, upper = 1> bp_init_xi;
  real<lower = 0, upper = 1> rho2_kappa;
  real<lower=0, upper = (1-rho2_kappa)> rho1_kappa;
  real<lower = 0, upper = 1> rho2_nu;
  real<lower=0, upper = (1-rho2_nu)> rho1_nu;
  real<lower = 0, upper = 1> rho2_xi;
  real<lower=0, upper = (1-rho2_xi)> rho1_xi;
}

transformed parameters {
  // real y[N];
  matrix[T, R] phi_kappa;
  matrix[T, R] phi_nu;
  matrix[T, R] phi_xi;
  matrix[T, R] reg_kappa;
  matrix[T, R] reg_nu;
  matrix[T, R] reg_xi;

  real<lower=0, upper = bp_init_kappa/2> bp_kappa = bp_init_kappa/2;
  real<lower=0, upper = bp_init_nu/2> bp_nu = bp_init_nu/2;
  real<lower=0, upper = bp_init_xi/2> bp_xi = bp_init_xi/2;
  real<lower=0, upper = tau_init_kappa/2> tau_kappa = tau_init_kappa/2;
  real<lower=0, upper = tau_init_nu/2> tau_nu = tau_init_nu/2;
  real<lower=0, upper = tau_init_xi/2> tau_xi = tau_init_xi/2;

  matrix[R, R] corr_kappa = l3 + rho2_kappa * l2 + rho1_kappa * l1;
  matrix[p, p] cov_ar1_kappa = equal + bp_kappa * bp_lin + bp_kappa^2 * bp_square + bp_kappa^3 * bp_cube + bp_kappa^4 * bp_quart;

  matrix[R, R] corr_nu = l3 + rho2_nu * l2 + rho1_nu * l1;
  matrix[p, p] cov_ar1_nu = equal + bp_nu * bp_lin + bp_nu^2 * bp_square + bp_nu^3 * bp_cube + bp_nu^4 * bp_quart;
  
  matrix[R, R] corr_xi = l3 + rho2_xi * l2 + rho1_xi * l1;
  matrix[p, p] cov_ar1_xi = equal + bp_xi * bp_lin + bp_xi^2 * bp_square + bp_xi^3 * bp_cube + bp_xi^4 * bp_quart;

  vector[R*T] kappa;
  vector[R*T] nu;
  vector[R*T] xi;
  vector[R*T] sigma;

  // y[ii_obs] = y_obs;
  // y[ii_miss] = y_miss;
  
  phi_kappa[1,] = (1/tau_kappa) * phi_init_kappa[1]';
  phi_nu[1,] = (1/tau_nu) * phi_init_nu[1]';
  phi_xi[1,] = (1/tau_xi) * phi_init_xi[1]';
  for (j in 2:T) {
    phi_kappa[j,] = eta_kappa * phi_kappa[j-1,] + (1/tau_kappa) * phi_init_kappa[j]';
    phi_nu[j,] = eta_nu * phi_nu[j-1,] + (1/tau_nu) * phi_init_nu[j]';
    phi_xi[j,] = eta_xi * phi_xi[j-1,] + (1/tau_xi) * phi_init_xi[j]';
  }
  
  for (i in 1:R) {
    reg_kappa[, i] = X[i] * beta_kappa[i, ]'/4 + phi_kappa[, i]/10;
    reg_nu[, i] = X[i] * beta_nu[i, ]'/4 + phi_nu[, i]/10;
    reg_xi[, i] = X[i] * beta_xi[i, ]'/8 + phi_xi[, i]/15;
  }

  kappa = to_vector(exp(reg_kappa));
  nu = to_vector(exp(reg_nu));
  xi = to_vector(exp(reg_xi));
  sigma = nu ./ (1 + xi);
}

model {
  vector[N_obs] kappa_trunc = kappa[ii_obs];
  vector[N_obs] nu_trunc = nu[ii_obs];
  vector[N_obs] xi_trunc = xi[ii_obs];
  vector[N_obs] sigma_trunc = sigma[ii_obs];
  
  // priors
  bp_init_kappa ~ uniform(0, 1);
  bp_init_nu ~ uniform(0, 1);
  bp_init_xi ~ uniform(0, 1);
  
  rho1_kappa ~ beta(1.5, 4);
  rho2_kappa ~ beta(3, 4);
  rho1_nu ~ beta(1.5, 4);
  rho2_nu ~ beta(3, 4);
  rho1_xi ~ beta(1.5, 4);
  rho2_xi ~ beta(3, 4);
  
  target += matnormal_lpdf(beta_kappa | corr_kappa, cov_ar1_kappa);
  target += matnormal_lpdf(beta_nu | corr_nu, cov_ar1_nu);
  target += matnormal_lpdf(beta_xi | corr_xi, cov_ar1_xi);
  
  // IAR prior
  eta_kappa ~ beta(2,8);
  eta_nu ~ beta(2,8);
  eta_xi ~ beta(2,8);
  tau_init_kappa ~ exponential(1);
  tau_init_nu ~ exponential(1);
  tau_init_xi ~ exponential(1);
  for (j in 1:T) {
    target += -.5 * dot_self(phi_init_kappa[j][node1] - phi_init_kappa[j][node2]);
    sum(phi_init_kappa[j]) ~ normal(0, 0.001*R);
    target += -.5 * dot_self(phi_init_nu[j][node1] - phi_init_nu[j][node2]);
    sum(phi_init_nu[j]) ~ normal(0, 0.001*R);
    target += -.5 * dot_self(phi_init_xi[j][node1] - phi_init_xi[j][node2]);
    sum(phi_init_xi[j]) ~ normal(0, 0.001*R);
  }
  // 
  // likelihood
  for (k in 1:N_obs) {
    target += egpd_g1_lpdf(y_obs[k] | sigma_trunc[k], xi_trunc[k], kappa_trunc[k]);
  }
}
