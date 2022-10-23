functions {
  real egpd_g1_lpdf(real y, real sigma, real xi, real kappa) {
    real lpdf;
    lpdf = log(kappa) - log(sigma) - (1/xi + 1) * log(1 + xi * (y/sigma)) + 
    (kappa-1) * log(1 - (1 + xi * (y/sigma))^(-1/xi));
    return lpdf;
  }
  
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
  int<lower = 1> t_tb; // # of timepoints
  int<lower = 1> r; // # of regions
  int<lower = 1> N_obs;
  int<lower = 1> N_mis;
  int<lower = 1> N_all;
  int<lower = 1> ii_obs[N_obs];
  int<lower = 1, upper = N_all> ii_mis[N_mis];
  int<lower = 1, upper = N_all> ii_all[N_all]; // for broadcasting

  int<lower=0> n_edges;
  int<lower=1, upper = r> node1[n_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper = r> node2[n_edges];  // and node1[i] < node2[i]
  
  matrix[t_tb, p] X[r]; // design matrix; 1-D array of size r with matrices t x p; this indexing is deprecated in version 2.32 and above
  real<lower = 0> y_obs[N_obs]; // observed response data
  
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
  real<lower = 0> y_mis[N_mis];
  // remember to change to ALL timepoints for the phi matrices
  vector[r] phi_init_kappa[t_tb];
  vector[r] phi_init_nu[t_tb];
  vector[r] phi_init_xi[t_tb];
  matrix[p, r] beta_kappa;
  matrix[p, r] beta_nu;
  matrix[p, r] beta_xi;
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
  real<lower = 0> y[N_all];
  matrix[t_tb, r] phi_kappa;
  matrix[t_tb, r] phi_nu;
  matrix[t_tb, r] phi_xi;
  matrix[t_tb, r] reg_kappa;
  matrix[t_tb, r] reg_nu;
  matrix[t_tb, r] reg_xi;

  real<lower=0, upper = bp_init_kappa/2> bp_kappa = bp_init_kappa/2;
  real<lower=0, upper = bp_init_nu/2> bp_nu = bp_init_nu/2;
  real<lower=0, upper = bp_init_xi/2> bp_xi = bp_init_xi/2;
  real<lower=0, upper = tau_init_kappa/2> tau_kappa = tau_init_kappa/2;
  real<lower=0, upper = tau_init_nu/2> tau_nu = tau_init_nu/2;
  real<lower=0, upper = tau_init_xi/2> tau_xi = tau_init_xi/2;

  matrix[r, r] corr_kappa = l3 + rho2_kappa * l2 + rho1_kappa * l1;
  matrix[p, p] cov_ar1_kappa = equal + bp_kappa * bp_lin + bp_kappa^2 * bp_square + bp_kappa^3 * bp_cube + bp_kappa^4 * bp_quart;

  matrix[r, r] corr_nu = l3 + rho2_nu * l2 + rho1_nu * l1;
  matrix[p, p] cov_ar1_nu = equal + bp_nu * bp_lin + bp_nu^2 * bp_square + bp_nu^3 * bp_cube + bp_nu^4 * bp_quart;
  
  matrix[r, r] corr_xi = l3 + rho2_xi * l2 + rho1_xi * l1;
  matrix[p, p] cov_ar1_xi = equal + bp_xi * bp_lin + bp_xi^2 * bp_square + bp_xi^3 * bp_cube + bp_xi^4 * bp_quart;

  vector<lower = 0>[N_all] kappa;
  vector<lower = 0>[N_all] nu;
  vector<lower = 0>[N_all] xi;
  vector<lower = 0>[N_all] sigma;

  y[ii_obs] = y_obs;
  y[ii_mis] = y_mis;

  phi_kappa[1,] = (1/tau_kappa) * phi_init_kappa[1]';
  phi_nu[1,] = (1/tau_nu) * phi_init_nu[1]';
  phi_xi[1,] = (1/tau_xi) * phi_init_xi[1]';
  for (j in 2:t_tb) {
    phi_kappa[j,] = eta_kappa * phi_kappa[j-1,] + (1/tau_kappa) * phi_init_kappa[j]';
    phi_nu[j,] = eta_nu * phi_nu[j-1,] + (1/tau_nu) * phi_init_nu[j]';
    phi_xi[j,] = eta_xi * phi_xi[j-1,] + (1/tau_xi) * phi_init_xi[j]';
  }
  
  for (i in 1:r) {
    reg_kappa[, i] = X[i] * beta_kappa[, i]/4 + phi_kappa[, i]/10;
    reg_nu[, i] = X[i] * beta_nu[, i]/4 + phi_nu[, i]/10;
    reg_xi[, i] = X[i] * beta_xi[, i]/8 + phi_xi[, i]/10;
  }

  kappa = exp(to_vector(reg_kappa))[ii_all];
  nu = exp(to_vector(reg_nu))[ii_all];
  xi = exp(to_vector(reg_xi))[ii_all];
  sigma = nu ./ (1 + xi);
}

model {
  // priors
  bp_init_kappa ~ uniform(0, 1);
  bp_init_nu ~ uniform(0, 1);
  bp_init_xi ~ uniform(0, 1);
  
  rho1_kappa ~ beta(3, 4);
  rho2_kappa ~ beta(1.5, 4);
  rho1_nu ~ beta(3, 4);
  rho2_nu ~ beta(1.5, 4);
  rho1_xi ~ beta(3, 4);
  rho2_xi ~ beta(1.5, 4);
  
  target += matnormal_lpdf(beta_kappa | cov_ar1_kappa, corr_kappa);
  target += matnormal_lpdf(beta_nu | cov_ar1_nu, corr_nu);
  target += matnormal_lpdf(beta_xi | cov_ar1_xi, corr_xi);
  
  // IAR prior
  eta_kappa ~ beta(2,8);
  eta_nu ~ beta(2,8);
  eta_xi ~ beta(2,8);
  tau_init_kappa ~ exponential(1);
  tau_init_nu ~ exponential(1);
  tau_init_xi ~ exponential(1);
  for (j in 1:t_tb) {
    target += -.5 * dot_self(phi_init_kappa[j][node1] - phi_init_kappa[j][node2]);
    sum(phi_init_kappa[j]) ~ normal(0, 0.001*r);
    target += -.5 * dot_self(phi_init_nu[j][node1] - phi_init_nu[j][node2]);
    sum(phi_init_nu[j]) ~ normal(0, 0.001*r);
    target += -.5 * dot_self(phi_init_xi[j][node1] - phi_init_xi[j][node2]);
    sum(phi_init_xi[j]) ~ normal(0, 0.001*r);
  }
  // 
  // likelihood
  for (k in 1:N_all) {
    target += egpd_g1_lpdf(y[k] | sigma[k], xi[k], kappa[k]);
  }
}
