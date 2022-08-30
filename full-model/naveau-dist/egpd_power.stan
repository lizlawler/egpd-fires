functions {
  real egpd_power_lpdf(real y, real sigma, real xi, real kappa) {
    if (!(sigma >= 0)) {
      reject("egpd_power_lpdf(y, sigma, xi, kappa): sigma must be positive; found sigma = ", sigma);
    }
    return log(kappa) - log(sigma) - (1/xi + 1) * log(1 + xi * (y/sigma)) + (kappa-1) * log(1 - (1 + xi * (y/sigma))^(-1/xi));
  }
}
data {
  int<lower = 1> N; // # spatial units
  int<lower = 1> T; // # timesteps
  int<lower = 1> p; // # columns in design matrix
 
  int n_fire;
  vector[n_fire] sizes;
  int<lower = 1, upper = N * T> burn_idx[n_fire];

  // full sparse matrix for all ecoregions X timesteps
  int<lower = 1> n_w;
  vector[n_w] w;
  int<lower = 1> v[n_w];
  int<lower = 1> u[N * T + 1];

  // design matrix for burn sizes
  int<lower = 1> n_w_tb;
  vector[n_w_tb] w_tb;
  int<lower = 1> v_tb[n_w_tb];
  int<lower = 1, upper = N * T + 1> n_u_tb;
  int<lower = 1> u_tb[n_u_tb];

  int<lower = 3, upper = 3> M; // num dimensions of betas

  int<lower = 0> n_edges;
  int<lower = 1, upper = N> node1[n_edges];
  int<lower = 1, upper = N> node2[n_edges];
  int<lower = 1, upper = N*T> b_idx[n_u_tb - 1]; // burn spatial idx
}

transformed data {
  int n_burn_mu = n_u_tb - 1;
}

parameters {
  matrix[M, p] betaR;
  vector[M] alpha; // intercept

  // double generalized pareto params
  vector<lower = 0>[M] sigmaR;     // DGP dist scale, un-scaled
  vector<lower = 0>[M] lambda_beta;    // exponential prior 

  // iar params
  // vector<lower = 0, upper = 1>[M] eta;   // autoregressive param
  // vector<lower = 0>[M] sigma_phi;       // time difference spatial sd
  // vector[N] phiR[M, T];                  // unscaled values
}


transformed parameters {
  matrix[M, p] beta;
  vector<lower = 0>[M] sigma_beta;
  vector[n_fire] kappa;
  vector[n_fire] nu;
  vector[n_fire] xi;
  vector[n_fire] sigma_egpd;

  // parameters for IAR on phi
  // matrix[T, N] phi[M];
  // vector[N * T] phi_vec[M];
  // 
  // for (i in 1:M) {
  //   phi[i][1] = phiR[i, 1]' * sigma_phi[i];
  //   for (t in 2:T) {
  //     // subsequent timesteps
  //     phi[i][t] = eta[i] * phi[i][t - 1] + phiR[i, t]' * sigma_phi[i];
  //   }
  //   phi_vec[i] = to_vector(phi[i]');
  // }

  // double generalized pareto prior on betas
  for (i in 1:M) {
    sigma_beta[i] = lambda_beta[i]^2/2 * sqrt(sigmaR[i]);
    beta[i, ] = betaR[i, ] * sigma_beta[i];
  }
    
  kappa = exp(alpha[1]
             + csr_matrix_times_vector(n_burn_mu, p, w_tb, v_tb, u_tb, beta[1, ]'));
             // + phi_vec[1][b_idx]);
  nu = exp(alpha[2] 
             + csr_matrix_times_vector(n_burn_mu, p, w_tb, v_tb, u_tb, beta[2, ]')); 
             // + phi_vec[2][b_idx]);
  xi = exp(alpha[3] 
             + csr_matrix_times_vector(n_burn_mu, p, w_tb, v_tb, u_tb, beta[3, ]')); 
             // + phi_vec[3][b_idx]);
  for (i in 1:n_fire) {
    sigma_egpd[i] = nu[i]/(1 + xi[i]);
  }
}

model {
  to_vector(betaR) ~ normal(0, 1);
  sigmaR ~ exponential(1);
  lambda_beta ~ gamma(1,1);
  
  alpha ~ normal(0, 5);
  // sigma_phi ~ normal(0, 1);
  // eta ~ beta(8, 2);
  // for (i in 1:M) {
  //   for (t in 1:T) {
  //     // IAR prior
  //     target += -0.5 * dot_self(phiR[i, t][node1] - phiR[i, t][node2]);
  //     sum(phiR[i, t]) ~ normal(0, .001 * N);
  //   }
  // }

  // fire sizes
  for (i in 1:n_fire) {
    target += egpd_power_lpdf(sizes[i] | sigma_egpd[i], xi[i], kappa[i]);
  }
}
