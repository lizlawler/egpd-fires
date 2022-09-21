functions {
  real egpd_g1_lpdf(real y, real sigma, real xi, real kappa) {
    real lpdf;
    lpdf = log(kappa) - log(sigma) - (1/xi + 1) * log(1 + xi * (y/sigma)) + 
    (kappa-1) * log(1 - (1 + xi * (y/sigma))^(-1/xi));
    return lpdf;
  }
}

data {
  int<lower = 1> p; // # of parameters
  int<lower = 1> T; // total # of timepoints
  int<lower = 1> R; // total # of regions
  int<lower=0> n_edges;
  int<lower=1, upper = R> node1[n_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper = R> node2[n_edges];  // and node1[i] < node2[i]

  // training data
  int<lower = 1> burn_er; // # of regions with nonzero burn areas
  int burn_er_idx[burn_er];
  int<lower = 1> burn_train; // # of total observations in training set
  int<lower = 1> cov_rows; // # of unique rows for X_train
  real<lower = 0> burn_size[burn_train]; // vector of burned areas
  matrix[cov_rows, p] X;
  int segment_length[burn_er];
  int<lower = 1, upper = cov_rows> rows_X[cov_rows];
  int<lower = 1> time_idx_er[cov_rows];
  int<lower = 1> burns_to_covar_idx[burn_train]; // match burn areas to covariates
  
  // indicator matrices for ecoregions - training
  matrix[burn_er, burn_er] l3;
  matrix[burn_er, burn_er] l2;
  matrix[burn_er, burn_er] l1;

  // indicator matrices for AR(1) process on betas
  matrix[p, p] equal;
  matrix[p, p] bp_lin;
  matrix[p, p] bp_square;
  matrix[p, p] bp_cube;
  matrix[p, p] bp_quart;
  
  int<lower=1, upper = 3> M; // dimensions for egpd params; if all 3 will include regression, then M = 3
}

// may want to add "transformed data" block to standardize X and y; remember that "./" us elementwise division

parameters {
  vector[R] phi_init[T, M];
  matrix[p, burn_er] Z[M];
  real<lower = 0> tau_init[M];
  real<lower = 0, upper = 1> eta[M];
  real<lower = 0, upper = 1> bp_init[M];
  real<lower = 0, upper = 1> rho2[M];
  real<lower=0, upper = (1-rho2[M])> rho1[M];
}

transformed parameters {
  matrix[T, R] phi[M];

  real<lower=0, upper = 0.5> bp[M];
  real<lower=0> tau[M];
  
  corr_matrix[burn_er] corr[M];
  cov_matrix[p] cov_ar1[M];
  matrix[p, burn_er] beta[M];

  for (i in 1:M) {
    bp[i] = bp_init[i]/2;
    tau[i] = tau_init[i]/2;
    corr[i] = l3 + rho2[i] * l2 + rho1[i] * l1;
    cov_ar1[i] = equal + bp[i] * bp_lin + bp[i]^2 * bp_square + bp[i]^3 * bp_cube + bp[i]^4 * bp_quart;
    beta[i] = cholesky_decompose(cov_ar1[i])' * Z[i] * cholesky_decompose(corr[i]);
    phi[i][1,] = (1/tau[i]) * phi_init[1, i]';
    for (j in 2:T) {
      phi[i][j, ] = eta[i] * phi[i][j-1,] + (1/tau[i]) * phi_init[j, i]';
    }
  }
}

model {
  vector[cov_rows] regress[M];
  vector[burn_train] kappa;
  vector[burn_train] nu;
  vector[burn_train] xi;
  vector[burn_train] sigma;
  
  int pos;
  for (i in 1:M) {
    pos = 1;
    for(k in 1:burn_er) {
      regress[i][segment(rows_X, pos, segment_length[k])] = X[segment(rows_X, pos, segment_length[k]), ] * beta[i][,k] + phi[i][segment(time_idx_er, pos, segment_length[k]), burn_er_idx[k]];
      pos += segment_length[k];
    }
  }
  
  kappa = exp(regress[1])[burns_to_covar_idx];
  nu = exp(regress[2])[burns_to_covar_idx];
  xi = exp(regress[3])[burns_to_covar_idx];
  sigma = nu ./ (1 + xi);
  
  // priors
  for (i in 1:M) {
    to_vector(Z[i]) ~ normal(0, 1);
    bp_init[i] ~ uniform(0, 1);
    rho1[i] ~ beta(1.5, 4);
    rho2[i] ~ beta(3, 4);
    
    // IAR prior
    eta[i] ~ beta(2,8);
    tau_init[i] ~ exponential(1);
    for (j in 1:T) {
      target += -.5 * dot_self(phi_init[j, i][node1] - phi_init[j, i][node2]);
      sum(phi_init[j, i]) ~ normal(0, 0.001*R);
    }
  }
  // likelihood
  for (k in 1:burn_train) {
    target += egpd_g1_lpdf(burn_size[k] | sigma[k], xi[k], kappa[k]);
  }
}

// generated quantities {
//   vector[burn_train] train_loglik;
//   
//     // training log likelihood
//   for (k in 1:burn_train) {
//     train_loglik[k] = egpd_g1_lpdf(burn_size[k] | sigma[k], xi[k], kappa[k]);
//   }
// }
