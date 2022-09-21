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
  int<lower = 1> T; // total # of timepoints
  int<lower = 1> R; // total # of regions
  int<lower=0> n_edges;
  int<lower=1, upper = R> node1[n_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper = R> node2[n_edges];  // and node1[i] < node2[i]

  matrix[T, p] X[R]; // design matrix; 1-D array of size r with matrices t x p; this indexing is deprecated in version 2.32 and above
  vector[T*R] y; // response data

  // indicator matrices for ecoregions - training
  matrix[R, R] l3;
  matrix[R, R] l2;
  matrix[R, R] l1;

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
  real<lower = 0> tau_init[M];
  real<lower = 0, upper = 1> eta[M];
  real<lower = 0, upper = 1> bp_init[M];
  real<lower = 0, upper = 1> rho2[M];
  real<lower=0, upper = (1-rho2[M])> rho1[M];
}

transformed parameters {
  matrix[T, R] phi[M];
  matrix[R, T] regress[M];
  vector[R * T] kappa;
  vector[R * T] nu;
  vector[R * T] xi;
  vector[R * T] sigma;

  real<lower=0, upper = 0.5> bp[M];
  real<lower=0> tau[M];
  
  corr_matrix[R] corr[M];
  cov_matrix[p] cov_ar1[M];
  matrix[R, p] betas[M];

  for (i in 1:M) {
    bp[i] = bp_init[i]/2;
    tau[i] = tau_init[i]/2;
    corr[i] = l3 + rho2[i] * l2 + rho1[i] * l1;
    cov_ar1[i] = equal + bp[i] * bp_lin + bp[i]^2 * bp_square + bp[i]^3 * bp_cube + bp[i]^4 * bp_quart;
    phi[i][1,] = (1/tau[i]) * phi_init[1, i]';
    for (j in 2:T) {
      phi[i][j, ] = eta[i] * phi[i][j-1,] + (1/tau[i]) * phi_init[j, i]';
    }
    for (k in 1:R) {
      regress[i][k, ] = (X[k] * betas[i][k, ]'/4 + phi[i][, k]/10)';
    }
  }
  kappa = to_vector(exp(regress[1]));
  nu = to_vector(exp(regress[2]));
  xi = to_vector(exp(regress[3]));
  sigma = nu ./ (1 + xi);
}

model {
  for (i in 1:M) {
    target += matnormal_lpdf(betas[i] | corr[i], cov_ar1[i]);
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
  for (k in 1:(T*R)) {
    target += egpd_g1_lpdf(y[k] | sigma[k], xi[k], kappa[k]);
  }
}
