functions {
  // real matnormal_lpdf(data matrix y, matrix u, matrix v) {
  //   real lpdf;
  //   real n;
  //   real p;
  //   n = rows(u);
  //   p = rows(v);
  //   lpdf = -(n*p/2) * log(2 * pi()) - (n/2)*log_determinant(v) - (p/2)*log_determinant(u) -
  //         0.5 * trace(mdivide_right_spd(mdivide_left_spd(v, y'), u) * y);
  //   return lpdf;
  // }

}

data {
  int<lower = 1> p; // # of parameters
  int<lower = 1> n; // # of observations
  int<lower = 1> r; // # of regions
  
  matrix[r, p] betas; // data
  
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

// should consider transforming X and y so they're centered and scaled

parameters {
  real<lower = 0, upper = 1> bp_init;
  real<lower = 0, upper = 1> rho2;
  real<lower=0, upper = (1-rho2)> rho1;
  matrix[r, p] Z;
}

transformed parameters {
  real<lower=0, upper = bp_init/2> bp = bp_init/2;
  
  matrix[r, r] corr = l3 + rho2 * l2 + rho1 * l1;
  matrix[p, p] cov_ar1 = equal + bp * bp_lin + bp^2 * bp_square + bp^3 * bp_cube + bp^4 * bp_quart;
}

model {
  matrix[r, p] Z_rev;
  // priors
  bp_init ~ uniform(0, 1);
  rho1 ~ beta(1.5, 4);
  rho2 ~ beta(3, 4);
  
  // likelihood
  Z_rev = inverse(cholesky_decompose(corr)') * betas * inverse(cholesky_decompose(cov_ar1));
  to_vector(Z_rev) ~ normal(0, 1);
}
