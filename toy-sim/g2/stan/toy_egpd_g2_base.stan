functions {
  real egpd_g2_lpdf(real y, real sigma, real xi, real kappa1, real kappa2, real prob) {
    real lpdf;
    lpdf =  -(1/xi + 1) * log(1 + xi * (y/sigma)) - log(sigma) + 
    log(kappa1 * prob * (1 - (1 + xi * (y/sigma))^(-1/xi))^(kappa1 - 1) + 
    kappa2 * (1-prob) * (1 - (1 + xi * (y/sigma))^(-1/xi))^(kappa2 - 1));
    return lpdf;
  }
}

data {
  int<lower = 1> p; // # of parameters
  int<lower = 1> n; // # of observations
  matrix[n, p] X; // design matrix
  vector[n] y; // response data
}

parameters {
  // need to add in constraint that beta_kappa1 < beta_kappa2
  vector[p] beta_kappa1;
  vector[p] beta_kappa2;
  vector[p] beta_nu;
  vector[p] beta_xi;
}

transformed parameters {
  vector[n] kappa1;
  vector[n] kappa2;
  vector[n] nu;
  vector[n] xi;
  vector[n] sigma;
  
  kappa1 = to_vector(exp(X * beta_kappa1));
  kappa2 = to_vector(exp(X * beta_kappa2));
  nu = to_vector(exp(X * beta_nu));
  xi = to_vector(exp(X * beta_xi));
  for (i in 1:n) {
    sigma[i] = nu[i] / (1 + xi[i]);
  }
}

model {
  // priors
  beta_kappa1 ~ normal(0, 1);
  beta_kappa2 ~ normal(0, 1);
  beta_nu ~ normal(0, 1);
  beta_xi ~ normal(0, 1);

  // likelihood
  for (i in 1:n) {
    target += egpd_g2_lpdf(y[i] | sigma[i], xi[i], kappa1[i], kappa2[i], 0.5);
  }
}
