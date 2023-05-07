functions{
#include /../../twcrps_matnorm_fcns.stan
  // truncated lognormal distribution
  real lognorm_trunc_lpdf(real y, real ymin, real mu, real sigma) {
    real cst = lognormal_lccdf(ymin | mu, sigma);
    real lpdf = lognormal_lpdf(y | mu, sigma);
    return lpdf - cst;
  }
  //forecast function
  // vector forecast_rng(int n_pred, real ymin, real sigma, real xi, real kappa) {
  //   vector[n_pred] forecast;
  //   vector[n_pred] a = rep_vector(0, n_pred);
  //   vector[n_pred] b = rep_vector(1, n_pred);
  //   array[n_pred] real u = uniform_rng(a, b);
  //   real cst = (1 - (1 + xi * (ymin / sigma)) ^ (-1 / xi)) ^ kappa;
  //   for (n in 1:n_pred) {
  //     real u_adj = u[n] * (1 - cst) + cst;
  //     forecast[n] = (sigma / xi) * ((1 - u_adj ^ (1 / kappa)) ^ -xi - 1);
  //   }
  //   return forecast;
  // }
}