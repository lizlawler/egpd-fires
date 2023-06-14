functions{
#include /../../twcrps_matnorm_fcns.stan
#include /../../gpd_fcns.stan
  // custom distribution functions of EGPD
  real egpd_lpdf(real y, real sigma, real xi, real kappa) {
    if (kappa > 1e-15) {
      return log(kappa) + (kappa-1) * gpareto_lcdf(y | sigma, xi) + 
              gpareto_lpdf(y | sigma, xi);
    }
    else {
      reject("kappa<=0; found kappa = ", kappa);
    }
  }  
  real egpd_cdf(real y, real sigma, real xi, real kappa) {
    if (kappa > 1e-15) {
      return exp(kappa * gpareto_lcdf(y | sigma, xi));
    }
    else {
      reject("kappa<=0; found kappa = ", kappa);
    }
  }
  real egpd_lcdf(real y, real sigma, real xi, real kappa) {
    if (kappa > 1e-15) {
      return kappa * gpareto_lcdf(y | sigma, xi);
    }
    else {
      reject("kappa<=0; found kappa = ", kappa);
    }
  }  
  real egpd_lccdf(real y, real sigma, real xi, real kappa) {
    if (kappa > 1e-15) {
      return log1m_exp(kappa * gpareto_lcdf(y | sigma, xi));
    }
    else {
      reject("kappa<=0; found kappa = ", kappa);
    }
  }  
  real egpd_icdf(real u_adj, real sigma, real xi, real kappa) {
    // input is 'u' adjusted for lower bound truncation
    if (kappa > 1e-15) {
      return gpareto_icdf(u_adj^(1/kappa), sigma, xi);
    }
    else {
      reject("kappa<=0; found kappa = ", kappa);
    }
  }  
  // truncated EGPD distribution
  real egpd_trunc_lpdf(real y, real ymin, real sigma, real xi, real kappa) {
    real lpdf = egpd_lpdf(y | sigma, xi, kappa);
    real cst = egpd_lccdf(ymin | sigma, xi, kappa);
    return lpdf - cst;
  }
  // forecast function
  // vector forecast_rng(int n_pred, real ymin, real sigma, real xi, real kappa) {
  //   vector[n_pred] forecast;
  //   vector[n_pred] a = rep_vector(0, n_pred);
  //   vector[n_pred] b = rep_vector(1, n_pred);
  //   array[n_pred] real u = uniform_rng(a, b);
  //   real cst = exp(egpd_lcdf(ymin | sigma, xi, kappa));
  //   for (n in 1:n_pred) {
  //     real u_adj = u[n] * (1 - cst) + cst;
  //     forecast[n] = egpd_icdf(u_adj, sigma, xi, kappa);
  //   }
  //   return forecast;
  // }
}