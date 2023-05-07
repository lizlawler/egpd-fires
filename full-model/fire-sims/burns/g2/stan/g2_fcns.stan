functions{
#include /../../twcrps_matnorm_fcns.stan
#include /../../gpd_fcns.stan
  // custom distribution functions of EGPD
  real egpd_lpdf(real y, real sigma, real xi, real kappa1, real kappa2, real prob) {
    if (prob < 0 || prob > 1) {
      reject("not a valid probability; found prob = ", prob);
    }
    // else if (kappa1 > kappa2) {
    //   reject("kappa1 > kappa2; found kappa1 = ", kappa1, ", kappa2 = ", kappa2);
    // }
    else if (kappa1 > 1e-15 && kappa2 > 1e-15) {
      return gpareto_lpdf(y | sigma, xi) + log_sum_exp(
              log(kappa1) + (kappa1-1) * (log(prob) + gpareto_lcdf(y | sigma, xi)),
              log(kappa2) + (kappa2-1) * (log1m(prob) + gpareto_lcdf(y | sigma, xi)));
    }
    else {
      reject("kappa1 or kappa2 <=0; found kappa1 = ", kappa1, ", kappa2 = ", kappa2);
    }
  }  
  real egpd_cdf(real y, real sigma, real xi, real kappa1, real kappa2, real prob) {
    if (prob < 0 || prob > 1) {
      reject("not a valid probability; found prob = ", prob);
    }
    // else if (kappa1 > kappa2) {
    //   reject("kappa1 > kappa2; found kappa1 = ", kappa1, ", kappa2 = ", kappa2);
    // }
    else if (kappa1 > 1e-15 && kappa2 > 1e-15) {
      return exp(log(prob) + kappa1 * gpareto_lcdf(y | sigma, xi)) + 
              exp(log1m(prob) + kappa2 * gpareto_lcdf(y | sigma, xi));
    }
    else {
      reject("kappa1 or kappa2 <=0; found kappa1 = ", kappa1, ", kappa2 = ", kappa2);
    }
  }
  real egpd_lcdf(real y, real sigma, real xi, real kappa1, real kappa2, real prob) {
    if (prob < 0 || prob > 1) {
      reject("not a valid probability; found prob = ", prob);
    }
    // else if (kappa1 > kappa2) {
    //   reject("kappa1 > kappa2; found kappa1 = ", kappa1, ", kappa2 = ", kappa2);
    // }
    else if (kappa1 > 1e-15 && kappa2 > 1e-15) {
      return log_sum_exp(
              log(prob) + kappa1 * gpareto_lcdf(y | sigma, xi), 
              log1m(prob) + kappa2 * gpareto_lcdf(y | sigma, xi));
    }
    else {
      reject("kappa1 or kappa2 <=0; found kappa1 = ", kappa1, ", kappa2 = ", kappa2);
    }
  }  
  real egpd_lccdf(real y, real sigma, real xi, real kappa1, real kappa2, real prob) {
    if (prob < 0 || prob > 1) {
      reject("not a valid probability; found prob = ", prob);
    }
    // else if (kappa1 > kappa2) {
    //   reject("kappa1 > kappa2; found kappa1 = ", kappa1, ", kappa2 = ", kappa2);
    // }
    else if (kappa1 > 1e-15 && kappa2 > 1e-15) {
      return log1m(exp(log(prob) + kappa1 * gpareto_lcdf(y | sigma, xi)) + 
              exp(log1m(prob) + kappa2 * gpareto_lcdf(y | sigma, xi)));
    }
    else {
      reject("kappa1 or kappa2 <=0; found kappa1 = ", kappa1, ", kappa2 = ", kappa2);
    }
  }  
  // truncated EGPD distribution
  real egpd_trunc_lpdf(real y, real ymin, real sigma, real xi, real kappa1, real kappa2, real prob) {
    real lpdf = egpd_lpdf(y | sigma, xi, kappa1, kappa2, prob);
    real cst = egpd_lccdf(ymin | sigma, xi, kappa1, kappa2, prob);
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