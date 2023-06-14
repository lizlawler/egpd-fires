  // twCRPS and matnormal_lpdf remain unchanged across models
  // real twCRPS(real y, vector forecast, real delta, real w_mean, real w_sd) {
  //   real score;
  //   real summand;
  //   int N = rows(forecast);
  //   summand = 0;
  //   for (n in 1:N) {
  //     summand += (forecast[n] - step(forecast[n] - y)) ^ 2
  //                * normal_cdf(forecast[n] | w_mean, w_sd);
  //   }
  //   score = summand * delta;
  //   return score;
  // }
  
  // probability forecast
  vector prob_forecast(int n_int, vector int_pts, real ymin, real sigma, real xi, real kappa) {
    vector[n_int] pred_probs;
    real numer_cst = egpd_cdf(ymin | sigma, xi, kappa);
    real denom_cst = exp(egpd_lccdf(ymin | sigma, xi, kappa));
    for (n in 1:n_int) {
      real unnorm_pred_prob = egpd_cdf(int_pts[n] | sigma, xi, kappa);
      pred_probs[n] = (unnorm_pred_prob - numer_cst)/denom_cst;
    }
    return pred_probs;
  }  
  
  real twCRPS(real y, int n_int, real y_int, vector int_pts, real ymin, real sigma, real xi, real kappa, real w_mean, real w_sd) {
    real score;
    real summand;
    real delta = y_int/(n_int - 1);
    vector[n_int] pred_probs = prob_forecast(n_int, int_pts, ymin, sigma, xi, kappa);
    summand = 0;
    for (n in 1:n_int) {
      summand += (pred_probs[n] - step(int_pts[n] - y)) ^ 2
                 * normal_cdf(int_pts[n] | w_mean, w_sd);
    }
    score = summand * delta;
    return score;
  }
  
  real matnormal_lpdf(matrix y, matrix cov, matrix corr) {
    real lpdf;
    real r = rows(corr);
    real p = rows(cov);
    lpdf = -(r * p / 2) * log(2 * pi()) - (p / 2) * log_determinant(corr)
           - (r / 2) * log_determinant(cov)
           - 0.5 * trace(mdivide_right_spd(mdivide_left_spd(corr, y'), cov) * y);
    return lpdf;
  }