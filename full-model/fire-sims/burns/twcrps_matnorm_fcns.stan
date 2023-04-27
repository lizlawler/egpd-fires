  // twCRPS and matnormal_lpdf remain unchanged across models
  real twCRPS(real y, vector forecast, real delta, real w_mean, real w_sd) {
    real score;
    real summand;
    int N = rows(forecast);
    summand = 0;
    for (n in 1:N) {
      summand += (forecast[n] - step(forecast[n] - y)) ^ 2
                 * normal_cdf(forecast[n] | w_mean, w_sd);
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