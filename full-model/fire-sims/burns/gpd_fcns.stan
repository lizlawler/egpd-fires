  // custom distribution functions for GPD
  real gpareto_lpdf(real y, real sigma, real xi) {
    real inv_xi = inv(xi);
    if (xi < 0 && y/sigma > -inv_xi) {
      reject("xi<0 and y/sigma > -1/xi; found xi = ", xi, " sigma = ", sigma);
    }
    else if (sigma < 1e-15) {
      reject("sigma<=0; found sigma = ", sigma);
    }
    else if (abs(xi) > 1e-15) {
      return -(1+inv_xi) * log1p(y * (xi/sigma)) - log(sigma);
    }
    else {
      return -log(sigma) - (y/sigma); // limit xi->0
    }
  }
  real gpareto_cdf(real y, real sigma, real xi) {
    real inv_xi = inv(xi);
    if (xi < 0 && y/sigma > -inv_xi) {
      reject("xi<0 and y/sigma > -1/xi; found xi = ", xi, " sigma = ", sigma);
    }
    else if (sigma < 1e-15) {
      reject("sigma<=0; found sigma = ", sigma);
    }
    else if (abs(xi) > 1e-15) {
      return exp(log1m_exp(-inv_xi * log1p(y * (xi/sigma))));
    }
    else {
      return exp(log1m_exp(-y/sigma)); // limit xi->0
    }
  }
  real gpareto_lcdf(real y, real sigma, real xi) {
    real inv_xi = inv(xi);
    if (xi < 0 && y/sigma > -inv_xi) {
      reject("xi<0 and y/sigma > -1/xi; found xi = ", xi, " sigma = ", sigma);
    }
    else if (sigma < 1e-15) {
      reject("sigma<=0; found sigma = ", sigma);
    }
    else if (abs(xi) > 1e-15) {
      return log1m_exp(-inv_xi * log1p(y * (xi/sigma)));
    }
    else {
      return log1m_exp(-y/sigma); // limit xi->0
    }
  }
  real gpareto_lccdf(real y, real sigma, real xi) {
    real inv_xi = inv(xi);
    if (xi < 0 && y/sigma > -inv_xi) {
      reject("xi<0 and y/sigma > -1/xi; found xi = ", xi, " sigma = ", sigma);
    }
    else if (sigma < 1e-15) {
      reject("sigma<=0; found sigma = ", sigma);
    }
    else if (abs(xi) > 1e-15) {
      return -inv_xi * log1p(y * (xi/sigma));
    }
    else {
      return -y/sigma; // limit xi->0
    }
  }
  real gpareto_icdf(real u, real sigma, real xi) {
    // input is adjusted 'u' with egpd carrier inverse applied
    if (sigma< 1e-15)
      reject("sigma<=0; found sigma =", sigma);
    else if (abs(xi) > 1e-15)
      return (sigma / xi) * ((1-u)^-xi - 1);
    else
      return -sigma * log1m(u); // limit xi->0
  }