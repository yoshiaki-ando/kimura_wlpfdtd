#include "wlp_fdtd2d.h"

#include <algorithm>

double weighting(const double x, const double w, const double t){
  return x * std::exp( -w - 0.5*t );
}

double wlp(const int p, const double t){
  /* Weighted Laguerre Polynomial 𝜑ₙ(z) */
  double ans;

  double *v = new double [p+1];
  double w[3] = { 0.0, 0.0, 0.0 };

  v[0] = 1.0;

  if ( p > 0 ){
    v[1] = 1.0 - t;

    for(int n = 2; n <= p; n++){
      w[n%3] = w[(n-1)%3];

      /* exp(10)を超える毎に exp(-10)倍する。何倍したか w に記憶しておく */
      if ( std::log(std::abs(v[n-1])) > 10.0 ){
        v[n-1] *= std::exp(-10.0);
        v[n-2] *= std::exp(-10.0);
        w[0] -= 10.0;
        w[1] -= 10.0;
        w[2] -= 10.0;
      }

      v[n] = (2.0*n - 1.0 - t)/n * v[n-1] - (n - 1.0)/n * v[n-2];
    }
  }

  /* exp(-10 * N)倍を戻すのと重み exp(-t/2) を同時に戻す */
  ans = weighting(v[p], w[p%3], t);
  delete [] v;

  return ans;
}

void calc_wlp_for_output(const int p, double **w, double **v, double *ret, const double dt, const int max_n){
  /* Weighted Laguerre Polynomial 𝜑ₙ(z) */
  // double ret* = new double [NT+1];
  const double s = TIME_SCALE;

  if(p == 0){
    for(int n = 0; n <= max_n; n++) v[0][n] = 1.0;
  }else if(p == 1) {
    for(int n = 0; n <= max_n; n++) v[1][n] = 1.0-s*n*dt;
  } else if(p >= 2) {
    //w[n] = w[n-1];
    for(int n = 0; n <= max_n; n++){
      w[p%3][n] = w[(p-1)%3][n];
      /* exp(10)を超える毎に exp(-10)倍する。何倍したか w に記憶しておく */
      if ( std::log(std::abs(v[(p-1)%3][n])) > 10.0 ){
        v[(p-1)%3][n] *= std::exp(-10.0);
        v[(p-2)%3][n] *= std::exp(-10.0);
        w[0][n] -= 10.0;
        w[1][n] -= 10.0;
        w[2][n] -= 10.0;
      }
      v[p%3][n] = (2.0*p - 1.0 - s*n*dt)/p * v[(p-1)%3][n] - (p - 1.0)/p * v[(p-2)%3][n];
    }
  }

  /* exp(-10 * N)倍を戻すのと重み exp(-t/2) を同時に戻す */
  for(int n = 0; n <= max_n; n++){
    ret[n] = weighting(v[p%3][n], w[p%3][n], s*n*dt);
  }

  return;
}