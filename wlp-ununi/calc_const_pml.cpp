#include "wlp_fdtd2d.h"


void calc_const_pml(double *dx, double *dy,
                    double *const_hx1, double *const_hx2,
		                double *const_hy1, double *const_hy2,
                    double *conduct_pmlx_ez, double *conduct_pmly_ez,
                    double *conduct_pmlx_h, double *conduct_pmly_h
                    ){

  const double s = TIME_SCALE;
  double sig = 0.0; // [/s] (= [S/m]/[F/m])

  for( int j = 0; j < NY7; j++ ){
    sig = conduct_pmly_h[j];

    // Eq.(55)
    const_hx1[j] = 1.0/MU0/dy[j] / (s/2.0 + sig);

    // Eq.(56)
    const_hx2[j] = s / (s/2.0 + sig);
  }

  for( int i = 0; i < NX7; i++ ){
    sig = conduct_pmlx_h[i];

    // Eq.(57)
    const_hy1[i] = 1.0/MU0/dx[i] / (s/2.0 + sig);

    // Eq.(58)
    const_hy2[i] = s / (s/2.0 + sig);
  }

  return;
}
