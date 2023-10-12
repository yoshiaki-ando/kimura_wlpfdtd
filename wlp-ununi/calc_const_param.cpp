#include "wlp_fdtd2d.h"

void calc_const_param(double *Dx, double *Dy, double *Cx, double *Cy,
                      double *dx, double *dy, double *dx_h, double *dy_h){

  const double s = TIME_SCALE;

  for(int i=1; i < NX7; i++){
  	Dx[i] = 2.0/EPS0/s/dx_h[i];
  }
  for(int j=1; j < NY7; j++){
  	Dy[j] = 2.0/EPS0/s/dy_h[j];
  }
  
  for(int i=0; i < NX7; i++){
  	Cx[i] = 2.0/MU0/s/dx[i];
  }
  for(int j=0; j < NY7; j++){
  	Cy[j] = 2.0/MU0/s/dy[j];
  }

  return;
}