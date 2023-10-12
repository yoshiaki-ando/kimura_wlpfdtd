/*
 * get_coordinate.cpp
 *
 *  Created on: 2020/09/27
 *      Author: ando
 */
#include "wlp_fdtd2d.h"

void set_coordinate(
    double *dx, const int nx, double *x, const double x_orig,
    double *dy, const int ny, double *y, const double y_orig){

  x[0] = x_orig - LX2; 
  y[0] = y_orig - LY2;

  for(int i = 1; i < nx; i++){
    x[i] = x[i-1] + dx[i-1];
  }
  for(int j = 1; j < ny; j++){
    y[j] = y[j-1] + dy[j-1];
  }
}



