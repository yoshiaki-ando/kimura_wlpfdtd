#include "wlp_fdtd2d.h"

void update_e(double **coef_ez,
              double **sum_coef_ez,
              double **sum_coef_ezx, double **sum_coef_ezy,
              double *x,
              double ***node,
              int p
              )
{
  
  for( int i = 1; i < NX7; i++ ){
    for( int j = 1; j < NY7; j++ ){
      coef_ez[i][j] = x[(int)node[i][j][EZ]];
      sum_coef_ez[i][j] += coef_ez[i][j];
      if( node[i][j][EZX] >= 0 ){
        sum_coef_ezx[i][j] += x[(int)node[i][j][EZX]];
        sum_coef_ezy[i][j] += x[(int)node[i][j][EZY]];
      }
    }
  }
  
  return;
}


void update_h(double *dx, double *dy,
              double **coef_hx, double **coef_hy,
              double **sum_coef_hx, double **sum_coef_hy,
              double **coef_ez,
              double *Cx, double *Cy,
              double *const_hx1, double *const_hx2,
              double *const_hy1, double *const_hy2,
              double ***node,
              int p
              ){
  
  const int l = L_PML;
  
  /*** hx(i,j+0.5) ***********************************************/
  for( int i = 0; i <= NX7; i++ ){
    for( int j = 0; j < NY7; j++ ){
      if( (i <= l) || (NX7-l <= i) ||
         (j < l) || (NY7-l <= j) ){
         // Eq.(82)
        coef_hx[i][j]
        = ( - const_hx1[j] * ( coef_ez[i][j+1]-coef_ez[i][j] )
           - const_hx2[j] * sum_coef_hx[i][j] );
      } else {
        // Eq.(16)
        coef_hx[i][j]
        = ( - Cy[j] * ( coef_ez[i][j+1]-coef_ez[i][j] )
           - 2.0*sum_coef_hx[i][j] );
      }
      sum_coef_hx[i][j] += coef_hx[i][j];
    }
  }
  /*** hx(i,j+0.5) (up to here)  *********************************/
  
  
  
  /*** hy(i+0.5,j) ***********************************************/
  for( int i = 0; i < NX7; i++ ){
    for( int j = 0; j <= NY7; j++ ){
      if( (i < l) || (NX7-l <= i) ||
         (j <= l) || (NY7-l <= j) ){
        // Eq.(83)
        coef_hy[i][j]
        = const_hy1[i] * ( coef_ez[i+1][j]-coef_ez[i][j] )
        - const_hy2[i] * sum_coef_hy[i][j];
      } else {
        // Eq.(17)
        coef_hy[i][j]
        = Cx[i] * (coef_ez[i+1][j] - coef_ez[i][j])
        - 2.0*sum_coef_hy[i][j];
      }
      sum_coef_hy[i][j] += coef_hy[i][j];
    }
  }
  /*** Hy(i+0.5,j) (up to here)  *********************************/
  
  return;
}

