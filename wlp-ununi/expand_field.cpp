#include "wlp_fdtd2d.h"

//p次の展開係数を時間波形に展開
void expand_field_p(double coef_field, double *field, 
                    double *wlp,
                    int p){

  for(int n = 0; n < NT+1; n++){ 
    field[n] += coef_field * wlp[n];
  }

  return;
}

//p次の展開係数を時間波形に展開(2次元)
void expand_field_2d_p(double **coef_field, double ***field,
                           int nx, int ny, int nt,
                           double *wlp,
                           int p){

  for(int i = 0; i < nx; i++){
    for(int j = 0; j < ny; j++){
      for(int n = 0; n < nt+1; n++){ 
        field[i][j][n] += coef_field[i+NX3][j+NY3] * wlp[n];
      }
    }
  }

  return;
}