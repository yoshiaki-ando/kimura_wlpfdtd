#include <cstdlib>
#include "wlp_fdtd2d.h"

double *allocate_memory_d1d(const int N_of_Array, const double Initial_value){

  double *mat = new double [N_of_Array];

  for(int i = 0; i < N_of_Array; i++){
    mat[i] = Initial_value;
  }

  return (mat);
}

double **allocate_memory_d2d(const int M, const int N, const double Initial_value){

  double **mat = new double* [M];
  double *pmat = new double [M * N];

  for(int i = 0; i < M; i++){
    mat[i] = pmat + i * N;  // ポインタの設定

    // 初期化
    for(int j = 0; j < N; j++){
      mat[i][j] = Initial_value;
    }

  }

  return (mat);
}

void free_memory_d2d(double **mat){
  delete [] mat[0];
  delete [] mat;
}

double ***allocate_memory_d3d(int ir, int jr, int kr, double Initial_value){

  double *pmat = (double*)malloc( sizeof(double) * ir*jr*kr );
  double **p2mat = (double**)malloc( sizeof(double*) * ir*jr );
  double ***mat = (double***)malloc( sizeof(double**) * ir );

  for(int i = 0; i < ir; i++){

    mat[i] = p2mat + i*jr;
    for(int j = 0; j < jr; j++){
      mat[i][j] = pmat + i*jr*kr + j*kr;
      for(int k = 0; k < kr; k++){
        mat[i][j][k] = Initial_value;
      }
    }
  }

  return(mat);
}

void free_memory_d3d(double ***mat){
  delete [] mat[0][0];
  delete [] mat[0];
  delete [] mat;
}


