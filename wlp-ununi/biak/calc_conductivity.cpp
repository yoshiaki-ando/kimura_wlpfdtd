#include "wlp_fdtd2d.h"

void calc_conductivity(double *x_pos, double *y_pos,
           double **conduct_media,
		       double *conduct_pmlx_ez, double *conduct_pmly_ez,
		       double *conduct_pmlx_h, double *conduct_pmly_h
           ){
  
/*** Media ********************************************/
  //half_conductivity[0][0]はi=0.5, j=0.5の位置の導電率
  double **half_conductivity = allocate_memory_d2d(NX7, NY7, 0.0);
  double **input_conductivity = allocate_memory_d2d(NX, 1000, 0.0);

  std::ifstream ifs("media_biak.dat");
  if(!ifs) {
    std::cout << "Media data does not exist." << std::endl;
    std::exit(0);
  }
  double lat, lon, x, flag, depth;
  for(int i = 0; i < NX; i++){
    for(int j = 0; j < 1000; j++){
      ifs >> lat  >> lon >> x >> depth >> input_conductivity[i][j] >> flag;
    }  
  }
  ifs.close();

  for( int i = NX3; i < NX4; i++ ){
    for( int j = NY3; j < NY4; j++ ){
      int jy = std::round(Y_SURFACE/DY)-1-(j-NY3);
      if(jy >= 0)
        half_conductivity[i][j] = input_conductivity[i-NX3][jy];
    } 
  }
  // 解析領域の下
  for(int i = NX3; i < NX4; i++){
    for(int j = 0; j < NY3; j++){
      half_conductivity[i][j] = half_conductivity[i][NY3];
    }
  }
  // 解析領域の上
  for(int i = NX3; i < NX4; i++){
    for(int j = NY4; j < NY7; j++){
      half_conductivity[i][j] = half_conductivity[i][NY4-1];
    }
  }
  // 解析領域の左
  for(int i = 0; i < NX3; i++){
    for(int j = 0; j < NY7; j++){
      half_conductivity[i][j] = half_conductivity[NX3][j];
    }
  }
  // 解析領域の右
  for(int i = NX4; i < NX7; i++){
    for(int j = 0; j < NY7; j++){
      half_conductivity[i][j] = half_conductivity[NX4-1][j];
    }
  }
  // 周囲４セルの平均
  for( int i = 1; i < NX7; i++ ){
    for( int j = 1; j < NY7; j++ ){
      conduct_media[i][j] = 
      (half_conductivity[i-1][j-1]+half_conductivity[i-1][j]
                               +half_conductivity[i][j-1]+half_conductivity[i][j])/4.0;
    }
  }
  free_memory_d2d(half_conductivity);


/*** Media (up to here) *******************************/


  
/* PML *********************************************/
  const int l = L_PML;
  const double m = M_PML;
  const double r_pml = R_PML;
  const double max_conduct_pmlx = - (m+1)*C0*std::log(r_pml) / 2.0/l/dx_max; // [/s] (= [S/m]/[F/m]) DXをdx_max
  const double max_conduct_pmly = - (m+1)*C0*std::log(r_pml) / 2.0/l/dy_max; // [/s]
  
  for( int i = 0; i <= l; i++ ){
    conduct_pmlx_ez[i] = max_conduct_pmlx * std::pow((double)(l-i)/l, m);
  }
  for( int i = NX7-l; i <= NX7; i++ ){
    conduct_pmlx_ez[i] = max_conduct_pmlx * std::pow((double)(i-(NX7-l))/l, m); //NXをNX7
  }
  for( int j = 0; j <= l; j++ ){
    conduct_pmly_ez[j] = max_conduct_pmly * std::pow((double)(l-j)/l, m);
  }
  for( int j = NY7-l; j <= NY7; j++ ){
    conduct_pmly_ez[j] = max_conduct_pmly * std::pow((double)(j-(NY7-l))/l, m); //NYをNY7
  }

  for( int i = 0; i < l; i++ ){
    double half_i = i + 0.50;
    conduct_pmlx_h[i] = max_conduct_pmlx * std::pow((double)(l-half_i)/l, m);
  }
  for( int i = NX7-l; i < NX7; i++ ){
    double half_i = i + 0.50;
    conduct_pmlx_h[i] = max_conduct_pmlx * std::pow((double)(half_i-(NX7-l))/l, m); //NXをNX7
  }
  for( int j = 0; j < l; j++ ){
    double half_j = j + 0.50;
    conduct_pmly_h[j] = max_conduct_pmly * std::pow((double)(l-half_j)/l, m);
  }
  for( int j = NY7-l; j < NY7; j++ ){
    double half_j = j + 0.50;
    conduct_pmly_h[j] = max_conduct_pmly * std::pow((double)(half_j-(NY7-l))/l, m); //NYをNY7
  }
/*** PML (up to here) *****************************/

  return;
}
