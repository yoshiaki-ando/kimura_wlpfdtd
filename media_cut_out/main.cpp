#include<iostream>
#include<cmath>
#include<string>
#include<fstream>
#include<vector>
#include<matio.h>

#include <memory_allocate.h>

#include "cut_out.h"

int main(void){

  /* パラメタのチェック */
  if(!(-90 < lat_epi && lat_epi < 90 || -180 < lon_epi && lon_epi < 180
  || -90 < lat_obs && lat_obs < 90 || -180 < lon_obs && lon_obs < 180)){
    std::cout << "[Error] Bad coordinate are given." << std::endl << 
    "(latitude: -90~90[deg], longitude: -180~180[deg])" << std::endl;
    std::exit(0);
  }

  coordinate_sph_orth a(lat_epi_rad, lon_epi_rad);
  coordinate_sph_orth b(lat_obs_rad, lon_obs_rad);
  std::cout << "Epicenter: "; a.print_params();
  std::cout << "Observatory: "; b.print_params();

  double dist_obs_epi = dist(a, b);
  const int nt = std::ceil(Rr/dr);
  std::cout << "obs-epi = " << dist_obs_epi << " [m]" << std::endl;
  if(dist_obs_epi > Rr) {
    std::cout << "[Error] Too small Rr" << std::endl;
    std::exit(0);
  }
  std::cout <<"arg(a, b) = "<< arg(a, b) / M_PI *180.0 << "[°]"<< std::endl;

  double *t = new double [nt+1]; // r(t) = r_s + t(r_e - r_s) の変数t
  divide_param(a, b, t, nt); // 球面上で等間隔になるようなtの集合を計算
  std::vector<coordinate_sph_orth> vt(nt+1); // 球面上の直線上を等間隔に指すベクトル
  for(int i = 0; i < nt+1; i++){
    vt[i].vec_line_on_sph(a, b, t[i]); // 球面上を等間隔に指すベクトル(直交座標)を計算
    vt[i].sph_coordinate(); // 直交座標から球座標(lat,lon)を計算
    if(i % (nt/20) == 0 || i == nt){
      std::cout << i << "/" << nt << " " << t[i] << " ";
      vt[i].print_params();
    }
  }
  delete [] t;

  // 等間隔にとれてるかチェック
  for(int i = 0; i < nt; i+=20){
    if(i+40 <= nt){
      if( std::abs(dist(vt[i], vt[i+20])-dist(vt[i+20], vt[i+40])) > 1.0 ) { 
        std::cout << "[Error] Not equal dist." << std::endl; 
        std::cout << int(dist(vt[10], vt[30])) << " " << int(dist(vt[30], vt[50])) << std::endl;
        std::exit(0); 
      }
    }
  }

  mat_t *matfp;
  matvar_t *matvar;
  matfp = Mat_Open("input_data/MODEL-VER-5.mat", MAT_ACC_RDONLY);
  matvar = Mat_VarRead(matfp,"MODEL");
  /*
  printf("0: %d\n", int(matvar->dims[0]));
  printf("1: %d\n", int(matvar->dims[1]));
  printf("2: %d\n", int(matvar->dims[2]));
  */
  const int nlat = int(matvar->dims[0]);
  const int nlon = int(matvar->dims[1]);
  const int nz = int(matvar->dims[2]);
  const double dlat = 180./nlat;
  const double dlon = 360./nlon;
  const double dz = 100.e3/nz;

  double *pd4 = (double*)matvar->data;
  
/*
  std::ofstream ofs_test("surface_media.dat");
  for(int k = 0; k < nlat; k++){
    for(int j = 0; j < nlon; j++){
      for(int i = 0; i < nz; i+=10){
         if(i==0) ofs_test << j*dlon-179.875 << " " << k*dlat-89.875 << " " << pd4[k+j*720+i*1440*720] << std::endl;
      }
    }
    ofs_test << std::endl;
  }
  ofs_test.close();
*/

  std::vector<coordinate_sph_orth> around_point(4);
  std::vector<std::pair<int, int>> p(4);
  p = {{0, 0}, {1, 0}, {0, 1}, {1, 1}}; // (lat, lon)のインクリメント
  double **line_media = allocate_memory2d(nt+1, nz+1, 0.0);

  // 線形補完
  for(int i = 0; i < nt+1; i++) {
    const int lat_i = int((rad_to_deg(vt[i].lat)+89.875)/dlat); // -89.875degを0番目とする
    const int lon_i = int((rad_to_deg(vt[i].lon)+179.875)/dlon);  // -179.875degを0番目とする
    if(i % (nt/10) == 0) {
      std::cout << "(dlat, dlon) = (" << dlat << ", " << dlon << "), (lat_i, lon_i) = (" << lat_i << ", " << lon_i << ")" << std::endl;
      vt[i].print_params();
    }

    for(int j = 0; j < 4; j++){
      around_point[j].lat = deg_to_rad((lat_i+p[j].first)*dlat) - deg_to_rad(89.875);
      around_point[j].lon = deg_to_rad((lon_i+p[j].second)*dlon) - deg_to_rad(179.875);
      around_point[j].orth_coordinate();
    }

    for(int k = 0; k < nz; k++){
      std::vector<double> coef(4);
      std::vector<double> log_sig(3);
      if(around_point[1].lat -1.0*(vt[i].lon - around_point[1].lon) > vt[i].lat) { // 左下
         log_sig = {
           std::log10(pd4[ (lat_i) + nlat*(lon_i) + nlat*nlon*k ]),
           std::log10(pd4[ (lat_i+p[1].first) + nlat*(lon_i+p[1].second) + nlat*nlon*k ]),
           std::log10(pd4[ (lat_i+p[2].first) + nlat*(lon_i+p[2].second) + nlat*nlon*k ])
         };
         coef = calc_coef_plane(around_point[0], around_point[1], around_point[2], log_sig);
      } else { // 右上
         log_sig = {
           std::log10(pd4[ (lat_i+p[1].first) + nlat*(lon_i+p[1].second) + nlat*nlon*k ]),
           std::log10(pd4[ (lat_i+p[3].first) + nlat*(lon_i+p[3].second) + nlat*nlon*k ]),
           std::log10(pd4[ (lat_i+p[2].first) + nlat*(lon_i+p[2].second) + nlat*nlon*k ])
         };
         coef = calc_coef_plane(around_point[1], around_point[3], around_point[2], log_sig);
      }
      line_media[i][k] = std::pow(10,linear_interpolaton(coef, vt[i]));
    } 
  }
  Mat_Close(matfp);

  const int t_epi = std::round( (Rr-dist_obs_epi)/2./dr );
  const int t_obs = t_epi + std::round( dist_obs_epi/dr );
  std::ofstream ofs("output_data/media_" + eq_name +".dat");
  for(int i = 0; i <= nt; i++){
    for(int j = 0; j < nz; j++){
      ofs << rad_to_deg(vt[i].lat) << " " << rad_to_deg(vt[i].lon) << " " << i*dr << " " << -j*dz-50. << " " << line_media[i][j] << " " << (i==t_epi || i==t_obs) <<std::endl;
    }
    ofs << std::endl;
  }
  ofs.close();

  return 0;
}
