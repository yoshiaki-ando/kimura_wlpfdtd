/***************************************************************
 * 観測点での最大値を計算
 * データ出力(1. 全領域, 2. 観測点)
 * gptファイルを出力(1. 全領域, 2. 観測点Ez, 3. 観測点Bx)
***************************************************************/
#include "wlp_fdtd2d.h"
//#include <experimental/filesystem>
#include <sys/stat.h>

void output_field_1d(double *field, int forN, double delta, std::string filename){
  
  std::ofstream ofs(filename);
  for(int i = 0; i < forN; i++){
    ofs << i*delta << " " << field[i] << std::endl;
  }
  ofs.close();

  return;
}

void output_field_2d(double ***field, int nx, int ny, int nt, double delta, int n_thin_out, std::string field_type, std::string directory_of_data){
  //std::experimental::filesystem::create_directory(directory_of_data + "/2d_data");
  mkdir( (directory_of_data + "/2d_data").c_str(), 0755 );
  for(int n = 0; n < nt+1; n++){
    std::ofstream ofs(directory_of_data + "/2d_data/" + field_type + "_" + std::to_string(n) + ".dat");
    for(int i = 0; i < nx; i+=n_thin_out){
      for(int j = 0; j < ny; j+=n_thin_out){
        ofs << i*DX << " " << j*DY << " " << field[i][j][n] << std::endl;
      }
      ofs << std::endl;
    }
    display_progress("２次元データ出力", n, nt);
    ofs.close();
  }
  
  return;
}
